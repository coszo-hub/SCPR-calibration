"""Drift-correct SCPR gauge pressure (LDO) of one station, e.g. OO.AXCC1 loc 40 = PARO1, 41 = PARO2.

corrected = raw - drift(t), written as data quality Q. Calibrations are removed: from
mask_before_s before the calibration flag (40.LEF = 1) or the pressure dropping below the
calibration threshold, to mask_after_s after it ends - the same rule, constants and code as
the processing pipeline. Samples after the drift model's valid_until are not written.
"""

from __future__ import annotations

from datetime import date, datetime, timedelta, timezone

import numpy as np
from obspy import Stream, Trace, UTCDateTime

from .drift import DriftModel


def _times(tr: Trace) -> np.ndarray:
    return tr.stats.starttime.timestamp + np.arange(tr.stats.npts) / tr.stats.sampling_rate


def _concat(st: Stream) -> tuple[np.ndarray, np.ndarray]:
    if not len(st):
        return np.array([]), np.array([])
    return (np.concatenate([_times(tr) for tr in st]),
            np.concatenate([np.asarray(tr.data, dtype=float) for tr in st]))


def mask_intervals(ft: np.ndarray, flag: np.ndarray, pt: np.ndarray, p: np.ndarray,
                   threshold: float, before: float, after: float) -> np.ndarray:
    """(start, end) epochs to drop: LEF = 1 or pressure below threshold, padded."""
    bad = np.sort(np.concatenate([ft[flag > 0], pt[p < threshold]])) if len(p) else np.sort(ft[flag > 0])
    if not len(bad):
        return np.empty((0, 2))
    breaks = np.flatnonzero(np.diff(bad) > 60.0)
    starts = np.r_[bad[0], bad[breaks + 1]] - before
    ends = np.r_[bad[breaks], bad[-1]] + after
    merged = [[starts[0], ends[0]]]
    for a, b in zip(starts[1:], ends[1:]):
        if a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return np.array(merged)


def in_intervals(t: np.ndarray, intervals: np.ndarray) -> np.ndarray:
    if not len(intervals):
        return np.zeros(len(t), dtype=bool)
    i = np.searchsorted(intervals[:, 0], t, side="right") - 1
    return (i >= 0) & (t <= intervals[np.maximum(i, 0), 1])


def correct_stream(st: Stream, model: DriftModel, settle_seconds: float | None = None) -> tuple[Stream, list[dict]]:
    """Corrected LDO traces for the model's station in `st` (its gauge locations, e.g. 40/41).

    `st` should also hold the calibration flag (LEF) and, for calibrations that cross the edges
    of the period, the neighbouring hours; without LEF the pressure threshold alone defines
    the calibrations. Traces of other stations are ignored. Returns (corrected stream,
    per-gauge statistics).
    """
    after = model.mask_after_s if settle_seconds is None else settle_seconds
    net, sta = model.network, model.station
    lef = st.select(network=net, station=sta, location=model.flag_location, channel="LEF")
    ft, flag = _concat(lef)
    out, stats = Stream(), []
    for loc, gauge in sorted(model.location_gauge.items()):
        ldo = st.select(network=net, station=sta, location=loc, channel="LDO")
        if not len(ldo):
            continue
        drift_model = model.gauges[gauge]
        pt, p = _concat(ldo)
        threshold = model.threshold_pa if model.threshold_pa is not None else float(np.median(p)) - model.auto_drop_pa
        masks = mask_intervals(ft, flag, pt, p, threshold, model.mask_before_s, after)
        n_in = n_mask = n_late = n_out = 0
        for tr in sorted(ldo, key=lambda x: x.stats.starttime):
            t = _times(tr)
            drift = drift_model.delta(t)
            masked = in_intervals(t, masks)
            keep = ~masked & np.isfinite(drift)
            n_in += len(t)
            n_mask += int(masked.sum())
            n_late += int((~masked & ~np.isfinite(drift)).sum())
            edges = np.flatnonzero(np.diff(np.r_[0, keep.astype(int), 0]))
            for a, b in zip(edges[::2], edges[1::2]):
                new = tr.copy()
                new.data = np.asarray(tr.data[a:b], dtype=float) - drift[a:b]
                new.stats.starttime = tr.stats.starttime + a / tr.stats.sampling_rate
                new.stats.mseed = {"dataquality": "Q"}
                out.append(new)
                n_out += b - a
        stats.append({"station": model.id, "gauge": gauge, "location": loc, "samples_in": n_in,
                      "masked_calibration": n_mask, "past_valid_until": n_late, "written": n_out,
                      "before_first_calibration": int(sum((_times(tr) < drift_model.t_first).sum() for tr in ldo)),
                      "provisional_samples": int(sum((_times(tr) > drift_model.t_last).sum() for tr in ldo))})
    return out, stats


def split_days(st: Stream) -> dict[tuple[str, str, str, date], Stream]:
    """{(network, station, location, UTC day): traces}, cut at midnight (no sample in two files)."""
    days: dict[tuple[str, date], Stream] = {}
    for tr in st:
        t = _times(tr)
        day_index = np.floor(t / 86400.0).astype(np.int64)
        for d in np.unique(day_index):
            idx = np.flatnonzero(day_index == d)
            piece = tr.copy()
            piece.data = tr.data[idx[0]:idx[-1] + 1]
            piece.stats.starttime = tr.stats.starttime + idx[0] / tr.stats.sampling_rate
            key = (tr.stats.network, tr.stats.station, tr.stats.location,
                   (datetime(1970, 1, 1, tzinfo=timezone.utc) + timedelta(days=int(d))).date())
            days.setdefault(key, Stream()).append(piece)
    return days


def sds_path(root, net: str, sta: str, loc: str, cha: str, day: date):
    from pathlib import Path
    doy = day.timetuple().tm_yday
    return Path(root) / f"{day.year}" / net / sta / f"{cha}.D" / f"{net}.{sta}.{loc}.{cha}.D.{day.year}.{doy:03d}"


def write_sds(st: Stream, root, days: set[date] | None = None, reclen: int = 512) -> list:
    """Write corrected traces as one file per gauge per UTC day (SDS layout, like the pipeline)."""
    written = []
    for (net, sta, loc, day), part in sorted(split_days(st).items()):
        if days is not None and day not in days:
            continue
        path = sds_path(root, net, sta, loc, "LDO", day)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + ".part")
        part.write(str(tmp), format="MSEED", reclen=reclen)
        tmp.replace(path)
        written.append(path)
    return written


def day_window(day: date) -> tuple[UTCDateTime, UTCDateTime]:
    lo = UTCDateTime(day.isoformat())
    return lo, lo + 86400
