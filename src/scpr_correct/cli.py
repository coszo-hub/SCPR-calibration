"""scpr-correct: SCPR gauge pressure (MiniSEED) -> drift-corrected MiniSEED, for any COSZO SCPR station.

  scpr-correct --station AXCC1 --day 2024-08-28                 # download from EarthScope, correct
  scpr-correct --station AXCC1 --day 2024-08-27 --end 2024-08-31
  scpr-correct /path/to/mseed                                   # local files/folders (station(s) from the data)
  scpr-correct --station AXCC1 --info                           # the calibration record that would be used
  ... --offsets masked                                          # the second calibration averaging
  ... --calibration my_drift_model.json                         # a specific calibration record

Each station's calibration record comes from the COSZO SCPR pipeline
(github.com/coszo-hub/self-calibrating-pressure, calibration/<NET>.<STA>/), updated daily.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import date, timedelta
import os
from pathlib import Path
import sys

import numpy as np
from obspy import Stream, UTCDateTime, read

from .correct import correct_stream, write_sds
from .drift import DriftModel, station_id


def read_inputs(paths: list[Path], start: date | None, end: date | None) -> Stream:
    """Every LDO and LEF trace in the given files/folders (all stations)."""
    files = []
    for p in paths:
        if p.is_dir():
            for root, _, names in os.walk(p):
                files += [Path(root) / n for n in names if not n.startswith(".") and not n.endswith(".part")]
        else:
            files.append(p)
    st = Stream()
    lo = UTCDateTime((start - timedelta(days=1)).isoformat()) if start else None
    hi = UTCDateTime((end + timedelta(days=2)).isoformat()) if end else None
    for f in sorted(files):
        if ".LDO." not in f.name and ".LEF." not in f.name and f.suffix.lower() not in (".mseed", ".ms", ".miniseed"):
            continue                                       # skip obviously unrelated files quickly
        try:
            tr_all = read(str(f), format="MSEED", starttime=lo, endtime=hi)
        except Exception:
            continue
        st += tr_all.select(channel="LDO") + tr_all.select(channel="LEF")
    return st


def fetch_fdsn(model: DriftModel, start: date, end: date, base_url: str) -> Stream:
    """The station's gauge pressure and calibration flag for [start, end] plus a day either side
    (so calibrations crossing midnight are removed in full)."""
    from obspy.clients.fdsn import Client
    from obspy.clients.fdsn.header import FDSNNoDataException
    client = Client(base_url)
    lo, hi = UTCDateTime((start - timedelta(days=1)).isoformat()), UTCDateTime((end + timedelta(days=2)).isoformat())
    st = Stream()
    wanted = [(loc, "LDO") for loc in sorted(model.location_gauge)] + [(model.flag_location, "LEF")]
    for loc, cha in wanted:
        try:
            st += client.get_waveforms(model.network, model.station, loc, cha, lo, hi)
        except FDSNNoDataException:
            print(f"  {model.id}.{loc}.{cha}: no data at {base_url} for {lo.date} .. {(hi - 1).date}")
    return st


def day_report(model: DriftModel, corrected: Stream, raw: Stream, day: date) -> list[str]:
    """What happened to each gauge on one day, in words."""
    lo = UTCDateTime(day.isoformat()).timestamp
    hi = lo + 86400
    first, last, valid = model.coverage()
    lines = []
    for loc, gauge in sorted(model.location_gauge.items()):
        n_raw = sum(int(((t >= lo) & (t < hi)).sum()) for t in
                    (tr.stats.starttime.timestamp + np.arange(tr.stats.npts) / tr.stats.sampling_rate
                     for tr in raw.select(network=model.network, station=model.station, location=loc, channel="LDO")))
        n_cor = sum(int(((t >= lo) & (t < hi)).sum()) for t in
                    (tr.stats.starttime.timestamp + np.arange(tr.stats.npts) / tr.stats.sampling_rate
                     for tr in corrected.select(location=loc)))
        head = f"  {day} {model.id}.{loc} {gauge}:"
        if n_raw == 0:
            lines.append(f"{head} NO DATA (nothing in the input / at the data centre for this day)")
            continue
        notes = []
        if n_raw < 86400:
            notes.append(f"{86400 - n_raw} s missing in the raw data (gaps)")
        if lo >= valid:
            notes.append(f"NOT CORRECTED: after {UTCDateTime(valid).date}, too long since the last calibration "
                         f"({UTCDateTime(last).date})")
        elif hi > valid:
            notes.append(f"corrected only until {UTCDateTime(valid).isoformat()[:19]}")
        elif hi > last:
            notes.append(f"PROVISIONAL: after the last calibration ({UTCDateTime(last).date}), drift extrapolated")
        if lo < first:
            notes.append(f"before the first calibration ({UTCDateTime(first).date}): drift held at its first value")
        masked = n_raw - n_cor - (0 if lo < valid else n_raw)
        if lo < valid and masked > 0:
            notes.append(f"{masked} s removed around a calibration")
        lines.append(f"{head} {n_cor} of {n_raw} samples corrected" + ("; " + "; ".join(notes) if notes else ""))
    return lines


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="scpr-correct", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="*", type=Path, help="local MiniSEED files/folders (else: download from EarthScope)")
    ap.add_argument("--station", help="station code, e.g. AXCC1 or OO.AXCC1 (required to download)")
    ap.add_argument("--day", "--start", dest="start", type=date.fromisoformat, help="(first) UTC day")
    ap.add_argument("--end", type=date.fromisoformat, help="last UTC day (default: --day)")
    ap.add_argument("--out", type=Path, default=Path("scpr_corrected"), help="output SDS root (default ./scpr_corrected)")
    ap.add_argument("--calibration", "--model", dest="calibration", default="latest",
                    help="calibration record: latest (default), bundled, or a path/URL to a drift_model*.json")
    ap.add_argument("--offsets", choices=("glenn", "masked"), default="glenn",
                    help="calibration averaging: glenn (default) or masked")
    ap.add_argument("--settle-seconds", type=float, help="removed after a calibration (default: the record's, 1800 s)")
    ap.add_argument("--fdsn-url", default="IRIS", help="FDSN data centre (default: IRIS = EarthScope)")
    ap.add_argument("--fdsn", action="store_true", help=argparse.SUPPRESS)      # old flag; downloading is the default
    ap.add_argument("--info", action="store_true", help="print the calibration record that would be used and exit")
    args = ap.parse_args(argv)
    end = args.end or args.start

    if args.info:
        if not args.station:
            ap.error("--info needs --station")
        model = DriftModel.load(args.station, args.calibration, args.offsets)
        print(model.describe())
        now = UTCDateTime().timestamp
        for g, d in model.gauges.items():
            drift = float(d.delta([now])[0])
            state = (f"drift now {drift:.0f} Pa (1 sigma {float(d.sigma([now])[0]):.0f} Pa)"
                     + (" PROVISIONAL (after the last calibration)" if now > d.t_last else "") if drift == drift
                     else f"not corrected now: too long since the last calibration (after {UTCDateTime(d.t_valid).date})")
            print(f"  {g} (loc {model.gauge_location[g]}): {d.kind}; {state}")
        return 0

    if args.inputs:
        st = read_inputs(args.inputs, args.start, end)
        stations = sorted({f"{tr.stats.network}.{tr.stats.station}" for tr in st.select(channel="LDO")})
        if args.station:
            stations = [s for s in stations if s == station_id(args.station)]
        if not stations:
            print("no LDO data" + (f" for {station_id(args.station)}" if args.station else "") + " in the input")
            return 1
    else:
        if not (args.station and args.start):
            ap.error("give local MiniSEED files, or --station and --day to download from EarthScope")
        stations, st = [station_id(args.station)], None

    status = 0
    for sid in stations:
        try:
            model = DriftModel.load(sid, args.calibration, args.offsets)
        except LookupError as exc:
            print(f"{sid}: {exc}")
            status = 1
            continue
        if model.id != sid:
            print(f"{sid}: the calibration record given is for {model.id}; skipped")
            status = 1
            continue
        print(model.describe())
        data = fetch_fdsn(model, args.start, end, args.fdsn_url) if st is None else \
            st.select(network=model.network, station=model.station)
        if not data.select(channel="LDO"):
            print(f"  no {sid} LDO data for {args.start} .. {end}")
            status = 1
            continue
        if not data.select(channel="LEF"):
            print(f"  note: no {sid}.{model.flag_location}.LEF calibration flag; calibrations found from the pressure alone")
        corrected, _ = correct_stream(data, model, args.settle_seconds)
        if args.start:
            days = [args.start + timedelta(days=i) for i in range((end - args.start).days + 1)]
        else:
            t = [tr.stats.starttime for tr in data.select(channel="LDO")] + [tr.stats.endtime for tr in data.select(channel="LDO")]
            first, last = min(t).date, max(t).date
            days = [first + timedelta(days=i) for i in range((last - first).days + 1)]
        written = write_sds(corrected, args.out, set(days))
        for d in days:
            for line in day_report(model, corrected, data, d):
                print(line)
        print(f"  {len(written)} file(s) written under {args.out}")
    return status


if __name__ == "__main__":
    sys.exit(main())
