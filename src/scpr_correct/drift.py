"""The SCPR gauge drift curves, evaluated from the published drift model.

Each station's drift model is fitted by the COSZO SCPR pipeline
(https://github.com/coszo-hub/self-calibrating-pressure) from every in-situ calibration of
that station and published daily as
`SCPR-data-collection/calibration/<NET>.<STA>/drift_model.json` (default calibration
averaging) and `drift_model_masked.json` (masked averaging). The file names its station and
which location holds which gauge. This module reads it and evaluates, for each gauge,
exactly what the pipeline subtracts:

  drift(t) = curve(t) - reference
    curve:     a + b*d + c*exp(-d/tau)   (d = days since the first calibration)
               after the last calibration: a weighted line through the last events,
               tapered in from the curve over the last calibration interval
    before the first calibration: held at its value there
    after valid_until (drift uncertainty > 100 Pa): undefined (NaN) - not corrected
  sigma(t) = sqrt(sigma0^2 + (rate * years since the last calibration)^2)
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from importlib import resources
from pathlib import Path
import urllib.request

import numpy as np
from obspy import UTCDateTime

# Where each station's calibration record is published (the pipeline repository).
CALIBRATION_URL = ("https://raw.githubusercontent.com/coszo-hub/self-calibrating-pressure/main/"
                   "SCPR-data-collection/calibration/{station_id}/{name}")
# AXCC1 was first published from the Axial Seamount repository; used if the above is unreachable.
LEGACY_URLS = {"OO.AXCC1": "https://raw.githubusercontent.com/Axial-Seamount/self-calibrating-pressure/main/"
                           "SCPR-data-collection/calibration/{name}"}
MODEL_FILES = {"glenn": "drift_model.json", "masked": "drift_model_masked.json"}
DEFAULT_NETWORK = "OO"


def station_id(station: str) -> str:
    """'AXCC1' or 'OO.AXCC1' -> 'OO.AXCC1'."""
    return station if "." in station else f"{DEFAULT_NETWORK}.{station}"


def _epoch(text: str) -> float:
    return UTCDateTime(text).timestamp


@dataclass
class GaugeDrift:
    gauge: str
    kind: str
    coef: tuple[float, float, float] | None
    tau_days: float | None
    t_first: float
    t_last: float
    t_valid: float
    ref: float
    tail_level: float
    tail_slope_per_day: float
    taper_days: float
    sigma0: float
    sigma_rate_per_yr: float
    event_t: np.ndarray
    event_y: np.ndarray

    @classmethod
    def from_report(cls, g: dict) -> "GaugeDrift":
        ref = float(g["reference"]["relative_pa"])
        events = g["events"]
        coef = g.get("coefficients_pa_per_day")
        return cls(
            gauge=g["gauge"], kind=g["model"],
            coef=None if coef is None else tuple(float(c) for c in coef),
            tau_days=None if g.get("tau_days") is None else float(g["tau_days"]),
            t_first=_epoch(events[0]["t0_utc"]), t_last=_epoch(g["last_event_utc"]),
            t_valid=_epoch(g["valid_until_utc"]), ref=ref,
            tail_level=float(g["tail"]["level_at_last_event_pa"]) + ref,
            tail_slope_per_day=float(g["tail"]["slope_pa_per_year"]) / 365.25,
            taper_days=float(g["tail"]["taper_days"]),
            sigma0=float(g["uncertainty"]["sigma0_pa"]), sigma_rate_per_yr=float(g["uncertainty"]["rate_pa_per_year"]),
            event_t=np.array([_epoch(e["t0_utc"]) for e in events]),
            event_y=np.array([float(e["drift_pa"]) + ref for e in events]),
        )

    def _interior(self, t: np.ndarray) -> np.ndarray:
        if self.coef is None:                                   # linear interpolation model
            return np.interp(t, self.event_t, self.event_y)
        a, b, c = self.coef
        d = (t - self.t_first) / 86400.0
        return a + b * d + (c * np.exp(-d / self.tau_days) if self.tau_days else 0.0)

    def _curve(self, t: np.ndarray) -> np.ndarray:
        out = self._interior(t)
        late = t > self.t_last
        if late.any():
            dt = (t[late] - self.t_last) / 86400.0
            a, b = self.tail_level, self.tail_slope_per_day
            step = float(self._interior(np.array([self.t_last]))[0]) - a
            taper = np.clip(1.0 - dt / self.taper_days, 0.0, 1.0) if self.taper_days > 0 else 0.0
            out[late] = a + b * dt + step * taper
        return out

    def delta(self, t) -> np.ndarray:
        """Drift (Pa) at epoch seconds t; NaN after valid_until (not corrected)."""
        t = np.asarray(t, dtype=float)
        out = self._curve(np.maximum(t, self.t_first)) - self.ref
        out[t > self.t_valid] = np.nan
        return out

    def sigma(self, t) -> np.ndarray:
        """1-sigma uncertainty (Pa) of the drift at t."""
        dt_yr = np.maximum(np.asarray(t, dtype=float) - self.t_last, 0.0) / 86400.0 / 365.25
        return np.sqrt(self.sigma0**2 + (self.sigma_rate_per_yr * dt_yr) ** 2)


@dataclass
class DriftModel:
    network: str
    station: str
    gauges: dict[str, GaugeDrift]
    gauge_location: dict[str, str]          # gauge -> location code
    flag_location: str                      # calibration flag (LEF) location
    table_version: str
    offsets: str
    source: str
    mask_before_s: float
    mask_after_s: float
    threshold_pa: float | None              # None = auto (median - auto_drop_pa)
    auto_drop_pa: float

    @property
    def id(self) -> str:
        return f"{self.network}.{self.station}"

    @property
    def location_gauge(self) -> dict[str, str]:
        return {loc: g for g, loc in self.gauge_location.items()}

    @classmethod
    def from_dict(cls, report: dict, source: str) -> "DriftModel":
        mask = report.get("mask", {})
        # The station block says what the file applies to; files from before it existed are AXCC1's.
        st = report.get("station") or {"network": "OO", "station": "AXCC1",
                                       "gauges": [{"gauge": "PARO1", "location": "40"},
                                                  {"gauge": "PARO2", "location": "41"}],
                                       "calibration_flag": {"location": "40", "channel": "LEF"}}
        thr = mask.get("threshold_pa", 15400e3)
        return cls(network=st["network"], station=st["station"],
                   gauges={g["gauge"]: GaugeDrift.from_report(g) for g in report["gauges"]},
                   gauge_location={g["gauge"]: str(g["location"]) for g in st["gauges"]},
                   flag_location=str(st.get("calibration_flag", {}).get("location", "40")),
                   table_version=report.get("table_version", "?"), offsets=report.get("offsets", "glenn"),
                   source=source, mask_before_s=float(mask.get("before_s", 60.0)),
                   mask_after_s=float(mask.get("after_s", 1800.0)),
                   threshold_pa=None if thr in (None, "auto") else float(thr),
                   auto_drop_pa=float(mask.get("auto_drop_pa", 100e3)))

    @classmethod
    def load(cls, station: str = "AXCC1", source: str = "latest", offsets: str = "glenn",
             timeout: float = 20.0) -> "DriftModel":
        """The calibration record of `station` ("AXCC1" or "OO.AXCC1").

        source: "latest" (from the pipeline repository; falls back to the copy bundled with this
        package when offline), "bundled", or a path / URL to a drift_model JSON file.
        """
        sid, name = station_id(station), MODEL_FILES[offsets]
        if source not in ("latest", "bundled"):
            if source.startswith(("http://", "https://")):
                with urllib.request.urlopen(source, timeout=timeout) as r:
                    return cls.from_dict(json.loads(r.read().decode()), source)
            return cls.from_dict(json.loads(Path(source).read_text()), str(source))
        if source == "latest":
            for template in (CALIBRATION_URL, LEGACY_URLS.get(sid)):
                if template is None:
                    continue
                url = template.format(station_id=sid, name=name)
                try:
                    with urllib.request.urlopen(url, timeout=timeout) as r:
                        return cls.from_dict(json.loads(r.read().decode()), url)
                except OSError:
                    continue                                     # not published there / offline
        bundled = resources.files("scpr_correct").joinpath("data", sid, name)
        if not bundled.is_file():
            raise LookupError(f"no calibration record for {sid}: none published yet at "
                              f"{CALIBRATION_URL.format(station_id=sid, name=name)} and none bundled. "
                              f"A station's record starts with its first calibration.")
        return cls.from_dict(json.loads(bundled.read_text()), f"bundled {sid}/{name}")

    def coverage(self) -> tuple[float, float, float]:
        """(first calibration, last calibration, corrected-until) epochs (first gauge)."""
        g = next(iter(self.gauges.values()))
        return g.t_first, g.t_last, g.t_valid

    def describe(self) -> str:
        first, last, valid = self.coverage()
        return (f"{self.id} calibration record {self.source} (table {self.table_version}, averaging "
                f"{self.offsets}; calibrations {UTCDateTime(first).isoformat()[:10]} .. "
                f"{UTCDateTime(last).isoformat()[:10]}, corrected until {UTCDateTime(valid).isoformat()[:10]})")
