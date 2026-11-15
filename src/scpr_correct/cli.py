"""scpr-correct: raw SCPR MiniSEED (OO.AXCC1) -> drift-corrected MiniSEED.

  scpr-correct /path/to/mseed --out corrected/                 # files or folders (e.g. an SDS archive)
  scpr-correct --fdsn --start 2024-08-27 --end 2024-08-29 --out corrected/   # from EarthScope
  scpr-correct ... --offsets masked                            # the second calibration averaging
  scpr-correct --info                                          # which drift model would be used
"""

from __future__ import annotations

import argparse
from datetime import date, timedelta
import os
from pathlib import Path
import sys

from obspy import Stream, UTCDateTime, read

from .correct import NETWORK, STATION, correct_stream, write_sds
from .drift import DriftModel


def read_inputs(paths: list[Path], start: date | None, end: date | None) -> Stream:
    """Every OO.AXCC1 LDO (loc 40/41) and 40.LEF trace in the given files/folders."""
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
        name = f.name
        if ".LDO." not in name and ".LEF." not in name and f.suffix.lower() not in (".mseed", ".ms", ".miniseed"):
            continue                                       # skip obviously unrelated files quickly
        try:
            tr_all = read(str(f), format="MSEED", starttime=lo, endtime=hi)
        except Exception:
            continue
        st += tr_all.select(network=NETWORK, station=STATION, channel="LDO").select(location="4[01]")
        st += tr_all.select(network=NETWORK, station=STATION, location="40", channel="LEF")
    return st


def fetch_fdsn(start: date, end: date, base_url: str) -> Stream:
    from obspy.clients.fdsn import Client
    client = Client(base_url)
    lo, hi = UTCDateTime((start - timedelta(days=1)).isoformat()), UTCDateTime((end + timedelta(days=2)).isoformat())
    st = Stream()
    for loc, cha in (("40", "LDO"), ("41", "LDO"), ("40", "LEF")):
        try:
            st += client.get_waveforms(NETWORK, STATION, loc, cha, lo, hi)
        except Exception as exc:                            # no data for that channel/period
            print(f"  {NETWORK}.{STATION}.{loc}.{cha}: {exc.__class__.__name__}: {str(exc)[:80]}")
    return st


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="scpr-correct", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="*", type=Path, help="MiniSEED files or folders with OO.AXCC1 LDO (+ 40.LEF)")
    ap.add_argument("--out", type=Path, default=Path("scpr_corrected"), help="output SDS root")
    ap.add_argument("--start", type=date.fromisoformat, help="first UTC day to write")
    ap.add_argument("--end", type=date.fromisoformat, help="last UTC day to write")
    ap.add_argument("--fdsn", action="store_true", help="download the raw data from an FDSN service")
    ap.add_argument("--fdsn-url", default="IRIS", help="FDSN service (default: IRIS = EarthScope)")
    ap.add_argument("--model", default="latest",
                    help="drift model: latest (GitHub, default), bundled, or a path/URL to drift_model*.json")
    ap.add_argument("--offsets", choices=("glenn", "masked"), default="glenn",
                    help="calibration averaging: glenn (default) or masked")
    ap.add_argument("--settle-seconds", type=float, help="mask after a calibration (default: the model's, 1800 s)")
    ap.add_argument("--info", action="store_true", help="print the drift model that would be used and exit")
    args = ap.parse_args(argv)

    model = DriftModel.load(args.model, args.offsets)
    print(model.describe())
    if args.info:
        now = UTCDateTime().timestamp
        for g, d in model.gauges.items():
            drift = float(d.delta([now])[0])
            state = (f"drift now {drift:.0f} Pa (1 sigma {float(d.sigma([now])[0]):.0f} Pa)" if drift == drift
                     else f"not corrected now: drift uncertainty above the limit since "
                          f"{UTCDateTime(d.t_valid).isoformat()[:10]}")
            print(f"  {g}: {d.kind}; {state}")
        return 0
    if args.fdsn:
        if not (args.start and args.end):
            ap.error("--fdsn needs --start and --end")
        st = fetch_fdsn(args.start, args.end, args.fdsn_url)
    else:
        if not args.inputs:
            ap.error("give MiniSEED files/folders, or --fdsn")
        st = read_inputs(args.inputs, args.start, args.end)
    if not st.select(channel="LDO"):
        print("no OO.AXCC1 LDO (loc 40/41) data found")
        return 1
    if not st.select(channel="LEF"):
        print("note: no 40.LEF calibration flag in the input; calibrations are found from the pressure alone")
    corrected, stats = correct_stream(st, model, args.settle_seconds)
    days = None
    if args.start or args.end:
        first = args.start or date(1970, 1, 1)
        last = args.end or date(2100, 1, 1)
        days = {first + timedelta(days=i) for i in range((min(last, date.today()) - first).days + 1)}
    paths = write_sds(corrected, args.out, days)
    for s in stats:
        note = f", {s['past_valid_until']} after valid_until (not corrected)" if s["past_valid_until"] else ""
        prov = f", {s['provisional_samples']} provisional (after the last calibration)" if s["provisional_samples"] else ""
        print(f"  {s['gauge']} (loc {s['location']}): {s['written']} of {s['samples_in']} samples corrected, "
              f"{s['masked_calibration']} in calibrations{note}{prov}")
    print(f"{len(paths)} file(s) written under {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
