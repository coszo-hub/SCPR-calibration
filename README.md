# SCPR calibration

Drift correction for COSZO **Self-Calibrating Pressure Recorders (SCPR)**: download a station's
gauge pressure from EarthScope (or use your own MiniSEED) and get calibrated, drift-corrected
MiniSEED.

Each SCPR has two Paroscientific gauges (channel `LDO`, Pa; location `40` = PARO1, `41` = PARO2)
that drift slowly. Every few months the instrument switches them to an internal dead-weight
tester; those in-situ calibrations give each gauge's drift curve, which this package subtracts:

```
corrected LDO (quality Q) = raw LDO - drift(t)
```

The calibration records are produced and updated daily, per station, by the COSZO SCPR
pipeline, [coszo-hub/self-calibrating-pressure](https://github.com/coszo-hub/self-calibrating-pressure)
(`SCPR-data-collection/calibration/<NET>.<STA>/`). This package only applies them.

| Station | Site | Calibration record |
|---|---|---|
| OO.AXCC1 | Axial Seamount central caldera | since 2018-07-09 |
| OO.CZMID | Mid Slope | starts with its first calibration (deployment planned summer 2027) |
| OO.CZOFF | Oregon Offshore | starts with its first calibration (deployment planned summer 2027) |

## Use it from a clone

```bash
git clone https://github.com/coszo-hub/SCPR-calibration && cd SCPR-calibration
python scpr-correct.py --station AXCC1 --day 2024-08-28
```

Needs Python ≥ 3.9 with NumPy and ObsPy (`conda install -c conda-forge obspy`). This

1. downloads `OO.AXCC1` `LDO` (loc 40, 41) and the calibration flag `40.LEF` for the day **and
   the day before and after** from EarthScope (so a calibration crossing midnight is removed in
   full);
2. loads the station's latest calibration record;
3. writes `scpr_corrected/2024/OO/AXCC1/LDO.D/OO.AXCC1.40.LDO.D.2024.241` (and `.41.`): one
   MiniSEED file per gauge per day, quality `Q`, Pa;
4. reports, per gauge and day, what it did — and what is missing:

```
OO.AXCC1 calibration record … (calibrations 2018-07-09 .. 2026-02-05, corrected until 2027-02-05)
  2020-05-04 OO.AXCC1.40 PARO1: 83197 of 86400 samples corrected; 3203 s removed around a calibration
  2026-05-10 OO.AXCC1.40 PARO1: 86400 of 86400 samples corrected; PROVISIONAL: after the last calibration …
  2027-03-01 OO.AXCC1.40 PARO1: 0 of 86400 samples corrected; NOT CORRECTED: more than a year after the last calibration …
  2024-08-28 OO.AXCC1.41 PARO2: NO DATA (nothing in the input / at the data centre for this day)
```

More:

```bash
python scpr-correct.py --station AXCC1 --day 2024-08-27 --end 2024-08-31 --out corrected/
python scpr-correct.py /path/to/mseed                  # your own files/folders; stations taken from the data
python scpr-correct.py --station AXCC1 --info          # which calibration record would be used
python scpr-correct.py ... --offsets masked            # the second calibration averaging
python scpr-correct.py ... --calibration drift_model.json   # a specific (e.g. archived) calibration record
```

Or install it (`pip install git+https://github.com/coszo-hub/SCPR-calibration`) and run
`scpr-correct` with the same options.

## In Python

```python
from obspy import read
from scpr_correct import DriftModel, correct_stream

model = DriftModel.load("AXCC1")              # latest calibration record of OO.AXCC1
st = read("raw/*.LDO.*") + read("raw/*.LEF.*")
corrected, stats = correct_stream(st, model)
drift = model.gauges["PARO1"].delta(times)    # drift (Pa) at epoch seconds; NaN = not corrected
sigma = model.gauges["PARO1"].sigma(times)    # its 1-sigma uncertainty (Pa)
```

## What it does

* **Calibration record.** `--calibration latest` (default) downloads
  `calibration/<NET>.<STA>/drift_model.json` (or `drift_model_masked.json`) from the pipeline
  repository; offline it falls back to the copy bundled with this package. The file names its
  station and which location holds which gauge, so traces of other stations are never touched.
  Every run prints the calibration-table version it used.
* **Drift curve** per gauge: `a + b·t + c·exp(−t/τ)` through the calibrations, zero at the first
  calibration. After the last calibration it follows a line through the last four, so recent
  data are **provisional** until the next calibration, with the drift uncertainty
  `σ = √(40² + (170 Pa/yr · years since the last calibration)²)` (~138 Pa after 9 months). More
  than a year after the last calibration the data are **not corrected** (left out, and reported). Before the first calibration the drift is
  held at its first value.
* **Calibrations removed:** from 60 s before the calibration flag (`LEF` = 1) or the pressure
  dropping below the station's calibration threshold, to 30 min after (`--settle-seconds`).
* **Calibration averaging:** `glenn` (default; Sasagawa et al. 2016) or `masked` (spin-up skipped,
  logged re-spins and valve operations removed). They differ mostly by a constant.
* **Timing.** Data at EarthScope is already timed by the pipeline, including the periods it
  re-timed when the logger wrote no timestamps; gaps are reported, never filled.

The output is identical, file for file, to the pipeline's own corrected MiniSEED (checked on a
calibration day and a normal day, both averagings), and the drift curves agree to 1e-12 Pa.

## Tests

```bash
pip install pytest && pytest
```

## License

MIT — COSZO: Cascadia Subduction Zone Observatory.
