# SCPR calibration

Drift correction for the **Self-Calibrating Pressure Recorder (SCPR)** at Axial Seamount
central caldera (OO.AXCC1): raw gauge-pressure MiniSEED in, calibrated (drift-corrected)
MiniSEED out.

The SCPR's two Paroscientific gauges (PARO1 = location `40`, PARO2 = location `41`, channel
`LDO`, Pa) drift slowly. Every few months the instrument switches them to an internal
dead-weight tester; those in-situ calibrations give each gauge's drift curve. This package
subtracts it:

```
corrected LDO (quality Q) = raw LDO - drift(t)
```

## Install

```bash
pip install git+https://github.com/coszo-hub/SCPR-calibration
```

Needs Python ≥ 3.9, NumPy and ObsPy.

## Use

```bash
scpr-correct /path/to/mseed --out corrected/                      # files or folders (e.g. an SDS archive)
scpr-correct /path/to/mseed --start 2024-08-27 --end 2024-08-29 --out corrected/
scpr-correct --fdsn --start 2024-08-27 --end 2024-08-29 --out corrected/   # raw data from EarthScope
scpr-correct /path/to/mseed --offsets masked --out corrected_masked/       # second calibration averaging
scpr-correct --info                                               # which drift model is used
```

Input: OO.AXCC1 `LDO` (loc 40 and/or 41) and, for removing calibrations, the calibration flag
`40.LEF`. Output: one MiniSEED file per gauge per UTC day in SDS layout
(`YEAR/OO/AXCC1/LDO.D/OO.AXCC1.LOC.LDO.D.YEAR.DOY`), data quality `Q`, FLOAT64 values in Pa.

In Python:

```python
from obspy import read
from scpr_correct import DriftModel, correct_stream

model = DriftModel.load()                     # latest drift model (see below)
st = read("raw/*.LDO.*") + read("raw/*.LEF.*")
corrected, stats = correct_stream(st, model)
drift = model.gauges["PARO1"].delta(times)    # drift (Pa) at epoch seconds; NaN = not corrected
sigma = model.gauges["PARO1"].sigma(times)    # its 1-sigma uncertainty
```

## What it does

* **Drift model.** Read from the calibration record that the processing pipeline
  ([Axial-Seamount/self-calibrating-pressure](https://github.com/Axial-Seamount/self-calibrating-pressure))
  fits from every calibration since 2018 and publishes daily
  (`calibration/drift_model.json`, and `drift_model_masked.json` for the masked averaging).
  `--model latest` (default) downloads the current one, falling back to the copy bundled
  with the package when offline; `--model bundled` or a file path/URL pins a version. The run
  prints the calibration-table version used.
* **Drift curve** per gauge: `a + b·t + c·exp(−t/τ)` through the calibrations (weighted, with the
  dead-weight-tester thermal correction), zero at the first calibration (2018-07-09). After the
  last calibration it follows a line through the last four, so recent data are **provisional**
  until the next calibration; once the drift uncertainty
  `σ = √(40² + (170 Pa/yr · years since the last calibration)²)` exceeds 100 Pa (~200 days) the
  data are **not corrected** (left out).
* **Calibrations removed:** from 60 s before the calibration flag (`40.LEF` = 1) or the pressure
  dropping below 15 400 kPa, to 30 min after (`--settle-seconds`). Include the day before your
  period (or `40.LEF`) so a calibration crossing midnight is fully removed.
* **Calibration averaging:** `glenn` (default; Sasagawa et al. 2016: 12 × 100-s bins from the
  pressure minimum, 1σ-trimmed means) or `masked` (spin-up skipped, logged re-spins and valve
  operations removed, bin medians). They differ mostly by a constant.

The output is identical, sample for sample, to the pipeline's own corrected MiniSEED (checked
on a calibration day and a normal day, both averagings), and the drift curves agree to 1e-12 Pa
from 2018 to the end of the record.

## Tests

```bash
pip install pytest && pytest
```

## License

MIT — COSZO: Cascadia Subduction Zone Observatory.
