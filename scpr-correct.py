#!/usr/bin/env python3
"""Run scpr-correct straight from a clone of this repository (no installation needed):

    git clone https://github.com/coszo-hub/SCPR-calibration && cd SCPR-calibration
    python scpr-correct.py --station AXCC1 --day 2024-08-28

Needs Python >= 3.9 with NumPy and ObsPy (e.g. `conda install -c conda-forge obspy`).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from scpr_correct.cli import main  # noqa: E402

sys.exit(main())
