"""Drift-correct the gauges of COSZO Self-Calibrating Pressure Recorders (SCPR), e.g. OO.AXCC1."""

from .correct import correct_stream, mask_intervals, write_sds
from .drift import DriftModel, GaugeDrift

__all__ = ["DriftModel", "GaugeDrift", "correct_stream", "mask_intervals", "write_sds"]
__version__ = "0.2.0"
