"""Drift-correct the Axial Seamount Self-Calibrating Pressure Recorder (OO.AXCC1) gauges."""

from .correct import correct_stream, mask_intervals, write_sds
from .drift import DriftModel, GaugeDrift

__all__ = ["DriftModel", "GaugeDrift", "correct_stream", "mask_intervals", "write_sds"]
__version__ = "0.1.0"
