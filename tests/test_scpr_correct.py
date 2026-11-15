import numpy as np
from obspy import Stream, Trace, UTCDateTime

from scpr_correct import DriftModel, correct_stream


def test_bundled_models_load():
    for offsets in ("glenn", "masked"):
        m = DriftModel.load("AXCC1", "bundled", offsets)
        assert set(m.gauges) == {"PARO1", "PARO2"}
        d = m.gauges["PARO1"]
        assert d.t_first < d.t_last < d.t_valid


def test_drift_shape():
    d = DriftModel.load("AXCC1", "bundled").gauges["PARO1"]
    t = np.array([d.t_first - 86400, d.t_first, d.t_last, d.t_valid - 1, d.t_valid + 1])
    v = d.delta(t)
    assert v[0] == v[1]                       # held before the first calibration
    assert abs(v[1]) < 1e-6                   # zero at the reference (first calibration)
    assert np.isfinite(v[3]) and np.isnan(v[4])   # not corrected after valid_until
    assert d.sigma([d.t_last])[0] == d.sigma0


def _trace(loc, cha, start, data):
    tr = Trace(np.asarray(data, dtype=float))
    tr.stats.update(dict(network="OO", station="AXCC1", location=loc, channel=cha, sampling_rate=1.0,
                         starttime=UTCDateTime(start)))
    return tr


def test_correction_and_calibration_mask():
    m = DriftModel.load("AXCC1", "bundled")
    start = "2020-05-04T00:00:00"
    p = np.full(86400, 15550e3)
    p[50000:51200] = 15269e3                  # a calibration: pressure drops to the DWT level
    flag = np.zeros(86400)
    flag[49998:51200] = 1
    st = Stream([_trace("40", "LDO", start, p), _trace("41", "LDO", start, p), _trace("40", "LEF", start, flag)])
    out, stats = correct_stream(st, m)
    s40 = [s for s in stats if s["location"] == "40"][0]
    # masked: 60 s before the flag to 1800 s after the calibration
    assert s40["masked_calibration"] == 60 + 1202 + 1800
    t0 = UTCDateTime(start).timestamp
    first = out.select(location="40")[0]
    expected = 15550e3 - m.gauges["PARO1"].delta([t0])[0]
    assert abs(first.data[0] - expected) < 1e-6
    assert first.stats.mseed.dataquality == "Q"


def test_station_block_and_unknown_station():
    m = DriftModel.load("OO.AXCC1", "bundled")
    assert m.id == "OO.AXCC1" and m.gauge_location == {"PARO1": "40", "PARO2": "41"}
    import pytest
    with pytest.raises(LookupError):
        DriftModel.load("CZMID", "bundled")      # no calibration record yet


def test_other_station_traces_ignored():
    m = DriftModel.load("AXCC1", "bundled")
    tr = _trace("40", "LDO", "2020-05-04T00:00:00", np.full(100, 15550e3))
    tr.stats.station = "CZMID"
    out, stats = correct_stream(Stream([tr]), m)
    assert len(out) == 0 and stats == []
