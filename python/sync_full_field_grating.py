import os
from typing import Tuple, Optional

import numpy as np

try:  # optional scipy dependency for reading MATLAB files
    from scipy.io import loadmat
except Exception:  # pragma: no cover - scipy may be unavailable
    loadmat = None


def linear_interp(x: np.ndarray, xp: np.ndarray, fp: np.ndarray) -> np.ndarray:
    """Linearly interpolate with extrapolation."""
    x = np.asarray(x)
    xp = np.asarray(xp)
    fp = np.asarray(fp)

    y = np.interp(x, xp, fp)

    left = x < xp[0]
    right = x > xp[-1]

    if np.any(left):
        slope = (fp[1] - fp[0]) / (xp[1] - xp[0])
        y[left] = fp[0] + slope * (x[left] - xp[0])

    if np.any(right):
        slope = (fp[-1] - fp[-2]) / (xp[-1] - xp[-2])
        y[right] = fp[-1] + slope * (x[right] - xp[-1])

    return y


def load_binary(file_path: str, n_channels: int, channel: int,
                 dtype: np.dtype = np.int16) -> np.ndarray:
    """Load a single channel from an interleaved binary file.

    Parameters
    ----------
    file_path : str
        Path to ``continuous.dat``.
    n_channels : int
        Total number of channels in the file.
    channel : int
        Channel index (0-based) to extract.
    dtype : np.dtype, optional
        Data type of the file, by default ``np.int16``.

    Returns
    -------
    np.ndarray
        1-D array containing the selected channel.
    """
    data = np.memmap(file_path, dtype=dtype)
    data = data.reshape(-1, n_channels)
    return np.array(data[:, channel], copy=True)


def threshold(time: np.ndarray, signal: np.ndarray, comparator: str,
              thresh: float, min_interval: float,
              max_interval: Optional[float] = None) -> np.ndarray:
    """Detect threshold crossings and return pulse periods.

    Parameters
    ----------
    time : np.ndarray
        Time vector corresponding to ``signal``.
    signal : np.ndarray
        Input signal.
    comparator : str
        ``'>'`` for pulses above ``thresh`` or ``'<'`` for pulses below
        ``thresh``.
    thresh : float
        Threshold value.
    min_interval : float
        Minimum pulse width (seconds).
    max_interval : float, optional
        Maximum pulse width (seconds).  If ``None`` no upper bound is used.

    Returns
    -------
    np.ndarray
        ``(n, 2)`` array with rising and falling edge times for each pulse.
    """
    if comparator == '>':
        mask = signal > thresh
    else:
        mask = signal < thresh

    edges = np.diff(mask.astype(np.int8))
    rising = np.where(edges == 1)[0] + 1
    falling = np.where(edges == -1)[0] + 1

    periods = []
    i = j = 0
    while i < len(rising) and j < len(falling):
        if rising[i] < falling[j]:
            start = time[rising[i]]
            end = time[falling[j]]
            width = end - start
            if width >= min_interval and (max_interval is None or width <= max_interval):
                periods.append((start, end))
            i += 1
            j += 1
        else:
            j += 1
    return np.array(periods)


def sync_full_field_grating(cf: str, stim_path: Optional[str],
                            adc_file: str, probe_file: str,
                            sr: int = 30000, sweep_time: float = 20.0,
                            nchan_adc: int = 12, nchan_probe: int = 385,
                            d385_start_time: float = 0.0) -> Tuple[Optional[dict], np.ndarray]:
    """Synchronize photodiode pulses with probe time.

    This function mirrors the behaviour of the MATLAB function
    ``Sync_FullFieldGrating`` and aligns photodiode events recorded on the
    ADC stream with the function generator pulses present in both the ADC
    and probe data streams.

    Parameters
    ----------
    cf : str
        Unused parameter kept for API compatibility.
    stim_path : str, optional
        Path to a stimulus ``.mat`` file containing a ``trials`` variable.
    adc_file : str
        Path to the ``continuous.dat`` file containing the ADC channels.
    probe_file : str
        Path to the ``continuous.dat`` file for the probe channels.
    sr : int, optional
        Sampling rate (Hz), by default ``30000``.
    sweep_time : float, optional
        Duration of each sweep in seconds, by default ``20``.
    nchan_adc : int, optional
        Number of ADC channels in ``adc_file``.
    nchan_probe : int, optional
        Number of channels in ``probe_file``.
    d385_start_time : float, optional
        Start time of the recording relative to a concatenated file.

    Returns
    -------
    trials : dict or ``None``
        Contents of ``trials`` variable from ``stim_path`` if present.
    Vq_periods : np.ndarray
        ``(n, 2)`` array with start and end times (in probe time) of each
        photodiode pulse.
    """
    trials = None
    if stim_path and os.path.exists(stim_path):
        if loadmat is None:
            raise ImportError('scipy is required to load MATLAB files')
        stim_data = loadmat(stim_path)
        trials = stim_data.get('trials', stim_data)

    d0 = load_binary(adc_file, nchan_adc, 0)
    d1 = load_binary(adc_file, nchan_adc, 1)
    d385 = load_binary(probe_file, nchan_probe, 384)

    d1_dur = len(d1) / sr

    td0 = np.arange(len(d0)) / sr
    td1 = np.arange(len(d1)) / sr
    td385 = np.arange(len(d385)) / sr

    periods0 = threshold(td0, d0.astype(float), '>', 10000, 0.1)
    ns = int(d1_dur / sweep_time)
    blocks1 = np.arange(0, ns * sweep_time, sweep_time)
    blocks2 = blocks1 + sweep_time
    p0_diff = periods0[:, 1] - periods0[:, 0]
    nperiods = [np.sum((periods0[:, 0] > b1) & (periods0[:, 0] <= b2))
                for b1, b2 in zip(blocks1, blocks2)]
    bf = np.cumsum(nperiods)
    bs = np.concatenate(([0], bf[:-1]))
    start_pts = []
    for b_start, b_end in zip(bs, bf):
        window = p0_diff[b_start:b_end]
        idx_local_min = np.argmin(window) + b_start
        idx = idx_local_min + 1
        w = p0_diff[idx]
        while idx + 1 < len(p0_diff) and p0_diff[idx + 1] > w:
            idx += 1
            w = p0_diff[idx]
        start_pts.append(idx)
    start_pts = np.array(sorted(set(start_pts)))
    st_adc0 = periods0[start_pts[:-1], 0]
    et_adc0 = periods0[start_pts[1:], 0]

    periods1 = threshold(td1, d1.astype(float), '<', 20000, 0.001, 0.001)
    widths = periods1[:, 1] - periods1[:, 0]
    periods1 = periods1[widths >= 0.0002]
    periods1_vector = periods1.reshape(-1)

    periods385 = threshold(td385, d385.astype(float), '>', 0.8, 0.1)
    p385_diff = periods385[:, 1] - periods385[:, 0]
    nperiods_probe = [np.sum((periods385[:, 0] > b1) & (periods385[:, 0] <= b2))
                      for b1, b2 in zip(blocks1, blocks2)]
    bf = np.cumsum(nperiods_probe)
    bs = np.concatenate(([0], bf[:-1]))
    start_pts_probe = []
    for b_start, b_end in zip(bs, bf):
        window = p385_diff[b_start:min(b_end + 1, len(p385_diff))]
        idx_local_min = np.argmin(window) + b_start
        idx = idx_local_min + 1
        w = p385_diff[idx]
        while idx + 1 < len(p385_diff) and p385_diff[idx + 1] > w:
            idx += 1
            w = p385_diff[idx]
        start_pts_probe.append(idx)
    start_pts_probe = np.array(sorted(set(start_pts_probe)))
    st_probe = periods385[start_pts_probe[:-1], 0]
    et_probe = periods385[start_pts_probe[1:], 0]

    npulses = min(len(st_adc0), len(st_probe))
    st_adc0 = st_adc0[:npulses]
    st_probe = st_probe[:npulses]
    et_adc0 = et_adc0[:npulses]
    et_probe = et_probe[:npulses]

    Vq = linear_interp(periods1_vector, st_adc0, st_probe)
    Vq = Vq + d385_start_time
    Vq_periods = np.column_stack([Vq[0:-1:2], Vq[1::2]])

    return trials, Vq_periods
