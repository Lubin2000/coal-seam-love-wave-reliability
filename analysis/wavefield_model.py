"""Reduced 2D SH finite-difference model for synthetic channel-wave gathers.

The model is intentionally modest in size so that E5 can compare features over
many geological structures. It preserves Love-wave trapping, multimode
interference, geometric spreading, scattering from layer boundaries, and
first-order intrinsic attenuation.
"""

from __future__ import annotations

import math

import numpy as np

from multilayer_model import Halfspace, Layer


def depth_profile(
    layers: list[Layer],
    seam_top: float,
    bottom_padding: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    depths = [0.0]
    velocities = []
    qualities = []
    current = 0.0
    for layer in layers:
        end = current + layer.thickness
        sample_count = max(1, int(round(layer.thickness / 0.25)))
        z = np.linspace(current, end, sample_count, endpoint=False)
        depths.extend(z.tolist())
        velocities.extend([layer.beta] * len(z))
        qualities.extend([layer.q] * len(z))
        current = end
    total = seam_top + current + bottom_padding
    return np.asarray(depths), np.asarray(velocities), np.asarray(qualities), total


def build_vertical_profile(
    layers: list[Layer],
    upper: Halfspace,
    lower: Halfspace,
    seam_top: float,
    z_max: float,
    dz: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z = np.arange(0.0, z_max + dz, dz)
    beta = np.full_like(z, upper.beta, dtype=float)
    quality = np.full_like(z, upper.q, dtype=float)
    current = seam_top
    for layer in layers:
        mask = (z >= current) & (z < current + layer.thickness)
        beta[mask] = layer.beta
        quality[mask] = layer.q
        current += layer.thickness
    mask = z >= current
    beta[mask] = lower.beta
    quality[mask] = lower.q
    beta = np.clip(beta, 700.0, 3000.0)
    quality = np.clip(quality, 8.0, 300.0)
    return z, beta, quality


def simulate_shot(
    layers: list[Layer],
    upper: Halfspace,
    lower: Halfspace,
    seam_top: float = 24.0,
    dx: float = 1.0,
    receiver_spacing: float = 3.0,
    receiver_count: int = 48,
    total_time: float = 0.44,
    source_frequency: float = 300.0,
    noise_level: float = 0.0,
    seed: int = 0,
    x_max: float = 174.0,
    source_x: float = 20.0,
    receiver_start_offset: float = 0.0,
    normalize: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return time, receiver offsets, and normalised receiver gathers."""
    max_velocity = max(upper.beta, lower.beta, *(layer.beta for layer in layers))
    dt = 0.20 * dx / max_velocity
    seam_thickness = sum(layer.thickness for layer in layers)
    bottom_padding = max(28.0, 1.5 * seam_thickness)
    z_max = seam_top + seam_thickness + bottom_padding
    nx = int(math.ceil(x_max / dx)) + 1
    z, beta_profile, quality_profile = build_vertical_profile(
        layers,
        upper,
        lower,
        seam_top,
        z_max,
        dx,
    )
    nz = len(z)
    beta_grid = np.tile(beta_profile[:, None], (1, nx))
    q_grid = np.tile(quality_profile[:, None], (1, nx))
    c2 = beta_grid**2

    x = np.arange(nx) * dx
    seam_center = seam_top + 0.5 * seam_thickness
    source_z = seam_center
    source_ix = int(round(source_x / dx))
    source_iz = int(round(source_z / dx))

    receiver_offsets = (
        source_x
        + receiver_start_offset
        + np.arange(1, receiver_count + 1) * receiver_spacing
    )
    receiver_ix = np.clip(np.round(receiver_offsets / dx).astype(int), 0, nx - 1)
    receiver_iz = int(round(seam_center / dx))

    # Sponge layers suppress boundary wraparound without distorting the central
    # receiver aperture.
    sponge_width = 12
    sponge_x = np.minimum(np.arange(nx), np.arange(nx)[::-1]) / sponge_width
    sponge_z = np.minimum(np.arange(nz), np.arange(nz)[::-1]) / sponge_width
    sponge = 1.0 - np.clip(np.minimum(sponge_x[None, :], sponge_z[:, None]), 0.0, 1.0)
    boundary_damping = 70.0 * sponge**2

    # A simple constant-Q mapping around the dominant source band. This is a
    # first-order approximation used only for relative feature discrimination.
    intrinsic_damping = math.pi * source_frequency / q_grid

    displacement = np.zeros((nz, nx), dtype=np.float64)
    velocity = np.zeros((nz, nx), dtype=np.float64)
    steps = int(total_time / dt)
    sample_stride = max(1, int(round(0.0008 / dt)))
    times = []
    receiver_data = []
    source_shift = 1.5 / source_frequency

    for step in range(steps):
        time = step * dt
        laplacian = np.zeros_like(displacement)
        laplacian[1:-1, 1:-1] = (
            displacement[1:-1, 2:]
            + displacement[1:-1, :-2]
            + displacement[2:, 1:-1]
            + displacement[:-2, 1:-1]
            - 4.0 * displacement[1:-1, 1:-1]
        ) / dx**2
        acceleration = c2 * laplacian
        acceleration -= 2.0 * (intrinsic_damping + boundary_damping) * velocity

        scaled_time = math.pi * source_frequency * (time - source_shift)
        ricker = (1.0 - 2.0 * scaled_time**2) * math.exp(-scaled_time**2)
        acceleration[source_iz, source_ix] += 5.0e4 * ricker

        velocity += dt * acceleration
        displacement += dt * velocity
        displacement *= np.exp(-boundary_damping * dt)
        velocity *= np.exp(-boundary_damping * dt)

        if step % sample_stride == 0:
            times.append(time)
            receiver_data.append(displacement[receiver_iz, receiver_ix].copy())

    data = np.asarray(receiver_data).T
    data -= np.mean(data, axis=1, keepdims=True)
    if normalize:
        scale = max(float(np.max(np.abs(data))), 1e-12)
        data /= scale

    if noise_level > 0.0:
        rng = np.random.default_rng(seed)
        white = rng.normal(size=data.shape)
        coherent = np.sin(
            np.linspace(0.0, 13.0 * math.pi, data.shape[1])[None, :]
            + np.linspace(0.0, 1.4 * math.pi, data.shape[0])[:, None]
        )
        noise = 0.7 * white + 0.3 * coherent
        noise *= noise_level / max(float(np.sqrt(np.mean(noise**2))), 1e-12)
        data += noise

    return np.asarray(times), receiver_offsets - source_x, data


def frequency_velocity_features(
    times: np.ndarray,
    offsets: np.ndarray,
    data: np.ndarray,
    frequency_bins: int = 10,
    velocity_bins: int = 28,
) -> np.ndarray:
    dt = float(np.median(np.diff(times)))
    dx = float(np.median(np.diff(offsets)))
    time_window = np.hanning(data.shape[1])[None, :]
    space_window = np.hanning(data.shape[0])[:, None]
    spectrum = np.fft.fftshift(
        np.fft.fft2(data * time_window * space_window),
        axes=0,
    )
    frequencies = np.fft.fftfreq(data.shape[1], dt)
    wavenumbers = np.fft.fftshift(np.fft.fftfreq(data.shape[0], dx))
    positive_frequency = frequencies > 80.0
    frequencies = frequencies[positive_frequency]
    spectrum = spectrum[:, positive_frequency]

    frequency_edges = np.linspace(100.0, min(720.0, float(np.max(frequencies))), frequency_bins + 1)
    velocity_edges = np.geomspace(850.0, 2400.0, velocity_bins + 1)
    image = np.zeros((frequency_bins, velocity_bins))

    nonzero_wavenumber = np.abs(wavenumbers) > 1e-8
    for frequency_index, frequency in enumerate(frequencies):
        mask_f = (
            (frequency >= frequency_edges[:-1])
            & (frequency < frequency_edges[1:])
        )
        if not np.any(mask_f):
            continue
        bin_index = int(np.argmax(mask_f))
        amplitude = np.abs(spectrum[:, frequency_index]) ** 2
        valid = (
            nonzero_wavenumber
            & (wavenumbers > 0)
            & (amplitude > 0)
        )
        phase_velocities = 2.0 * math.pi * frequency / np.maximum(
            np.abs(wavenumbers[valid]),
            1e-12,
        )
        values = amplitude[valid]
        hist, _ = np.histogram(
            phase_velocities,
            bins=velocity_edges,
            weights=values,
        )
        image[bin_index] += hist

    image /= max(float(np.sum(image)), 1e-12)
    return np.log1p(image * 1e4).ravel()


def wavefield_features(
    times: np.ndarray,
    offsets: np.ndarray,
    data: np.ndarray,
) -> np.ndarray:
    fv = frequency_velocity_features(times, offsets, data)
    singular_values = np.linalg.svd(data, compute_uv=False)
    singular_values /= max(float(np.sum(singular_values)), 1e-12)
    singular = singular_values[:12]
    trace_norm = np.linalg.norm(data, axis=1)
    trace_coherence = np.abs(
        np.corrcoef(data[:, len(data[0]) // 4 : 3 * len(data[0]) // 4])
    )
    trace_coherence = np.nan_to_num(trace_coherence, nan=0.0)
    coherence_features = np.asarray(
        [
            np.mean(np.abs(trace_coherence[np.triu_indices_from(trace_coherence, 1)])),
            np.std(trace_norm),
            np.max(trace_norm),
            np.mean(trace_norm),
        ]
    )
    return np.concatenate((fv, singular, coherence_features))
