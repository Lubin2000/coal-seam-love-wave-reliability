"""E23: exploratory forward check of attenuation using paired SH models.

Two otherwise identical wavefields are generated with low and high coal Q.
Their frequency-domain fundamental-mode amplitudes are compared in sliding
offset windows. The logarithmic amplitude ratio removes source amplitude and
geometric spreading and is compared with modal-energy weighted Q theory.

The present time-domain damping is a single-frequency constant-Q
approximation. It does not yet include a broadband relaxation mechanism and
does not calibrate source coupling. Results are therefore diagnostic rather
than a validated Q inversion.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import e13_acquisition_information as e13
import multilayer_model as ml
import pilot_joint_inversion as pilot
import wavefield_model as wm


H = 3.8
BETA_COAL = 1120.0
RHO_COAL = 1380.0
Q_LOW = 40.0
Q_HIGH = 300.0
BETA_ROOF = 2400.0
RHO_ROOF = 2450.0
Q_ROOF = 100.0
BETA_FLOOR = 2650.0
RHO_FLOOR = 2600.0
Q_FLOOR = 120.0
SOURCE_FREQUENCIES = (220.0, 320.0)
DX = 0.5
RECEIVER_COUNT = 100
RECEIVER_SPACING = 2.0
RECEIVER_START_OFFSET = 20.0
X_MAX = 300.0
SOURCE_X = 30.0
TOTAL_TIME = 0.32


def simulate(q_coal: float, source_frequency: float):
    layers = [ml.Layer(H, BETA_COAL, RHO_COAL, q_coal)]
    upper = ml.Halfspace(BETA_ROOF, RHO_ROOF, Q_ROOF)
    lower = ml.Halfspace(BETA_FLOOR, RHO_FLOOR, Q_FLOOR)
    return wm.simulate_shot(
        layers,
        upper,
        lower,
        seam_top=24.0,
        dx=DX,
        receiver_spacing=RECEIVER_SPACING,
        receiver_count=RECEIVER_COUNT,
        total_time=TOTAL_TIME,
        source_frequency=source_frequency,
        noise_level=0.0,
        x_max=X_MAX,
        source_x=SOURCE_X,
        receiver_start_offset=RECEIVER_START_OFFSET,
        normalize=False,
    )


def windowed_amplitudes(
    times: np.ndarray,
    offsets: np.ndarray,
    data: np.ndarray,
    frequency: float,
    phase_velocity: float,
    window_centers: np.ndarray,
    half_width: float = 18.0,
) -> np.ndarray:
    sample_dt = float(np.median(np.diff(times)))
    trace_count, sample_count = data.shape
    nfft = 1 << int(math.ceil(math.log2(max(4096, sample_count * 8))))
    spectrum = np.fft.rfft(
        data.T * np.hanning(sample_count)[:, None],
        n=nfft,
        axis=0,
    )
    frequencies = np.fft.rfftfreq(nfft, sample_dt)
    spectrum_index = int(np.argmin(np.abs(frequencies - frequency)))
    values = []
    for center in window_centers:
        mask = np.abs(offsets - center) <= half_width
        local_offsets = offsets[mask]
        local_spectrum = spectrum[spectrum_index, mask]
        phasor = np.exp(
            1j * 2.0 * math.pi * frequency * local_offsets / phase_velocity
        )
        values.append(float(abs(np.sum(local_spectrum * phasor))))
    return np.asarray(values)


def mode_filtered_amplitudes(
    times: np.ndarray,
    offsets: np.ndarray,
    data: np.ndarray,
    source_frequency: float,
) -> np.ndarray:
    """Isolate the fundamental mode in f-k space and return its trace RMS."""
    sample_dt = float(np.median(np.diff(times)))
    spacing = float(np.median(np.diff(offsets)))
    trace_count, sample_count = data.shape
    nfft_x = 1 << int(math.ceil(math.log2(max(128, trace_count * 2))))
    nfft_t = 1 << int(math.ceil(math.log2(max(4096, sample_count * 8))))
    spectrum = np.fft.fft2(data, s=(nfft_x, nfft_t))
    frequencies = np.fft.fftfreq(nfft_t, sample_dt)
    wavenumbers = 2.0 * math.pi * np.fft.fftfreq(nfft_x, spacing)
    mask = np.zeros_like(spectrum)

    for frequency_index, frequency in enumerate(np.abs(frequencies)):
        if frequency < 80.0 or frequency > min(500.0, 1.55 * source_frequency):
            continue
        roots = e13.all_phase_velocities(
            float(frequency),
            H,
            BETA_COAL,
        )
        if not roots:
            continue
        phase_velocity = roots[0]
        target_k = 2.0 * math.pi * frequency / phase_velocity
        sigma_k = max(0.08 * target_k, 2.0 * math.pi / (trace_count * spacing))
        frequency_weight = math.exp(
            -0.5
            * (
                (frequency - source_frequency)
                / (0.32 * source_frequency)
            )
            ** 2
        )
        mask[:, frequency_index] = frequency_weight * np.exp(
            -0.5 * ((np.abs(wavenumbers) - target_k) / sigma_k) ** 2
        )

    filtered = np.real(
        np.fft.ifft2(spectrum * mask)
    )[:trace_count, :sample_count]
    return np.sqrt(np.mean(filtered**2, axis=1))


def analytical_effective_q(frequency: float, phase_velocity: float):
    roof_fraction, coal_fraction, floor_fraction = (
        pilot.modal_energy_fractions(
            frequency,
            phase_velocity,
            H,
            BETA_COAL,
        )
    )
    q_low_effective_inverse = (
        roof_fraction / Q_ROOF
        + coal_fraction / Q_LOW
        + floor_fraction / Q_FLOOR
    )
    q_high_effective_inverse = (
        roof_fraction / Q_ROOF
        + coal_fraction / Q_HIGH
        + floor_fraction / Q_FLOOR
    )
    coal_energy = float(coal_fraction)
    return (
        1.0 / q_low_effective_inverse,
        1.0 / q_high_effective_inverse,
        coal_energy,
    )


def main() -> None:
    rows = []
    window_centers = np.arange(50.0, 181.0, 15.0)
    for source_frequency in SOURCE_FREQUENCIES:
        print(f"E23 source {source_frequency:.0f} Hz", flush=True)
        times_high, offsets, data_high = simulate(Q_HIGH, source_frequency)
        times_low, offsets_low, data_low = simulate(Q_LOW, source_frequency)
        if not np.allclose(offsets, offsets_low):
            raise RuntimeError("Receiver geometry changed between Q runs")

        phase = e13.all_phase_velocities(
            source_frequency,
            H,
            BETA_COAL,
        )[0]
        amplitude_high_full = mode_filtered_amplitudes(
            times_high,
            offsets,
            data_high,
            source_frequency,
        )
        amplitude_low_full = mode_filtered_amplitudes(
            times_low,
            offsets,
            data_low,
            source_frequency,
        )
        amplitude_high = np.asarray(
            [
                np.sqrt(
                    np.mean(
                        amplitude_high_full[
                            np.abs(offsets - center) <= 10.0
                        ]
                        ** 2
                    )
                )
                for center in window_centers
            ]
        )
        amplitude_low = np.asarray(
            [
                np.sqrt(
                    np.mean(
                        amplitude_low_full[
                            np.abs(offsets - center) <= 10.0
                        ]
                        ** 2
                    )
                )
                for center in window_centers
            ]
        )
        log_ratio = np.log(
            np.maximum(amplitude_low, 1.0e-20)
            / np.maximum(amplitude_high, 1.0e-20)
        )
        # Windows beyond 125 m are increasingly affected by the absorbing
        # boundary and by interference between body and guided arrivals.
        stable = window_centers <= 140.0
        effective_distance = float(np.mean(window_centers[stable]))
        effective_log_ratio = float(np.mean(log_ratio[stable]))
        slope = effective_log_ratio / effective_distance

        q_low_effective, q_high_effective, coal_energy = analytical_effective_q(
            source_frequency,
            phase,
        )
        expected_slope = (
            -math.pi
            * source_frequency
            / phase
            * (1.0 / q_low_effective - 1.0 / q_high_effective)
        )
        estimated_low_inverse = (
            1.0 / q_high_effective
            - slope * phase / (math.pi * source_frequency)
        )
        estimated_q_low = 1.0 / max(estimated_low_inverse, 1.0e-12)

        rows.append(
            {
                "source_frequency_hz": source_frequency,
                "h_m": H,
                "phase_velocity_mps": phase,
                "q_low_coal_input": Q_LOW,
                "q_high_coal_input": Q_HIGH,
                "coal_energy_fraction": coal_energy,
                "q_low_effective_theory": q_low_effective,
                "q_high_effective_theory": q_high_effective,
                "expected_log_slope": expected_slope,
                "observed_log_slope": slope,
                "naive_estimated_q_low_effective": estimated_q_low,
                "relative_q_error_percent": 100.0
                * abs(estimated_q_low - q_low_effective)
                / q_low_effective,
                "amplitude_high": json.dumps(amplitude_high.tolist()),
                "amplitude_low": json.dumps(amplitude_low.tolist()),
                "log_ratio": json.dumps(log_ratio.tolist()),
                "window_centers_m": json.dumps(window_centers.tolist()),
            }
        )
        print(json.dumps(rows[-1], ensure_ascii=False, indent=2), flush=True)

    with (ROOT / "e23_viscoelastic_attenuation_validation.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "method": {
            "wavefield": "2D SH finite difference with first-order constant-Q damping",
            "mode_estimation": "fundamental-mode f-k ridge filtering followed by trace-RMS amplitude fitting",
            "low_q_coal": Q_LOW,
            "high_q_coal": Q_HIGH,
            "source_frequencies_hz": list(SOURCE_FREQUENCIES),
            "window_centers_m": window_centers.tolist(),
            "window_half_width_m": 18.0,
            "interpretation": "Exploratory only: the 220-Hz effective-Q estimate is approximate, while 320-Hz mode interference remains significant.",
            "grid_spacing_m": DX,
            "receiver_count": RECEIVER_COUNT,
            "receiver_spacing_m": RECEIVER_SPACING,
        },
        "results": rows,
    }
    (ROOT / "e23_viscoelastic_attenuation_validation_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print("E23 validation complete.", flush=True)


if __name__ == "__main__":
    main()
