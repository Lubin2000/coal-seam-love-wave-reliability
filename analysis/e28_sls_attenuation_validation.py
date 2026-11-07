"""E28: broadband viscoelastic Love-wave attenuation validation with SLS.

A five-mechanism standard-linear-solid SH solver is used to generate paired
low-Q and high-Q fields. The coal mode is isolated in f-k space and its
amplitude decay is compared with energy-weighted modal Q theory.
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
import pilot_joint_inversion as pilot


H = 3.8
BETA_COAL = 1120.0
RHO_COAL = 1380.0
VP_COAL = 2000.0
BETA_ROOF = 2400.0
RHO_ROOF = 2450.0
VP_ROOF = 4000.0
BETA_FLOOR = 2650.0
RHO_FLOOR = 2600.0
VP_FLOOR = 4300.0
Q_LOW = 40.0
Q_HIGH = 300.0
Q_ROOF = 200.0
Q_FLOOR = 240.0
SEAM_TOP = 24.0
DX = 0.5
SOURCE_X = 30.0
X_MAX = 230.0
RECEIVER_START = 22.0
RECEIVER_COUNT = 96
RECEIVER_SPACING = 2.0
TOTAL_TIME = 0.30
SOURCE_FREQUENCY = 280.0
N_MECHANISMS = 5
SLS_MIN_FREQUENCY = 80.0
SLS_MAX_FREQUENCY = 550.0


def fit_sls_parameters(q_target: float):
    frequencies = np.geomspace(SLS_MIN_FREQUENCY, SLS_MAX_FREQUENCY, 100)
    centres = np.geomspace(
        SLS_MIN_FREQUENCY,
        SLS_MAX_FREQUENCY,
        N_MECHANISMS,
    )
    tau_sigma = 1.0 / (2.0 * math.pi * centres)
    weights = np.full(N_MECHANISMS, 1.0 / N_MECHANISMS)
    strengths = np.full(N_MECHANISMS, 2.0 / q_target)
    target_inverse = 1.0 / q_target

    def predicted_q_inverse(values):
        modulus = np.zeros_like(frequencies, dtype=complex)
        omega = 2.0 * math.pi * frequencies
        for weight, tau, strength in zip(
            weights,
            tau_sigma,
            values,
            strict=True,
        ):
            modulus += (
                weight
                * (1.0 + 1j * omega * tau * (1.0 + strength))
                / (1.0 + 1j * omega * tau)
            )
        return modulus.imag / modulus.real

    score = float(
        np.mean((predicted_q_inverse(strengths) - target_inverse) ** 2)
    )
    for _ in range(180):
        improved = False
        for mechanism in range(N_MECHANISMS):
            for factor in (0.85, 0.93, 0.97, 1.03, 1.08, 1.18):
                trial = strengths.copy()
                trial[mechanism] *= factor
                trial_score = float(
                    np.mean(
                        (predicted_q_inverse(trial) - target_inverse) ** 2
                    )
                )
                if trial_score < score:
                    strengths = trial
                    score = trial_score
                    improved = True
        if not improved:
            break

    predicted_q = 1.0 / predicted_q_inverse(strengths)
    omega_scale = 2.0 * math.pi * SOURCE_FREQUENCY
    reference_modulus = 0.0j
    for weight, tau, strength in zip(
        weights,
        tau_sigma,
        strengths,
        strict=True,
    ):
        reference_modulus += (
            weight
            * (1.0 + 1j * omega_scale * tau * (1.0 + strength))
            / (1.0 + 1j * omega_scale * tau)
        )
    return {
        "tau_sigma": tau_sigma,
        "tau_epsilon": tau_sigma * (1.0 + strengths),
        "weights": weights,
        "modulus_scale": 1.0 / float(reference_modulus.real),
        "q_min": float(np.min(predicted_q)),
        "q_max": float(np.max(predicted_q)),
        "q_mean": float(np.mean(predicted_q)),
        "relative_error": float(
            np.mean(np.abs(predicted_q - q_target) / q_target)
        ),
    }


SLS_LOW = fit_sls_parameters(Q_LOW)
SLS_HIGH = fit_sls_parameters(Q_HIGH)
SLS_ROOF = fit_sls_parameters(Q_ROOF)
SLS_FLOOR = fit_sls_parameters(Q_FLOOR)
UNIQUE_SLS = {
    Q_LOW: SLS_LOW,
    Q_HIGH: SLS_HIGH,
    Q_ROOF: SLS_ROOF,
    Q_FLOOR: SLS_FLOOR,
}


def harmonic_q_probe(q_target: float, frequency: float = SOURCE_FREQUENCY):
    """Measure the discrete constitutive response under harmonic strain."""
    parameters = UNIQUE_SLS[q_target]
    dt = 0.13 * DX / max(BETA_COAL, BETA_ROOF, BETA_FLOOR)
    duration = 0.08
    steps = int(duration / dt)
    time = np.arange(steps) * dt
    omega = 2.0 * math.pi * frequency
    strain = np.sin(omega * time)
    memory = np.zeros(N_MECHANISMS)
    stress = np.zeros(steps)
    for step in range(steps):
        strain_previous = strain[step - 1] if step > 0 else 0.0
        strain_current = strain[step]
        total_stress = 0.0
        for mechanism, weight in enumerate(parameters["weights"]):
            tau_sigma = float(parameters["tau_sigma"][mechanism])
            tau_epsilon = float(parameters["tau_epsilon"][mechanism])
            auxiliary_modulus = weight * (1.0 - tau_epsilon / tau_sigma)
            factor = dt / (2.0 * tau_sigma)
            memory[mechanism] = (
                memory[mechanism] * (1.0 - factor)
                + factor
                * auxiliary_modulus
                * (strain_previous + strain_current)
            ) / (1.0 + factor)
            total_stress += (
                memory[mechanism]
                + weight * tau_epsilon / tau_sigma * strain_current
            )
        stress[step] = total_stress

    start = steps // 2
    design = np.column_stack(
        (
            np.sin(omega * time[start:]),
            np.cos(omega * time[start:]),
        )
    )
    real_coefficient, imaginary_coefficient = np.linalg.lstsq(
        design,
        stress[start:],
        rcond=None,
    )[0]
    measured_q = float(real_coefficient / imaginary_coefficient)
    return {
        "input_q": q_target,
        "frequency_hz": frequency,
        "discrete_q": measured_q,
        "relative_error_percent": 100.0
        * abs(measured_q - q_target)
        / q_target,
    }


def build_profile(z: np.ndarray, q_coal: float):
    beta = np.full_like(z, BETA_ROOF)
    density = np.full_like(z, RHO_ROOF)
    q = np.full_like(z, Q_ROOF)
    seam_mask = (z >= SEAM_TOP) & (z < SEAM_TOP + H)
    beta[seam_mask] = BETA_COAL
    density[seam_mask] = RHO_COAL
    q[seam_mask] = q_coal
    floor_mask = z >= SEAM_TOP + H
    beta[floor_mask] = BETA_FLOOR
    density[floor_mask] = RHO_FLOOR
    q[floor_mask] = Q_FLOOR
    return beta, density, q


def stress_average(values: np.ndarray, axis: int):
    if axis == 1:
        return 0.5 * (values[:, 1:] + values[:, :-1])
    return 0.5 * (values[1:, :] + values[:-1, :])


def mechanism_coefficients(weight: float, tau_sigma: float, shear_modulus, q_values):
    coefficient = np.zeros_like(shear_modulus)
    tau_eps = np.zeros_like(shear_modulus)
    for q_value, parameters in UNIQUE_SLS.items():
        mask = np.isclose(q_values, q_value)
        index = int(
            np.argmin(np.abs(parameters["tau_sigma"] - tau_sigma))
        )
        tau_eps[mask] = parameters["tau_epsilon"][index]
    relaxed = shear_modulus * weight
    return relaxed, tau_eps


def simulate(q_coal: float):
    z_max = SEAM_TOP + H + 34.0
    x = np.arange(0.0, X_MAX + DX, DX)
    z = np.arange(0.0, z_max + DX, DX)
    nx, nz = len(x), len(z)
    beta_profile, density_profile, q_profile = build_profile(z, q_coal)
    beta_grid = np.tile(beta_profile[:, None], (1, nx))
    density_grid = np.tile(density_profile[:, None], (1, nx))
    q_grid = np.tile(q_profile[:, None], (1, nx))
    mu_grid = density_grid * beta_grid**2

    beta_x = stress_average(beta_grid, axis=1)
    mu_x = stress_average(mu_grid, axis=1)
    q_x = stress_average(q_grid, axis=1)
    beta_z = stress_average(beta_grid, axis=0)
    mu_z = stress_average(mu_grid, axis=0)
    q_z = stress_average(q_grid, axis=0)

    sxz_memory = np.zeros((N_MECHANISMS, nz, nx - 1))
    syz_memory = np.zeros((N_MECHANISMS, nz - 1, nx))
    ex_strain = np.zeros((nz, nx - 1))
    ez_strain = np.zeros((nz - 1, nx))

    velocity = np.zeros((nz, nx))
    max_velocity = max(BETA_COAL, BETA_ROOF, BETA_FLOOR)
    dt = 0.13 * DX / max_velocity
    sample_stride = max(1, int(round(0.0008 / dt)))
    steps = int(TOTAL_TIME / dt)
    times = []
    receiver_data = []

    seam_center = SEAM_TOP + 0.5 * H
    source_ix = int(round(SOURCE_X / DX))
    source_iz = int(round(seam_center / DX))
    offsets = RECEIVER_START + np.arange(RECEIVER_COUNT) * RECEIVER_SPACING
    receiver_x = SOURCE_X + offsets
    receiver_ix = np.clip(np.round(receiver_x / DX).astype(int), 0, nx - 1)
    receiver_iz = source_iz

    sponge_width = 18
    sponge_x = np.minimum(np.arange(nx), np.arange(nx)[::-1]) / sponge_width
    sponge_z = np.minimum(np.arange(nz), np.arange(nz)[::-1]) / sponge_width
    damping = 65.0 * (
        1.0
        - np.clip(np.minimum(sponge_x[None, :], sponge_z[:, None]), 0.0, 1.0)
    ) ** 2

    q_uniques = np.asarray(list(UNIQUE_SLS.keys()))
    tau_eps_x_by_mechanism = []
    tau_eps_z_by_mechanism = []
    auxiliary_x_by_mechanism = []
    auxiliary_z_by_mechanism = []
    instantaneous_x_by_mechanism = []
    instantaneous_z_by_mechanism = []
    nearest_x = np.argmin(
        np.abs(q_x[..., None] - q_uniques),
        axis=-1,
    )
    nearest_z = np.argmin(
        np.abs(q_z[..., None] - q_uniques),
        axis=-1,
    )
    for mechanism in range(N_MECHANISMS):
        tau_eps_x = np.zeros_like(mu_x)
        tau_eps_z = np.zeros_like(mu_z)
        for q_index, q_value in enumerate(q_uniques):
            tau_value = UNIQUE_SLS[float(q_value)]["tau_epsilon"][
                mechanism
            ]
            tau_eps_x[nearest_x == q_index] = tau_value
            tau_eps_z[nearest_z == q_index] = tau_value
        tau_eps_x_by_mechanism.append(tau_eps_x)
        tau_eps_z_by_mechanism.append(tau_eps_z)
        tau_sigma = float(SLS_LOW["tau_sigma"][mechanism])
        weight = float(SLS_LOW["weights"][mechanism])
        q_scale_x = np.ones_like(mu_x)
        q_scale_z = np.ones_like(mu_z)
        for q_index, q_value in enumerate(q_uniques):
            q_scale_x[nearest_x == q_index] = UNIQUE_SLS[
                float(q_value)
            ]["modulus_scale"]
            q_scale_z[nearest_z == q_index] = UNIQUE_SLS[
                float(q_value)
            ]["modulus_scale"]
        instantaneous_x = (
            q_scale_x
            * mu_x
            * weight
            * tau_eps_x
            / tau_sigma
        )
        instantaneous_z = (
            q_scale_z
            * mu_z
            * weight
            * tau_eps_z
            / tau_sigma
        )
        auxiliary_x_by_mechanism.append(
            q_scale_x * mu_x * weight - instantaneous_x
        )
        auxiliary_z_by_mechanism.append(
            q_scale_z * mu_z * weight - instantaneous_z
        )
        instantaneous_x_by_mechanism.append(instantaneous_x)
        instantaneous_z_by_mechanism.append(instantaneous_z)

    source_shift = 1.5 / SOURCE_FREQUENCY
    sxz_total = np.zeros((nz, nx - 1))
    syz_total = np.zeros((nz - 1, nx))
    acceleration = np.zeros_like(velocity)
    for step in range(steps):
        if step > 0 and step % 1000 == 0:
            print(
                f"  step {step}/{steps} "
                f"({100.0 * step / steps:.1f}%)",
                flush=True,
            )
        time = step * dt
        strain_rate_x = (velocity[:, 1:] - velocity[:, :-1]) / DX
        strain_rate_z = (velocity[1:, :] - velocity[:-1, :]) / DX
        ex_strain += dt * strain_rate_x
        ez_strain += dt * strain_rate_z

        sxz_total.fill(0.0)
        syz_total.fill(0.0)
        for mechanism, weight in enumerate(SLS_LOW["weights"]):
            tau_sigma = float(SLS_LOW["tau_sigma"][mechanism])
            auxiliary_x = auxiliary_x_by_mechanism[mechanism]
            auxiliary_z = auxiliary_z_by_mechanism[mechanism]
            midpoint_strain_x = 2.0 * ex_strain - dt * strain_rate_x
            midpoint_strain_z = 2.0 * ez_strain - dt * strain_rate_z
            sxz_memory[mechanism] = (
                sxz_memory[mechanism]
                * (1.0 - dt / (2.0 * tau_sigma))
                + dt
                / (2.0 * tau_sigma)
                * auxiliary_x
                * midpoint_strain_x
            ) / (1.0 + dt / (2.0 * tau_sigma))
            syz_memory[mechanism] = (
                syz_memory[mechanism]
                * (1.0 - dt / (2.0 * tau_sigma))
                + dt
                / (2.0 * tau_sigma)
                * auxiliary_z
                * midpoint_strain_z
            ) / (1.0 + dt / (2.0 * tau_sigma))
            sxz_total += (
                sxz_memory[mechanism]
                + instantaneous_x_by_mechanism[mechanism] * ex_strain
            )
            syz_total += (
                syz_memory[mechanism]
                + instantaneous_z_by_mechanism[mechanism] * ez_strain
            )

        acceleration.fill(0.0)
        acceleration[:, 1:-1] = (
            sxz_total[:, 1:] - sxz_total[:, :-1]
        ) / DX
        acceleration[1:-1, :] += (
            syz_total[1:, :] - syz_total[:-1, :]
        ) / DX
        acceleration /= density_grid

        scaled_time = math.pi * SOURCE_FREQUENCY * (time - source_shift)
        ricker = (1.0 - 2.0 * scaled_time**2) * math.exp(
            -scaled_time**2
        )
        acceleration[source_iz, source_ix] += 8.0e4 * ricker
        velocity += dt * acceleration
        velocity *= np.exp(-damping * dt)
        sxz_memory *= np.exp(-damping[:, :-1] * dt)[None, :, :]
        syz_memory *= np.exp(-damping[:-1, :] * dt)[None, :, :]
        ex_strain *= np.exp(-damping[:, :-1] * dt)
        ez_strain *= np.exp(-damping[:-1, :] * dt)

        if step % sample_stride == 0:
            times.append(time)
            receiver_data.append(
                velocity[receiver_iz, receiver_ix].copy()
            )

    data = np.asarray(receiver_data).T
    return np.asarray(times), offsets, data


def mode_filtered_amplitudes(times, offsets, data):
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
        if frequency < 80.0 or frequency > 500.0:
            continue
        roots = e13.all_phase_velocities(
            float(frequency),
            H,
            BETA_COAL,
        )
        if not roots:
            continue
        target_k = 2.0 * math.pi * frequency / roots[0]
        sigma_k = max(
            0.10 * target_k,
            2.0 * math.pi / (trace_count * spacing),
        )
        frequency_weight = math.exp(
            -0.5
            * (
                (frequency - SOURCE_FREQUENCY)
                / (0.35 * SOURCE_FREQUENCY)
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


def theoretical_q(frequency, q_coal):
    phase = e13.all_phase_velocities(
        frequency,
        H,
        BETA_COAL,
    )[0]
    e0, e1, e2 = pilot.modal_energy_fractions(
        frequency,
        phase,
        H,
        BETA_COAL,
    )
    return 1.0 / (
        e0 / Q_ROOF + e1 / q_coal + e2 / Q_FLOOR
    ), phase


def main():
    print("simulating low-Q SLS field", flush=True)
    times_low, offsets, data_low = simulate(Q_LOW)
    print("simulating high-Q SLS field", flush=True)
    times_high, offsets_high, data_high = simulate(Q_HIGH)
    if not np.allclose(offsets, offsets_high):
        raise RuntimeError("Receiver geometry mismatch")

    amplitude_low = mode_filtered_amplitudes(
        times_low,
        offsets,
        data_low,
    )
    amplitude_high = mode_filtered_amplitudes(
        times_high,
        offsets,
        data_high,
    )
    centers = np.arange(50.0, 161.0, 10.0)
    low_windows = []
    high_windows = []
    for center in centers:
        mask = np.abs(offsets - center) <= 12.0
        low_windows.append(np.sqrt(np.mean(amplitude_low[mask] ** 2)))
        high_windows.append(np.sqrt(np.mean(amplitude_high[mask] ** 2)))
    low_windows = np.asarray(low_windows)
    high_windows = np.asarray(high_windows)
    log_ratio = np.log(
        np.maximum(low_windows, 1e-20)
        / np.maximum(high_windows, 1e-20)
    )
    fit, covariance = np.polyfit(centers, log_ratio, 1, cov=True)
    slope = float(fit[0])
    intercept = float(fit[1])
    fitted = intercept + slope * centers
    residual_sum_squares = float(np.sum((log_ratio - fitted) ** 2))
    total_sum_squares = float(
        np.sum((log_ratio - np.mean(log_ratio)) ** 2)
    )
    slope_standard_error = float(math.sqrt(covariance[0, 0]))

    q_low_theory, phase = theoretical_q(SOURCE_FREQUENCY, Q_LOW)
    q_high_theory, _ = theoretical_q(SOURCE_FREQUENCY, Q_HIGH)
    expected_slope = (
        -math.pi
        * SOURCE_FREQUENCY
        / phase
        * (1.0 / q_low_theory - 1.0 / q_high_theory)
    )
    estimated_q_low = 1.0 / (
        1.0 / q_high_theory
        - slope * phase / (math.pi * SOURCE_FREQUENCY)
    )

    summary = {
        "method": "five-mechanism standard linear solid SH velocity-stress finite difference",
        "source_frequency_hz": SOURCE_FREQUENCY,
        "sls_frequency_range_hz": [
            SLS_MIN_FREQUENCY,
            SLS_MAX_FREQUENCY,
        ],
        "sls_fit": {
            "q40_min_max_mean_error_percent": [
                SLS_LOW["q_min"],
                SLS_LOW["q_max"],
                SLS_LOW["q_mean"],
                100.0 * SLS_LOW["relative_error"],
            ],
            "q300_min_max_mean_error_percent": [
                SLS_HIGH["q_min"],
                SLS_HIGH["q_max"],
                SLS_HIGH["q_mean"],
                100.0 * SLS_HIGH["relative_error"],
            ],
        },
        "coal_q_low_input": Q_LOW,
        "coal_q_high_input": Q_HIGH,
        "effective_q_low_theory": q_low_theory,
        "effective_q_high_theory": q_high_theory,
        "phase_velocity_mps": phase,
        "harmonic_constitutive_check": {
            "q40": harmonic_q_probe(Q_LOW),
            "q300": harmonic_q_probe(Q_HIGH),
        },
        "observed_log_slope": slope,
        "observed_log_slope_standard_error": slope_standard_error,
        "log_slope_fit_r_squared": 1.0
        - residual_sum_squares / max(total_sum_squares, 1e-30),
        "expected_log_slope": expected_slope,
        "estimated_effective_q_low": estimated_q_low,
        "relative_error_percent": 100.0
        * abs(estimated_q_low - q_low_theory)
        / q_low_theory,
        "window_centers_m": centers.tolist(),
        "log_ratio": log_ratio.tolist(),
        "amplitude_low": low_windows.tolist(),
        "amplitude_high": high_windows.tolist(),
    }
    (ROOT / "e28_sls_attenuation_validation_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    with (ROOT / "e28_sls_attenuation_validation.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            (
                "offset_m",
                "amplitude_low",
                "amplitude_high",
                "log_ratio",
            )
        )
        writer.writerows(
            zip(centers, low_windows, high_windows, log_ratio, strict=True)
        )
    np.savez_compressed(
        ROOT / "e28_sls_wavefields.npz",
        times=times_low,
        offsets=offsets,
        data_low=data_low,
        data_high=data_high,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
