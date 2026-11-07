"""E13: information budget of channel-wave acquisition designs.

The experiment uses linearised Fisher information on exact three-layer Love
dispersion, modal attenuation and group velocity. It is an optimistic bound
because the same physical model is used for data and inversion.
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

import pilot_joint_inversion as pilot


SEED = 20261105
N_MODELS = 140
FREQUENCIES = np.arange(60.0, 920.0, 20.0)
PHASE_NOISE = 0.008
ATTENUATION_NOISE = 0.08
GROUP_NOISE = 0.015
Z90 = 1.6448536269514722


def all_phase_velocities(
    frequency: float,
    h: float,
    beta: float,
) -> list[float]:
    low = beta + 1e-4
    high = min(pilot.TRUE_MODEL["beta_roof"], pilot.TRUE_MODEL["beta_floor"]) - 1e-4
    linear = np.linspace(low, high, 520)
    near_high = high - np.geomspace(1e-6, min(5.0, high - low), 120)
    grid = np.unique(np.concatenate((linear, near_high)))
    omega = 2.0 * math.pi * frequency
    k = omega / grid
    nu = np.sqrt(np.maximum(omega**2 / beta**2 - k**2, 0.0))
    gamma_roof = np.sqrt(
        np.maximum(k**2 - omega**2 / pilot.TRUE_MODEL["beta_roof"] ** 2, 0.0)
    )
    gamma_floor = np.sqrt(
        np.maximum(k**2 - omega**2 / pilot.TRUE_MODEL["beta_floor"] ** 2, 0.0)
    )
    mu_roof = (
        pilot.TRUE_MODEL["rho_roof"] * pilot.TRUE_MODEL["beta_roof"] ** 2
    )
    mu_floor = (
        pilot.TRUE_MODEL["rho_floor"] * pilot.TRUE_MODEL["beta_floor"] ** 2
    )
    mu_coal = pilot.TRUE_MODEL["rho_coal"] * beta**2
    safe_nu = np.maximum(nu, 1e-12)
    values = nu * h - np.arctan(mu_floor * gamma_floor / (mu_coal * safe_nu)) - np.arctan(
        mu_roof * gamma_roof / (mu_coal * safe_nu)
    )
    roots = []
    for mode_order in range(3):
        mode_values = values - mode_order * math.pi
        for index in range(len(grid) - 1):
            left = float(mode_values[index])
            right = float(mode_values[index + 1])
            if not math.isfinite(left) or not math.isfinite(right):
                continue
            if left == 0.0:
                roots.append(float(grid[index]))
                break
            if left * right < 0.0:
                a = float(grid[index])
                b = float(grid[index + 1])
                fa = left
                for _ in range(48):
                    mid = 0.5 * (a + b)
                    fm = (
                        pilot.secular_value(mid, frequency, h, beta)
                        - mode_order * math.pi
                    )
                    if fa * fm <= 0.0:
                        b = mid
                    else:
                        a = mid
                        fa = fm
                roots.append(0.5 * (a + b))
                break
    # The secular equation is solved separately for mode orders 0, 1 and 2.
    # Keep mode 0 as the physical fundamental branch.  Sorting descending would
    # relabel newly appearing higher branches as "mode 0" near their cutoffs.
    roots.sort()
    return roots


def modal_response(
    h: float,
    beta: float,
    q_coal: float,
) -> dict[str, np.ndarray]:
    phase_modes = []
    attenuation_modes = []
    for frequency in FREQUENCIES:
        velocities = all_phase_velocities(float(frequency), h, beta)
        phase_row = np.full(3, np.nan)
        attenuation_row = np.full(3, np.nan)
        for mode_index, velocity in enumerate(velocities[:3]):
            phase_row[mode_index] = velocity
            energy_roof, _, energy_floor = pilot.modal_energy_fractions(
                float(frequency),
                velocity,
                h,
                beta,
            )
            energy_coal = max(1.0 - energy_roof - energy_floor, 0.0)
            attenuation_row[mode_index] = (
                energy_roof / pilot.TRUE_MODEL["q_roof"]
                + energy_coal / q_coal
                + energy_floor / pilot.TRUE_MODEL["q_floor"]
            )
        phase_modes.append(phase_row)
        attenuation_modes.append(attenuation_row)
    phase_modes = np.asarray(phase_modes)
    attenuation_modes = np.asarray(attenuation_modes)
    group_modes = np.full_like(phase_modes, np.nan)
    for mode_index in range(3):
        mode_phase = phase_modes[:, mode_index]
        for index in range(1, len(FREQUENCIES) - 1):
            if not (
                np.isfinite(mode_phase[index - 1])
                and np.isfinite(mode_phase[index])
                and np.isfinite(mode_phase[index + 1])
            ):
                continue
            derivative = (mode_phase[index + 1] - mode_phase[index - 1]) / (
                FREQUENCIES[index + 1] - FREQUENCIES[index - 1]
            )
            denominator = 1.0 - FREQUENCIES[index] * derivative / mode_phase[index]
            if abs(denominator) > 1e-6:
                group_modes[index, mode_index] = (
                    mode_phase[index] / denominator
                )
    return {
        "phase": phase_modes,
        "attenuation": attenuation_modes,
        "group": group_modes,
    }


def response_for_parameters(parameters: np.ndarray):
    return modal_response(
        float(parameters[0]),
        float(parameters[1]),
        float(math.exp(parameters[2])),
    )


def design_components(name: str) -> list[tuple[str, int]]:
    designs = {
        "c0": [("phase", 0)],
        "c0_c1": [("phase", 0), ("phase", 1)],
        "c0_c1_c2": [("phase", 0), ("phase", 1), ("phase", 2)],
        "c0_q0": [("phase", 0), ("attenuation", 0)],
        "c0_c1_q": [
            ("phase", 0),
            ("phase", 1),
            ("attenuation", 0),
            ("attenuation", 1),
        ],
        "c0_c1_group_q": [
            ("phase", 0),
            ("phase", 1),
            ("group", 0),
            ("group", 1),
            ("attenuation", 0),
            ("attenuation", 1),
        ],
    }
    return designs[name]


def design_information(
    center_response,
    plus_responses,
    minus_responses,
    components: list[tuple[str, int]],
) -> tuple[np.ndarray, int]:
    jacobian_parts = []
    for component, mode_index in components:
        center = center_response[component][:, mode_index]
        values = []
        for response in (plus_responses + minus_responses):
            values.append(response[component][:, mode_index])
        mask = np.isfinite(center)
        for value in values:
            mask &= np.isfinite(value)
        if not np.any(mask):
            continue
        if component == "phase":
            sigma = PHASE_NOISE * np.maximum(np.abs(center[mask]), 1.0)
        elif component == "attenuation":
            sigma = ATTENUATION_NOISE * np.maximum(np.abs(center[mask]), 1e-6)
        else:
            sigma = GROUP_NOISE * np.maximum(np.abs(center[mask]), 1.0)
        derivative_columns = []
        steps = np.asarray([0.004, 1.0, 0.008])
        for parameter_index, step in enumerate(steps):
            plus = plus_responses[parameter_index][component][:, mode_index][mask]
            minus = minus_responses[parameter_index][component][:, mode_index][mask]
            derivative_columns.append((plus - minus) / (2.0 * step) / sigma)
        jacobian_parts.append(np.column_stack(derivative_columns))
    if not jacobian_parts:
        return np.zeros((3, 3)), 0
    jacobian = np.vstack(jacobian_parts)
    return jacobian.T @ jacobian, len(jacobian)


def main() -> None:
    rng = np.random.default_rng(SEED)
    designs = (
        "c0",
        "c0_c1",
        "c0_c1_c2",
        "c0_q0",
        "c0_c1_q",
        "c0_c1_group_q",
    )
    rows = []
    prior_std = np.asarray([0.60, 100.0, 0.65])
    prior_precision = np.diag(1.0 / prior_std**2)

    for model_index in range(N_MODELS):
        while True:
            h = float(rng.uniform(1.3, 3.8))
            beta = float(rng.uniform(950.0, 1320.0))
            q = float(math.exp(rng.uniform(math.log(10.0), math.log(80.0))))
            parameters = np.asarray([h, beta, math.log(q)])
            try:
                center = response_for_parameters(parameters)
            except RuntimeError:
                continue
            if np.isfinite(center["phase"][:, 0]).any():
                break
        steps = np.asarray([0.004, 1.0, 0.008])
        plus_responses = []
        minus_responses = []
        for index, step in enumerate(steps):
            plus = parameters.copy()
            minus = parameters.copy()
            plus[index] += step
            minus[index] -= step
            plus_responses.append(response_for_parameters(plus))
            minus_responses.append(response_for_parameters(minus))

        for design_name in designs:
            information, observation_count = design_information(
                center,
                plus_responses,
                minus_responses,
                design_components(design_name),
            )
            posterior_covariance = np.linalg.pinv(information + prior_precision)
            posterior_std = np.sqrt(np.maximum(np.diag(posterior_covariance), 0.0))
            h_width = 2.0 * Z90 * posterior_std[0]
            beta_width = 2.0 * Z90 * posterior_std[1]
            log_q_std = posterior_std[2]
            q_relative_width = 2.0 * Z90 * log_q_std
            rows.append(
                {
                    "model": model_index,
                    "h_true": h,
                    "beta_true": beta,
                    "q_true": q,
                    "design": design_name,
                    "observation_count": observation_count,
                    "h_width_90": h_width,
                    "beta_width_90": beta_width,
                    "q_relative_width_90": q_relative_width,
                    "h_width_accept": h_width <= 0.5,
                }
            )
        if (model_index + 1) % 20 == 0:
            print(f"processed {model_index + 1}/{N_MODELS}")

    (ROOT / "e13_acquisition_information.csv").write_text(
        "model,h_true,beta_true,q_true,design,observation_count,"
        "h_width_90,beta_width_90,q_relative_width_90,h_width_accept\n"
        + "\n".join(
            ",".join(str(row[key]) for key in row)
            for row in rows
        )
        + "\n",
        encoding="utf-8",
    )

    summary = []
    for design in designs:
        subset = [row for row in rows if row["design"] == design]
        summary.append(
            {
                "design": design,
                "median_h_width_90": float(
                    np.median([row["h_width_90"] for row in subset])
                ),
                "p80_h_width_90": float(
                    np.quantile([row["h_width_90"] for row in subset], 0.8)
                ),
                "median_beta_width_90": float(
                    np.median([row["beta_width_90"] for row in subset])
                ),
                "median_q_relative_width_90": float(
                    np.median([row["q_relative_width_90"] for row in subset])
                ),
                "h_width_accept_percent": float(
                    np.mean([row["h_width_accept"] for row in subset]) * 100.0
                ),
                "mean_observation_count": float(
                    np.mean([row["observation_count"] for row in subset])
                ),
            }
        )
    (ROOT / "e13_acquisition_information_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    with (ROOT / "e13_acquisition_information_summary.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
