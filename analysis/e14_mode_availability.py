"""E14: information retention under higher-mode dropout."""

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


SEED = 20261105
N_MODELS = 100
MASK_REPLICATES = 4
Z90 = e13.Z90


SCENARIOS = {
    "fundamental_only": {
        "components": (("phase", 0), ("attenuation", 0), ("group", 0)),
        "probabilities": (0.90, 0.0, 0.0),
    },
    "ideal_all_modes": {
        "components": (
            ("phase", 0),
            ("phase", 1),
            ("phase", 2),
            ("attenuation", 0),
            ("attenuation", 1),
            ("attenuation", 2),
            ("group", 0),
            ("group", 1),
            ("group", 2),
        ),
        "probabilities": (1.0, 1.0, 1.0),
    },
    "good_field_multimode": {
        "components": (
            ("phase", 0),
            ("phase", 1),
            ("phase", 2),
            ("attenuation", 0),
            ("attenuation", 1),
            ("attenuation", 2),
            ("group", 0),
            ("group", 1),
            ("group", 2),
        ),
        "probabilities": (0.95, 0.75, 0.40),
    },
    "poor_field_multimode": {
        "components": (
            ("phase", 0),
            ("phase", 1),
            ("phase", 2),
            ("attenuation", 0),
            ("attenuation", 1),
            ("attenuation", 2),
        ),
        "probabilities": (0.82, 0.45, 0.15),
    },
    "poor_phase_only": {
        "components": (("phase", 0), ("phase", 1), ("phase", 2)),
        "probabilities": (0.82, 0.45, 0.15),
    },
}


def masked_information(
    center_response,
    plus_responses,
    minus_responses,
    components,
    probabilities,
    rng: np.random.Generator,
) -> tuple[np.ndarray, int, dict[int, int]]:
    jacobian_parts = []
    mode_counts = {0: 0, 1: 0, 2: 0}
    steps = np.asarray([0.004, 1.0, 0.008])

    for component, mode_index in components:
        center = center_response[component][:, mode_index]
        mask = np.isfinite(center)
        random_mask = rng.random(len(center)) < probabilities[mode_index]
        mask &= random_mask
        for response in plus_responses + minus_responses:
            mask &= np.isfinite(response[component][:, mode_index])
        if not np.any(mask):
            continue

        if component == "phase":
            sigma = e13.PHASE_NOISE * np.maximum(np.abs(center[mask]), 1.0)
        elif component == "attenuation":
            sigma = e13.ATTENUATION_NOISE * np.maximum(
                np.abs(center[mask]),
                1e-6,
            )
        else:
            sigma = e13.GROUP_NOISE * np.maximum(
                np.abs(center[mask]),
                1.0,
            )

        columns = []
        for parameter_index, step in enumerate(steps):
            plus = plus_responses[parameter_index][component][:, mode_index][mask]
            minus = minus_responses[parameter_index][component][:, mode_index][mask]
            columns.append((plus - minus) / (2.0 * step) / sigma)
        jacobian_parts.append(np.column_stack(columns))
        mode_counts[mode_index] += int(np.sum(mask))

    if not jacobian_parts:
        return np.zeros((3, 3)), 0, mode_counts
    jacobian = np.vstack(jacobian_parts)
    return jacobian.T @ jacobian, len(jacobian), mode_counts


def main() -> None:
    rng = np.random.default_rng(SEED)
    prior_std = np.asarray([0.60, 100.0, 0.65])
    prior_precision = np.diag(1.0 / prior_std**2)
    rows = []
    mode_fraction_rows = []

    for model_index in range(N_MODELS):
        while True:
            h = float(rng.uniform(1.3, 3.8))
            beta = float(rng.uniform(950.0, 1320.0))
            q = float(math.exp(rng.uniform(math.log(10.0), math.log(80.0))))
            parameters = np.asarray([h, beta, math.log(q)])
            center = e13.response_for_parameters(parameters)
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
            plus_responses.append(e13.response_for_parameters(plus))
            minus_responses.append(e13.response_for_parameters(minus))

        available_mode_counts = [
            int(np.sum(np.isfinite(center["phase"][:, mode])))
            for mode in range(3)
        ]
        mode_fraction_rows.append(
            {
                "model": model_index,
                "mode0_available_percent": available_mode_counts[0]
                / len(e13.FREQUENCIES)
                * 100.0,
                "mode1_available_percent": available_mode_counts[1]
                / len(e13.FREQUENCIES)
                * 100.0,
                "mode2_available_percent": available_mode_counts[2]
                / len(e13.FREQUENCIES)
                * 100.0,
            }
        )

        for scenario_name, scenario in SCENARIOS.items():
            for mask_replicate in range(MASK_REPLICATES):
                information, observation_count, counts = masked_information(
                    center,
                    plus_responses,
                    minus_responses,
                    scenario["components"],
                    scenario["probabilities"],
                    rng,
                )
                covariance = np.linalg.pinv(information + prior_precision)
                std = np.sqrt(np.maximum(np.diag(covariance), 0.0))
                rows.append(
                    {
                        "model": model_index,
                        "scenario": scenario_name,
                        "mask_replicate": mask_replicate,
                        "h_true": h,
                        "beta_true": beta,
                        "q_true": q,
                        "observation_count": observation_count,
                        "mode0_observations": counts[0],
                        "mode1_observations": counts[1],
                        "mode2_observations": counts[2],
                        "h_width_90": 2.0 * Z90 * std[0],
                        "beta_width_90": 2.0 * Z90 * std[1],
                        "q_relative_width_90": 2.0 * Z90 * std[2],
                    }
                )
        if (model_index + 1) % 20 == 0:
            print(f"processed {model_index + 1}/{N_MODELS}")

    summary = []
    for scenario in SCENARIOS:
        subset = [row for row in rows if row["scenario"] == scenario]
        summary.append(
            {
                "scenario": scenario,
                "median_h_width_90": float(
                    np.median([row["h_width_90"] for row in subset])
                ),
                "p90_h_width_90": float(
                    np.quantile([row["h_width_90"] for row in subset], 0.90)
                ),
                "h_accept_percent": float(
                    np.mean([row["h_width_90"] <= 0.5 for row in subset])
                    * 100.0
                ),
                "median_beta_width_90": float(
                    np.median([row["beta_width_90"] for row in subset])
                ),
                "median_q_relative_width_90": float(
                    np.median([row["q_relative_width_90"] for row in subset])
                ),
                "mean_mode0_observations": float(
                    np.mean([row["mode0_observations"] for row in subset])
                ),
                "mean_mode1_observations": float(
                    np.mean([row["mode1_observations"] for row in subset])
                ),
                "mean_mode2_observations": float(
                    np.mean([row["mode2_observations"] for row in subset])
                ),
            }
        )

    with (ROOT / "e14_mode_availability_runs.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    with (ROOT / "e14_mode_availability_summary.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)
    (ROOT / "e14_mode_availability_summary.json").write_text(
        json.dumps(
            {
                "settings": {
                    "n_models": N_MODELS,
                    "mask_replicates": MASK_REPLICATES,
                    "frequency_hz": [
                        float(e13.FREQUENCIES[0]),
                        float(e13.FREQUENCIES[-1]),
                    ],
                },
                "mean_mode_availability_percent": {
                    "mode0": float(
                        np.mean([row["mode0_available_percent"] for row in mode_fraction_rows])
                    ),
                    "mode1": float(
                        np.mean([row["mode1_available_percent"] for row in mode_fraction_rows])
                    ),
                    "mode2": float(
                        np.mean([row["mode2_available_percent"] for row in mode_fraction_rows])
                    ),
                },
                "summary": summary,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "mode_availability": {
                    "mode0": float(
                        np.mean([row["mode0_available_percent"] for row in mode_fraction_rows])
                    ),
                    "mode1": float(
                        np.mean([row["mode1_available_percent"] for row in mode_fraction_rows])
                    ),
                    "mode2": float(
                        np.mean([row["mode2_available_percent"] for row in mode_fraction_rows])
                    ),
                },
                "summary": summary,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
