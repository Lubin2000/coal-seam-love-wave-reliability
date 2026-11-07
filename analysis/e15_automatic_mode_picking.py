"""E15: synthetic automated higher-mode ridge picking and information loss."""

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


SEED = 20261105
N_MODELS = 50
VELOCITY_GRID = np.linspace(850.0, 2400.0, 240)
PICK_TOLERANCE = 35.0
Z90 = e13.Z90


def make_dispersion_image(
    phase_modes: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    image = np.zeros((len(e13.FREQUENCIES), len(VELOCITY_GRID)))
    amplitudes = (1.0, 0.58, 0.33)
    widths = (19.0, 23.0, 27.0)
    for frequency_index in range(len(e13.FREQUENCIES)):
        frequency_visibility = rng.uniform(0.65, 1.0)
        if rng.random() < 0.10:
            frequency_visibility *= 0.35
        for mode_index in range(3):
            velocity = phase_modes[frequency_index, mode_index]
            if not np.isfinite(velocity):
                continue
            ridge = np.exp(
                -0.5 * ((VELOCITY_GRID - velocity) / widths[mode_index]) ** 2
            )
            image[frequency_index] += (
                amplitudes[mode_index]
                * frequency_visibility
                * ridge
            )
        # A frequency-coherent body/air-wave-like ridge creates plausible false
        # picks without placing it near the Love modal branches.
        false_ridge = np.exp(
            -0.5
            * (
                (VELOCITY_GRID - (1450.0 + 60.0 * math.sin(frequency_index / 8.0)))
                / 28.0
            )
            ** 2
        )
        image[frequency_index] += 0.20 * false_ridge
        image[frequency_index] += rng.normal(0.0, 0.025, len(VELOCITY_GRID))
    image = np.maximum(image, 0.0)
    image /= np.maximum(np.max(image, axis=1, keepdims=True), 1e-9)
    return image


def pick_ridges(image: np.ndarray) -> list[list[float]]:
    picks = []
    for frequency_index in range(image.shape[0]):
        row = image[frequency_index]
        candidates = []
        for index in range(1, len(row) - 1):
            if row[index] >= row[index - 1] and row[index] >= row[index + 1]:
                if row[index] >= 0.16:
                    candidates.append((VELOCITY_GRID[index], row[index]))
        candidates.sort(key=lambda item: item[0], reverse=True)
        selected = []
        for velocity, amplitude in candidates:
            if all(abs(velocity - previous) > 45.0 for previous in selected):
                selected.append(velocity)
            if len(selected) == 3:
                break
        picks.append(sorted(selected, reverse=True))
    return picks


def match_picks(
    true_modes: np.ndarray,
    picks: list[list[float]],
) -> list[dict]:
    matches = []
    for frequency_index, frequency in enumerate(e13.FREQUENCIES):
        for picked_mode_index, picked_velocity in enumerate(
            picks[frequency_index]
        ):
            if not np.isfinite(picked_velocity):
                continue
            valid_modes = [
                (mode_index, true_modes[frequency_index, mode_index])
                for mode_index in range(3)
                if np.isfinite(true_modes[frequency_index, mode_index])
            ]
            if not valid_modes:
                matches.append(
                    {
                        "frequency_index": frequency_index,
                        "frequency": float(frequency),
                        "picked_mode": picked_mode_index,
                        "picked_velocity": picked_velocity,
                        "true_mode": None,
                        "error": None,
                    }
                )
                continue
            true_mode, true_velocity = min(
                valid_modes,
                key=lambda item: abs(item[1] - picked_velocity),
            )
            error = picked_velocity - true_velocity
            matches.append(
                {
                    "frequency_index": frequency_index,
                    "frequency": float(frequency),
                    "picked_mode": picked_mode_index,
                    "picked_velocity": picked_velocity,
                    "true_mode": int(true_mode),
                    "error": float(error),
                    "correct": abs(error) <= PICK_TOLERANCE,
                }
            )
    return matches


def information_from_picks(
    center,
    plus_responses,
    minus_responses,
    matches,
) -> tuple[np.ndarray, int]:
    columns = []
    for match in matches:
        if not match.get("correct", False):
            continue
        frequency_index = match["frequency_index"]
        mode_index = match["true_mode"]
        center_velocity = center["phase"][frequency_index, mode_index]
        sigma = e13.PHASE_NOISE * center_velocity
        derivative = []
        for parameter_index, step in enumerate((0.004, 1.0, 0.008)):
            plus = plus_responses[parameter_index]["phase"][
                frequency_index, mode_index
            ]
            minus = minus_responses[parameter_index]["phase"][
                frequency_index, mode_index
            ]
            derivative.append((plus - minus) / (2.0 * step) / sigma)
        if not np.all(np.isfinite(derivative)):
            continue
        columns.append(derivative)
    if not columns:
        return np.zeros((3, 3)), 0
    jacobian = np.asarray(columns)
    return jacobian.T @ jacobian, len(columns)


def main() -> None:
    rng = np.random.default_rng(SEED)
    rows = []
    mode_stats = {
        0: {"possible": 0, "detected": 0, "errors": []},
        1: {"possible": 0, "detected": 0, "errors": []},
        2: {"possible": 0, "detected": 0, "errors": []},
    }
    false_picks = 0
    total_picks = 0
    prior_std = np.asarray([0.60, 100.0, 0.65])
    prior_precision = np.diag(1.0 / prior_std**2)

    for model_index in range(N_MODELS):
        while True:
            h = float(rng.uniform(1.3, 3.8))
            beta = float(rng.uniform(950.0, 1320.0))
            q = float(math.exp(rng.uniform(math.log(10.0), math.log(80.0))))
            parameters = np.asarray([h, beta, math.log(q)])
            center = e13.response_for_parameters(parameters)
            if np.isfinite(center["phase"][:, 0]).any():
                break
        plus_responses = []
        minus_responses = []
        for parameter_index, step in enumerate((0.004, 1.0, 0.008)):
            plus = parameters.copy()
            minus = parameters.copy()
            plus[parameter_index] += step
            minus[parameter_index] -= step
            plus_responses.append(e13.response_for_parameters(plus))
            minus_responses.append(e13.response_for_parameters(minus))

        image = make_dispersion_image(center["phase"], rng)
        picks = pick_ridges(image)
        matches = match_picks(center["phase"], picks)
        information, detected_count = information_from_picks(
            center,
            plus_responses,
            minus_responses,
            matches,
        )
        covariance = np.linalg.pinv(information + prior_precision)
        std = np.sqrt(np.maximum(np.diag(covariance), 0.0))

        total_picks += len(matches)
        false_picks += sum(1 for match in matches if not match.get("correct", False))
        for frequency_index in range(len(e13.FREQUENCIES)):
            true_modes = center["phase"][frequency_index]
            frequency_matches = [
                match for match in matches if match["frequency_index"] == frequency_index
            ]
            for mode_index in range(3):
                if not np.isfinite(true_modes[mode_index]):
                    continue
                mode_stats[mode_index]["possible"] += 1
                candidate = min(
                    frequency_matches,
                    key=lambda item: abs(
                        item["picked_velocity"] - true_modes[mode_index]
                    ),
                    default=None,
                )
                if candidate is not None and abs(
                    candidate["picked_velocity"] - true_modes[mode_index]
                ) <= PICK_TOLERANCE:
                    mode_stats[mode_index]["detected"] += 1
                    mode_stats[mode_index]["errors"].append(
                        float(candidate["picked_velocity"] - true_modes[mode_index])
                    )

        rows.append(
            {
                "model": model_index,
                "h_true": h,
                "beta_true": beta,
                "q_true": q,
                "detected_picks": detected_count,
                "total_picks": len(matches),
                "false_picks": sum(
                    1 for match in matches if not match.get("correct", False)
                ),
                "h_width_90": 2.0 * Z90 * std[0],
                "beta_width_90": 2.0 * Z90 * std[1],
                "q_relative_width_90": 2.0 * Z90 * std[2],
            }
        )
        if (model_index + 1) % 10 == 0:
            print(f"processed {model_index + 1}/{N_MODELS}")

    mode_summary = {}
    for mode_index, stats in mode_stats.items():
        mode_summary[f"mode{mode_index}"] = {
            "possible_count": stats["possible"],
            "detected_count": stats["detected"],
            "detection_percent": stats["detected"]
            / max(stats["possible"], 1)
            * 100.0,
            "error_mae_mps": float(np.mean(np.abs(stats["errors"])))
            if stats["errors"]
            else None,
            "error_rmse_mps": float(np.sqrt(np.mean(np.square(stats["errors"]))))
            if stats["errors"]
            else None,
        }

    summary = {
        "settings": {
            "n_models": N_MODELS,
            "pick_tolerance_mps": PICK_TOLERANCE,
            "image_frequency_points": len(e13.FREQUENCIES),
            "image_velocity_bins": len(VELOCITY_GRID),
        },
        "picking": {
            "mode_statistics": mode_summary,
            "false_pick_percent": false_picks / max(total_picks, 1) * 100.0,
        },
        "information_after_picking": {
            "median_h_width_90": float(
                np.median([row["h_width_90"] for row in rows])
            ),
            "p90_h_width_90": float(
                np.quantile([row["h_width_90"] for row in rows], 0.90)
            ),
            "h_accept_percent": float(
                np.mean([row["h_width_90"] <= 0.5 for row in rows]) * 100.0
            ),
            "median_beta_width_90": float(
                np.median([row["beta_width_90"] for row in rows])
            ),
            "median_q_relative_width_90": float(
                np.median([row["q_relative_width_90"] for row in rows])
            ),
        },
    }
    with (ROOT / "e15_automatic_picking_runs.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (ROOT / "e15_automatic_picking_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
