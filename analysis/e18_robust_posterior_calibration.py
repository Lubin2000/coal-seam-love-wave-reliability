"""E18: calibration of robust inversion intervals under mode-pick outliers."""

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
import e15_automatic_mode_picking as e15
import e16_robust_mode_tracking as e16
import e17_robust_likelihood as e17


SEED = 20261105
N_CALIBRATION = 35
N_TEST = 35
OUTLIER_FRACTIONS = (0.0, 0.10, 0.20, 0.30)
Z90 = e13.Z90


def robust_fit(jacobian, residual):
    delta = np.zeros(3)
    for _ in range(30):
        standardized = residual + jacobian @ delta
        weights = 5.0 / (4.0 + standardized**2)
        weighted_jacobian = jacobian * np.sqrt(weights[:, None])
        information = weighted_jacobian.T @ weighted_jacobian + np.diag(
            1.0 / e17.PRIOR_STD**2
        )
        rhs = -weighted_jacobian.T @ (np.sqrt(weights) * residual)
        new_delta = np.linalg.solve(information, rhs)
        if np.linalg.norm(new_delta - delta) < 1e-9:
            delta = new_delta
            break
        delta = new_delta
    standardized = residual + jacobian @ delta
    weights = 5.0 / (4.0 + standardized**2)
    weighted_jacobian = jacobian * np.sqrt(weights[:, None])
    covariance = np.linalg.pinv(
        weighted_jacobian.T @ weighted_jacobian
        + np.diag(1.0 / e17.PRIOR_STD**2)
    )
    std = np.sqrt(np.maximum(np.diag(covariance), 1e-12))
    outlier_score = 1.0 - float(np.mean(weights))
    return delta, std, outlier_score, standardized


def build_observation_set(
    rng,
    model_count,
    seed_offset,
):
    rng = np.random.default_rng(SEED + seed_offset)
    rows = []
    for model_index in range(model_count):
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

        image = e15.make_dispersion_image(center["phase"], rng)
        candidates = e16.detect_candidate_peaks(image)
        picks = e16.viterbi_track(candidates)
        observations = e17.build_observation_rows(center, picks, rng)
        if len(observations) < 6:
            continue
        true_params = np.asarray([h, beta, q])
        for fraction in OUTLIER_FRACTIONS:
            perturbed = e17.inject_outliers(rng, observations, fraction)
            jacobian, residual, _ = e17.build_design_matrices(
                center,
                plus_responses,
                minus_responses,
                perturbed,
                include_q=True,
            )
            if len(residual) < 3:
                continue
            delta, std, outlier_score, standardized = robust_fit(
                jacobian,
                residual,
            )
            estimate = parameters + delta
            estimate[2] = math.exp(estimate[2])
            rows.append(
                {
                    "model": model_index,
                    "outlier_fraction": fraction,
                    "h_true": true_params[0],
                    "beta_true": true_params[1],
                    "q_true": true_params[2],
                    "h_est": estimate[0],
                    "beta_est": estimate[1],
                    "q_est": estimate[2],
                    "h_std": std[0],
                    "beta_std": std[1],
                    "logq_std": std[2],
                    "outlier_score": outlier_score,
                    "residual_scale": float(
                        np.sqrt(np.mean(standardized**2))
                    ),
                }
            )
    return rows


def conformal_quantiles(rows, indices):
    z_h = np.abs([rows[index]["h_est"] - rows[index]["h_true"] for index in indices]) / np.maximum(
        np.asarray([rows[index]["h_std"] for index in indices]),
        1e-9,
    )
    z_beta = np.abs([rows[index]["beta_est"] - rows[index]["beta_true"] for index in indices]) / np.maximum(
        np.asarray([rows[index]["beta_std"] for index in indices]),
        1e-9,
    )
    z_q = np.abs(
        np.log(
            np.maximum(
                np.asarray([rows[index]["q_est"] for index in indices]),
                1e-9,
            )
            / np.asarray([rows[index]["q_true"] for index in indices])
        )
    ) / np.maximum(
        np.asarray([rows[index]["logq_std"] for index in indices]),
        1e-9,
    )
    return np.asarray(
        [
            np.quantile(z_h, 0.90, method="higher"),
            np.quantile(z_beta, 0.90, method="higher"),
            np.quantile(z_q, 0.90, method="higher"),
        ]
    )


def metrics(rows, indices, quantiles, outlier_score_bins=None, quantile_map=None):
    cover_h = []
    cover_beta = []
    cover_q = []
    width_h = []
    width_beta = []
    relative_width_q = []
    for index in indices:
        row = rows[index]
        q_values = quantiles
        if quantile_map is not None and outlier_score_bins is not None:
            q_values = quantile_map.get(int(outlier_score_bins[index]), quantiles)
        h_low = row["h_est"] - q_values[0] * row["h_std"]
        h_high = row["h_est"] + q_values[0] * row["h_std"]
        beta_low = row["beta_est"] - q_values[1] * row["beta_std"]
        beta_high = row["beta_est"] + q_values[1] * row["beta_std"]
        lower_logq = math.log(max(row["q_est"], 1e-9)) - q_values[2] * row["logq_std"]
        upper_logq = math.log(max(row["q_est"], 1e-9)) + q_values[2] * row["logq_std"]
        cover_h.append(h_low <= row["h_true"] <= h_high)
        cover_beta.append(beta_low <= row["beta_true"] <= beta_high)
        cover_q.append(lower_logq <= math.log(row["q_true"]) <= upper_logq)
        width_h.append(h_high - h_low)
        width_beta.append(beta_high - beta_low)
        relative_width_q.append(math.exp(upper_logq) / max(math.exp(lower_logq), 1e-9) - 1.0)
    return {
        "count": len(indices),
        "coverage_percent": [
            float(np.mean(cover_h) * 100.0),
            float(np.mean(cover_beta) * 100.0),
            float(np.mean(cover_q) * 100.0),
        ],
        "mean_width": [
            float(np.mean(width_h)),
            float(np.mean(width_beta)),
            float(np.mean(relative_width_q)),
        ],
        "mae": [
            float(np.mean([abs(rows[index]["h_est"] - rows[index]["h_true"]) for index in indices])),
            float(np.mean([abs(rows[index]["beta_est"] - rows[index]["beta_true"]) for index in indices])),
            float(np.mean([abs(rows[index]["q_est"] - rows[index]["q_true"]) / rows[index]["q_true"] for index in indices])),
        ],
    }


def main() -> None:
    calibration_rows = build_observation_set(
        np.random.default_rng(SEED + 100),
        N_CALIBRATION,
        100,
    )
    test_rows = build_observation_set(
        np.random.default_rng(SEED + 200),
        N_TEST,
        200,
    )
    all_calibration_indices = np.arange(len(calibration_rows))
    all_test_indices = np.arange(len(test_rows))
    raw_quantiles = np.asarray([Z90, Z90, Z90])
    global_quantiles = conformal_quantiles(
        calibration_rows,
        all_calibration_indices,
    )

    outlier_scores = np.asarray(
        [row["outlier_score"] for row in calibration_rows]
    )
    score_thresholds = np.quantile(outlier_scores, [1 / 3, 2 / 3])
    calibration_bins = np.digitize(outlier_scores, score_thresholds)

    test_outlier_scores = np.asarray(
        [row["outlier_score"] for row in test_rows]
    )
    test_bins = np.digitize(test_outlier_scores, score_thresholds)
    stratified_quantiles = {}
    for bin_index in (0, 1, 2):
        members = np.where(calibration_bins == bin_index)[0]
        stratified_quantiles[bin_index] = (
            conformal_quantiles(calibration_rows, members)
            if len(members) >= 8
            else global_quantiles
        )

    calibration_h_median = float(
        np.median([row["h_true"] for row in calibration_rows])
    )
    calibration_beta_median = float(
        np.median([row["beta_true"] for row in calibration_rows])
    )
    calibration_score_median = float(
        np.median([row["outlier_score"] for row in calibration_rows])
    )
    cross_domain = {}
    scenarios = {
        "thin_to_thick": (
            lambda row: row["h_true"] < calibration_h_median,
            lambda row: row["h_true"] >= calibration_h_median,
        ),
        "thick_to_thin": (
            lambda row: row["h_true"] >= calibration_h_median,
            lambda row: row["h_true"] < calibration_h_median,
        ),
        "low_beta_to_high_beta": (
            lambda row: row["beta_true"] < calibration_beta_median,
            lambda row: row["beta_true"] >= calibration_beta_median,
        ),
        "high_beta_to_low_beta": (
            lambda row: row["beta_true"] >= calibration_beta_median,
            lambda row: row["beta_true"] < calibration_beta_median,
        ),
        "low_outlier_to_high_outlier": (
            lambda row: row["outlier_score"] < calibration_score_median,
            lambda row: row["outlier_score"] >= calibration_score_median,
        ),
        "high_outlier_to_low_outlier": (
            lambda row: row["outlier_score"] >= calibration_score_median,
            lambda row: row["outlier_score"] < calibration_score_median,
        ),
    }
    for name, (calibration_rule, test_rule) in scenarios.items():
        calibration_indices = [
            index
            for index, row in enumerate(calibration_rows)
            if calibration_rule(row)
        ]
        test_indices = [
            index
            for index, row in enumerate(test_rows)
            if test_rule(row)
        ]
        if len(calibration_indices) < 12 or len(test_indices) < 12:
            cross_domain[name] = {"status": "insufficient_samples"}
            continue
        quantiles = conformal_quantiles(
            calibration_rows,
            calibration_indices,
        )
        cross_domain[name] = {
            "n_calibration": len(calibration_indices),
            "n_test": len(test_indices),
            "quantiles": quantiles.tolist(),
            "raw": metrics(test_rows, test_indices, raw_quantiles),
            "conformal": metrics(test_rows, test_indices, quantiles),
        }

    summary = {
        "settings": {
            "n_calibration": len(calibration_rows),
            "n_test": len(test_rows),
            "outlier_fractions": OUTLIER_FRACTIONS,
            "raw_interval_z": Z90,
            "global_quantiles": global_quantiles.tolist(),
            "score_thresholds": score_thresholds.tolist(),
            "stratified_quantiles": {
                str(key): value.tolist()
                for key, value in stratified_quantiles.items()
            },
        },
        "overall": {
            "raw": metrics(
                test_rows,
                all_test_indices,
                raw_quantiles,
            ),
            "global_conformal": metrics(
                test_rows,
                all_test_indices,
                global_quantiles,
            ),
            "outlier_stratified": metrics(
                test_rows,
                all_test_indices,
                global_quantiles,
                outlier_score_bins=test_bins,
                quantile_map=stratified_quantiles,
            ),
        },
        "cross_domain": cross_domain,
        "by_outlier_fraction": {},
    }
    for fraction in OUTLIER_FRACTIONS:
        selection = np.where(
            [row["outlier_fraction"] == fraction for row in test_rows]
        )[0]
        summary["by_outlier_fraction"][str(fraction)] = {
            "raw": metrics(test_rows, selection, raw_quantiles),
            "global_conformal": metrics(
                test_rows,
                selection,
                global_quantiles,
            ),
            "outlier_stratified": metrics(
                test_rows,
                selection,
                global_quantiles,
                outlier_score_bins=test_bins,
                quantile_map=stratified_quantiles,
            ),
        }

    (ROOT / "e18_robust_posterior_calibration_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    combined_rows = [
        {"split": "calibration", **row}
        for row in calibration_rows
    ] + [
        {"split": "test", **row}
        for row in test_rows
    ]
    with (ROOT / "e18_robust_posterior_calibration_runs.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(combined_rows[0].keys()))
        writer.writeheader()
        writer.writerows(combined_rows)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
