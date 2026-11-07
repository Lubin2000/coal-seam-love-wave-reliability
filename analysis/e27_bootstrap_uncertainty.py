"""E27: model-level bootstrap confidence intervals for headline metrics."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parent
N_BOOTSTRAP = 2000
SEED = 20261106
Z90 = 1.6448536269514722


def load_csv(name):
    with (ROOT / name).open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def percentile_interval(values):
    return [
        float(np.quantile(values, 0.025)),
        float(np.quantile(values, 0.975)),
    ]


def bootstrap_picking():
    rows = load_csv("e26_additional_picker_runs.csv")
    models = sorted({int(row["model"]) for row in rows})
    rng = np.random.default_rng(SEED)
    methods = (
        "independent_peaks",
        "greedy_tracking",
        "viterbi_order_tracking",
        "component_segmentation",
    )

    def aggregate(sample_rows):
        result = {}
        for method in methods:
            subset = [
                row for row in sample_rows if row["method"] == method
            ]
            false_picks = sum(int(row["false_picks"]) for row in subset)
            total_picks = sum(int(row["total_picks"]) for row in subset)
            result[method] = {
                "mode_detection_percent": {
                    f"mode{mode}": 100.0
                    * sum(int(row[f"mode{mode}_detected"]) for row in subset)
                    / max(
                        sum(int(row[f"mode{mode}_possible"]) for row in subset),
                        1,
                    )
                    for mode in range(3)
                },
                "false_pick_percent": 100.0
                * false_picks
                / max(total_picks, 1),
                "median_h_width_90": float(
                    np.median(
                        [float(row["h_width_90"]) for row in subset]
                    )
                ),
            }
        return result

    rows_by_model = {
        model: [row for row in rows if int(row["model"]) == model]
        for model in models
    }
    point = aggregate(rows)
    samples = []
    for _ in range(N_BOOTSTRAP):
        sample_models = rng.choice(models, size=len(models), replace=True)
        sample_rows = [
            row
            for model in sample_models
            for row in rows_by_model[int(model)]
        ]
        samples.append(aggregate(sample_rows))

    intervals = point
    for method in methods:
        intervals[method]["mode_detection_percent_ci"] = {}
        for mode in range(3):
            key = f"mode{mode}"
            intervals[method]["mode_detection_percent_ci"][key] = (
                percentile_interval(
                    [
                        sample[method]["mode_detection_percent"][key]
                        for sample in samples
                    ]
                )
            )
        for key in ("false_pick_percent", "median_h_width_90"):
            intervals[method][f"{key}_ci"] = percentile_interval(
                [sample[method][key] for sample in samples]
            )
    return intervals


def bootstrap_inversion():
    rows = load_csv("e17_robust_likelihood_runs.csv")
    models = sorted({int(row["model"]) for row in rows})
    rng = np.random.default_rng(SEED + 1)

    result = {}
    for fraction in (0.0, 0.2, 0.3):
        for include_q in (False, True):
            for likelihood in ("gaussian", "huber", "student_t"):
                subset = [
                    row
                    for row in rows
                    if float(row["outlier_fraction"]) == fraction
                    and row["include_attenuation"] == str(include_q)
                    and row["likelihood"] == likelihood
                ]
                if not subset:
                    continue

                def aggregate(sample_rows):
                    h_errors = np.asarray(
                        [float(row["h_error"]) for row in sample_rows]
                    )
                    beta_errors = np.asarray(
                        [float(row["beta_error"]) for row in sample_rows]
                    )
                    q_errors = np.asarray(
                        [float(row["q_relative_error"]) for row in sample_rows]
                    )
                    return {
                        "h_rmse": float(np.sqrt(np.mean(h_errors**2))),
                        "beta_rmse": float(np.sqrt(np.mean(beta_errors**2))),
                        "h_coverage_percent": float(
                            np.mean(
                                np.abs(h_errors)
                                <= Z90
                                * np.asarray(
                                    [
                                        float(row["h_std"])
                                        for row in sample_rows
                                    ]
                                )
                            )
                            * 100.0
                        ),
                        "q_relative_mae": (
                            float(np.mean(np.abs(q_errors)))
                            if include_q
                            else None
                        ),
                    }

                point = aggregate(subset)
                subset_by_model = {
                    model: [
                        row
                        for row in subset
                        if int(row["model"]) == model
                    ]
                    for model in models
                }
                sample_values = {"h_rmse": [], "beta_rmse": [], "h_coverage_percent": []}
                if include_q:
                    sample_values["q_relative_mae"] = []
                for _ in range(N_BOOTSTRAP):
                    sampled_models = rng.choice(
                        models,
                        size=len(models),
                        replace=True,
                    )
                    sample_rows = [
                        row
                        for model in sampled_models
                        for row in subset_by_model[int(model)]
                    ]
                    if not sample_rows:
                        continue
                    values = aggregate(sample_rows)
                    for key in sample_values:
                        sample_values[key].append(values[key])
                key = f"fraction_{fraction}_attenuation_{include_q}_{likelihood}"
                result[key] = {
                    **point,
                    **{
                        f"{name}_ci": percentile_interval(values)
                        for name, values in sample_values.items()
                    },
                }
    return result


def bootstrap_calibration():
    rows = [
        row
        for row in load_csv("e18_robust_posterior_calibration_runs.csv")
        if row["split"] == "test"
    ]
    summary = json.loads(
        (
            ROOT / "e18_robust_posterior_calibration_summary.json"
        ).read_text(encoding="utf-8")
    )
    quantiles = np.asarray(summary["settings"]["global_quantiles"])
    models = sorted({int(row["model"]) for row in rows})
    rng = np.random.default_rng(SEED + 2)

    def aggregate(sample_rows):
        h_covered = []
        beta_covered = []
        q_covered = []
        for row in sample_rows:
            h_covered.append(
                abs(float(row["h_est"]) - float(row["h_true"]))
                <= quantiles[0] * float(row["h_std"])
            )
            beta_covered.append(
                abs(float(row["beta_est"]) - float(row["beta_true"]))
                <= quantiles[1] * float(row["beta_std"])
            )
            q_covered.append(
                abs(
                    math.log(float(row["q_est"]))
                    - math.log(float(row["q_true"]))
                )
                <= quantiles[2] * float(row["logq_std"])
            )
        return [
            100.0 * float(np.mean(h_covered)),
            100.0 * float(np.mean(beta_covered)),
            100.0 * float(np.mean(q_covered)),
        ]

    rows_by_model = {
        model: [row for row in rows if int(row["model"]) == model]
        for model in models
    }
    point = aggregate(rows)
    samples = []
    for _ in range(N_BOOTSTRAP):
        sampled_models = rng.choice(models, size=len(models), replace=True)
        sample_rows = [
            row
            for model in sampled_models
            for row in rows_by_model[int(model)]
        ]
        samples.append(aggregate(sample_rows))
    samples = np.asarray(samples)
    return {
        "coverage_percent": point,
        "coverage_percent_ci": {
            "thickness": percentile_interval(samples[:, 0]),
            "velocity": percentile_interval(samples[:, 1]),
            "q": percentile_interval(samples[:, 2]),
        },
    }


def main() -> None:
    output = {
        "settings": {
            "n_bootstrap": N_BOOTSTRAP,
            "resampling_unit": "model",
            "confidence_level": 0.95,
        },
        "picking": bootstrap_picking(),
        "inversion": bootstrap_inversion(),
        "calibration": bootstrap_calibration(),
    }
    (ROOT / "e27_bootstrap_uncertainty_summary.json").write_text(
        json.dumps(output, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
