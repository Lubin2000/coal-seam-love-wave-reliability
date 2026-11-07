"""E17: robust local inversion of automatically tracked mode picks."""

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


SEED = 20261105
N_MODELS = 30
PRIOR_STD = np.asarray([0.60, 100.0, 0.65])
PHASE_SIGMA = 0.012
ATTENUATION_SIGMA = 0.10
Z90 = e13.Z90


def build_observation_rows(center, picks, rng):
    rows = []
    for frequency_index, picked_velocities in enumerate(picks):
        for assigned_mode, picked_velocity in enumerate(picked_velocities[:3]):
            if not math.isfinite(picked_velocity):
                continue
            true_velocity = center["phase"][frequency_index, assigned_mode]
            true_attenuation = center["attenuation"][frequency_index, assigned_mode]
            if not np.isfinite(true_velocity) or not np.isfinite(true_attenuation):
                continue
            q_observed = true_attenuation * math.exp(
                rng.normal(0.0, ATTENUATION_SIGMA)
            )
            rows.append(
                {
                    "frequency_index": frequency_index,
                    "mode": assigned_mode,
                    "phase_observed": picked_velocity,
                    "phase_truth": true_velocity,
                    "phase_sigma": max(PHASE_SIGMA * true_velocity, 1.0),
                    "attenuation_truth": true_attenuation,
                    "q_observed": q_observed,
                    "q_sigma": max(ATTENUATION_SIGMA * true_attenuation, 1e-6),
                }
            )
    return rows


def build_design_matrices(
    center,
    plus_responses,
    minus_responses,
    observations,
    include_q: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    jacobian_phase = []
    residual_phase = []
    sigmas_phase = []
    jacobian_q = []
    residual_q = []
    sigmas_q = []
    for observation in observations:
        frequency_index = observation["frequency_index"]
        mode = observation["mode"]
        phase_values = []
        for parameter_index, step in enumerate((0.004, 1.0, 0.008)):
            plus = plus_responses[parameter_index]["phase"][frequency_index, mode]
            minus = minus_responses[parameter_index]["phase"][frequency_index, mode]
            phase_values.append((plus - minus) / (2.0 * step))
        q_values = None
        if include_q:
            q_values = []
            for parameter_index, step in enumerate((0.004, 1.0, 0.008)):
                plus = plus_responses[parameter_index]["attenuation"][
                    frequency_index, mode
                ]
                minus = minus_responses[parameter_index]["attenuation"][
                    frequency_index, mode
                ]
                q_values.append((plus - minus) / (2.0 * step))
        if not np.all(np.isfinite(phase_values)):
            continue
        if include_q and not np.all(np.isfinite(q_values)):
            continue
        jacobian_phase.append(phase_values)
        residual_phase.append(
            observation["phase_truth"] - observation["phase_observed"]
        )
        sigmas_phase.append(observation["phase_sigma"])
        if include_q:
            jacobian_q.append(q_values)
            residual_q.append(
                observation["attenuation_truth"] - observation["q_observed"]
            )
            sigmas_q.append(observation["q_sigma"])
    if len(jacobian_phase) == 0 or (include_q and len(jacobian_q) == 0):
        return (
            np.zeros((0, 3)),
            np.zeros(0),
            np.zeros(0),
        )
    if include_q:
        jacobian = np.vstack(
            (
                np.asarray(jacobian_phase) / np.asarray(sigmas_phase)[:, None],
                np.asarray(jacobian_q) / np.asarray(sigmas_q)[:, None],
            )
        )
        residual = np.concatenate(
            (
                np.asarray(residual_phase) / np.asarray(sigmas_phase),
                np.asarray(residual_q) / np.asarray(sigmas_q),
            )
        )
    else:
        jacobian = np.asarray(jacobian_phase) / np.asarray(sigmas_phase)[:, None]
        residual = np.asarray(residual_phase) / np.asarray(sigmas_phase)
    jacobian = np.nan_to_num(jacobian, nan=0.0, posinf=0.0, neginf=0.0)
    residual = np.nan_to_num(residual, nan=0.0, posinf=0.0, neginf=0.0)
    return jacobian, residual, np.asarray(sigmas_phase)


def robust_update(jacobian, residual, likelihood):
    delta = np.zeros(3)
    weights = np.ones(len(residual))
    for _ in range(30):
        standardized = residual + jacobian @ delta
        if likelihood == "gaussian":
            weights = np.ones(len(standardized))
        elif likelihood == "huber":
            threshold = 1.5
            weights = np.minimum(
                1.0,
                threshold / np.maximum(np.abs(standardized), 1e-8),
            )
        elif likelihood == "student_t":
            nu = 4.0
            weights = (nu + 1.0) / (
                nu + standardized**2
            )
        else:
            raise ValueError(likelihood)
        weighted_jacobian = jacobian * np.sqrt(weights[:, None])
        information = weighted_jacobian.T @ weighted_jacobian + np.diag(
            1.0 / PRIOR_STD**2
        )
        rhs = -weighted_jacobian.T @ (
            np.sqrt(weights) * residual
        )
        new_delta = np.linalg.solve(information, rhs)
        if np.linalg.norm(new_delta - delta) < 1e-9:
            delta = new_delta
            break
        delta = new_delta
    standardized = residual + jacobian @ delta
    weighted_jacobian = jacobian * np.sqrt(weights[:, None])
    covariance = np.linalg.pinv(
        weighted_jacobian.T @ weighted_jacobian
        + np.diag(1.0 / PRIOR_STD**2)
    )
    return delta, np.sqrt(np.maximum(np.diag(covariance), 0.0))


def inject_outliers(rng, observations, fraction):
    perturbed = [dict(observation) for observation in observations]
    count = int(round(fraction * len(perturbed)))
    if count == 0:
        return perturbed
    indices = rng.choice(len(perturbed), size=count, replace=False)
    for index in indices:
        observation = perturbed[index]
        sign = rng.choice([-1.0, 1.0])
        observation["phase_observed"] += sign * rng.uniform(80.0, 250.0)
        observation["q_observed"] *= math.exp(sign * rng.uniform(0.4, 1.0))
    return perturbed


def main() -> None:
    rng = np.random.default_rng(SEED)
    rows = []
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

        image = e15.make_dispersion_image(center["phase"], rng)
        candidates = e16.detect_candidate_peaks(image)
        picks = e16.viterbi_track(candidates)
        observations = build_observation_rows(center, picks, rng)
        if len(observations) < 6:
            continue
        for fraction in (0.0, 0.10, 0.20, 0.30):
            perturbed = inject_outliers(rng, observations, fraction)
            for include_q in (False, True):
                jacobian, residual, _ = build_design_matrices(
                    center,
                    plus_responses,
                    minus_responses,
                    perturbed,
                    include_q,
                )
                if len(residual) < 3:
                    continue
                for likelihood in ("gaussian", "huber", "student_t"):
                    delta, std = robust_update(
                        jacobian,
                        residual,
                        likelihood,
                    )
                    estimated = parameters + delta
                    estimated[2] = math.exp(estimated[2])
                    true_physical = np.asarray([h, beta, q])
                    error = estimated - true_physical
                    row = {
                        "model": model_index,
                        "outlier_fraction": fraction,
                        "include_attenuation": include_q,
                        "likelihood": likelihood,
                        "h_error": float(error[0]),
                        "beta_error": float(error[1]),
                        "q_relative_error": float(abs(error[2]) / q),
                        "h_std": float(std[0]),
                        "beta_std": float(std[1]),
                        "q_std_relative": float(std[2]),
                    }
                    rows.append(row)
        if (model_index + 1) % 10 == 0:
            print(f"processed {model_index + 1}/{N_MODELS}")

    summary = []
    for fraction in (0.0, 0.10, 0.20, 0.30):
        for include_q in (False, True):
            for likelihood in ("gaussian", "huber", "student_t"):
                subset = [
                    row
                    for row in rows
                    if row["outlier_fraction"] == fraction
                    and row["include_attenuation"] == include_q
                    and row["likelihood"] == likelihood
                ]
                summary.append(
                    {
                        "outlier_fraction": fraction,
                        "include_attenuation": include_q,
                        "likelihood": likelihood,
                        "h_rmse": float(
                            np.sqrt(np.mean([row["h_error"] ** 2 for row in subset]))
                        ),
                        "h_bias": float(
                            np.mean([row["h_error"] for row in subset])
                        ),
                        "beta_rmse": float(
                            np.sqrt(np.mean([row["beta_error"] ** 2 for row in subset]))
                        ),
                        "q_relative_mae": float(
                            np.mean(
                                [abs(row["q_relative_error"]) for row in subset]
                            )
                        )
                        if include_q
                        else None,
                        "h_coverage_percent": float(
                            np.mean(
                                [
                                    abs(row["h_error"]) <= Z90 * row["h_std"]
                                    for row in subset
                                ]
                            )
                            * 100.0
                        ),
                    }
                )

    with (ROOT / "e17_robust_likelihood_runs.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (ROOT / "e17_robust_likelihood_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
