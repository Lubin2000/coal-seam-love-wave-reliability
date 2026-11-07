"""E24: track and invert dispersion images extracted from independent SH fields.

The same order-constrained tracker is applied to:
1. analytical synthetic images with known modal ridges;
2. phase-shift images extracted from the E22 SH finite-difference fields.
Thickness and shear velocity are then estimated from the tracked picks.
"""

from __future__ import annotations

import csv
from itertools import combinations, permutations
import json
import math
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import e13_acquisition_information as e13


FREQUENCIES = np.arange(100.0, 561.0, 20.0)
CASES = (
    {"name": "thin_1p8m", "h": 1.8, "beta": 1120.0},
    {"name": "medium_2p8m", "h": 2.8, "beta": 1120.0},
    {"name": "thick_3p8m", "h": 3.8, "beta": 1120.0},
)
PEAK_THRESHOLD = 0.08
START_AMPLITUDE_THRESHOLD = 0.14
MAX_JUMP = 120.0
NU = 4.0


def true_modes(h: float, beta: float) -> np.ndarray:
    modes = np.full((len(FREQUENCIES), 3), np.nan)
    for frequency_index, frequency in enumerate(FREQUENCIES):
        values = e13.all_phase_velocities(
            float(frequency),
            h,
            beta,
        )
        modes[frequency_index, : min(3, len(values))] = values[:3]
    return modes


def analytical_image(
    modes: np.ndarray,
    velocities: np.ndarray,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    image = np.zeros((len(FREQUENCIES), len(velocities)))
    amplitudes = (1.0, 0.58, 0.33)
    widths = (24.0, 28.0, 34.0)
    for frequency_index in range(len(FREQUENCIES)):
        visibility = rng.uniform(0.75, 1.0)
        if rng.random() < 0.10:
            visibility *= 0.4
        for mode_index in range(3):
            velocity = modes[frequency_index, mode_index]
            if not np.isfinite(velocity):
                continue
            image[frequency_index] += (
                amplitudes[mode_index]
                * visibility
                * np.exp(
                    -0.5
                    * ((velocities - velocity) / widths[mode_index]) ** 2
                )
            )
        image[frequency_index] += 0.20 * np.exp(
            -0.5
            * (
                (
                    velocities
                    - (1450.0 + 60.0 * math.sin(frequency_index / 8.0))
                )
                / 30.0
            )
            ** 2
        )
        image[frequency_index] += rng.normal(
            0.0,
            0.025,
            len(velocities),
        )
    image = np.maximum(image, 0.0)
    image /= np.maximum(
        np.max(image, axis=1, keepdims=True),
        1e-12,
    )
    return image


def candidate_rows(
    image: np.ndarray,
    velocities: np.ndarray,
) -> list[list[tuple[float, float]]]:
    normalised = image / np.maximum(
        np.max(image, axis=1, keepdims=True),
        1e-12,
    )
    rows = []
    for row in normalised:
        candidates = []
        for index in range(1, len(row) - 1):
            if (
                row[index] >= row[index - 1]
                and row[index] >= row[index + 1]
                and row[index] >= PEAK_THRESHOLD
            ):
                candidates.append(
                    (float(velocities[index]), float(row[index]))
                )
        candidates.sort(key=lambda item: item[1], reverse=True)
        separated = []
        for velocity, amplitude in candidates:
            if all(
                abs(velocity - previous_velocity) > 45.0
                for previous_velocity, _ in separated
            ):
                separated.append((velocity, amplitude))
            if len(separated) == 6:
                break
        rows.append(separated)
    return rows


def independent_peaks(
    candidates: list[list[tuple[float, float]]],
) -> list[list[float]]:
    picks = []
    for candidates_at_frequency in candidates:
        selected = []
        for velocity, amplitude in candidates_at_frequency:
            if amplitude < 0.12:
                continue
            if all(abs(velocity - value) > 40.0 for value in selected):
                selected.append(velocity)
            if len(selected) == 3:
                break
        picks.append(sorted(selected, reverse=True))
    return picks


def greedy_picks(
    candidates: list[list[tuple[float, float]]],
) -> list[list[float]]:
    tracks = [math.nan] * 3
    picks = [[] for _ in candidates]
    start = None
    for frequency_index in range(len(candidates) - 1, -1, -1):
        strong = [
            item
            for item in candidates[frequency_index]
            if item[1] >= START_AMPLITUDE_THRESHOLD
        ]
        if strong:
            start = frequency_index
            seed = sorted(
                [item[0] for item in strong[:3]],
                reverse=True,
            )
            for mode_index, velocity in enumerate(seed):
                tracks[mode_index] = velocity
                picks[frequency_index].append(velocity)
            break
    if start is None:
        return picks

    for frequency_index in range(start - 1, -1, -1):
        used = set()
        selected = []
        for mode_index, previous in enumerate(tracks):
            if not math.isfinite(previous):
                continue
            viable = []
            for candidate_index, (velocity, amplitude) in enumerate(
                candidates[frequency_index]
            ):
                if candidate_index in used:
                    continue
                jump = velocity - previous
                if -10.0 <= jump <= MAX_JUMP:
                    score = amplitude - 0.45 * (abs(jump) / MAX_JUMP) ** 2
                    viable.append((score, candidate_index, velocity))
            if viable:
                _, candidate_index, velocity = max(viable)
                used.add(candidate_index)
                tracks[mode_index] = velocity
                selected.append(velocity)
        picks[frequency_index] = sorted(selected, reverse=True)
    return picks


def states_at_frequency(
    candidates: list[tuple[float, float]],
) -> list[tuple[tuple[float | None, ...], float]]:
    states = []
    for mode_count in range(min(3, len(candidates)) + 1):
        for selected in combinations(candidates, mode_count):
            for active_modes in combinations(range(3), mode_count):
                sorted_velocities = sorted(
                    (item[0] for item in selected),
                )
                slots: list[float | None] = [None, None, None]
                for mode, velocity in zip(
                    active_modes,
                    sorted_velocities,
                    strict=True,
                ):
                    slots[mode] = velocity
                emission = (
                    sum(item[1] for item in selected) - 0.08 * mode_count
                )
                states.append((tuple(slots), emission))
    unique = {}
    for velocities, emission in states:
        unique[velocities] = max(
            emission,
            unique.get(velocities, -math.inf),
        )
    return [(velocities, emission) for velocities, emission in unique.items()]


def global_tracking(
    candidates: list[list[tuple[float, float]]],
) -> list[list[float]]:
    states = [states_at_frequency(row) for row in candidates]
    scores = {
        velocities: emission
        for velocities, emission in states[-1]
    }
    predecessors = [None] * len(states)
    predecessors[-1] = {}

    for frequency_index in range(len(states) - 2, -1, -1):
        next_scores = {}
        next_predecessors = {}
        for current, emission in states[frequency_index]:
            best_score = -math.inf
            best_next = ()
            for following, following_score in scores.items():
                smoothness = sum(
                    (
                        (current[index] - following[index]) / 100.0
                    )
                    ** 2
                    for index in range(3)
                    if current[index] is not None
                    and following[index] is not None
                )
                appearance = (
                    sum(
                        current[index] is not None
                        and following[index] is None
                        for index in range(3)
                    )
                    * 0.5
                )
                score = (
                    following_score
                    + emission
                    - 0.20 * smoothness
                    - appearance
                )
                if score > best_score:
                    best_score = score
                    best_next = following
            next_scores[current] = best_score
            next_predecessors[current] = best_next
        scores = next_scores
        predecessors[frequency_index] = next_predecessors

    state = max(scores, key=scores.get)
    tracked = []
    for frequency_index in range(len(states)):
        tracked.append(
            [
                float(value) if value is not None else math.nan
                for value in state
            ]
        )
        if frequency_index < len(states) - 1:
            state = predecessors[frequency_index][state]
    return tracked


def tolerance(frequency: float, velocity: float) -> float:
    aperture = 198.0
    return max(50.0, 1.5 * velocity**2 / (frequency * aperture))


def optimal_matching(
    mode_values: list[float],
    picks: list[float],
    frequency: float,
) -> tuple[dict[int, tuple[int, float]], int]:
    """One-to-one matching with unmatched modes and false picks allowed."""
    best_key = (-1, -math.inf)
    best_matches = {}
    choices = list(range(-1, len(picks)))
    for assignment in permutations(choices, len(mode_values)):
        if len(set(value for value in assignment if value >= 0)) < len(
            [value for value in assignment if value >= 0]
        ):
            continue
        matches = {}
        error_sum = 0.0
        valid = True
        for mode_index, (mode_value, pick_index) in enumerate(
            zip(mode_values, assignment)
        ):
            if pick_index < 0:
                continue
            error = abs(picks[pick_index] - mode_value)
            if error > tolerance(frequency, mode_value):
                valid = False
                break
            matches[mode_index] = (pick_index, error)
            error_sum += error
        if not valid:
            continue
        key = (len(matches), -error_sum)
        if key > best_key:
            best_key = key
            best_matches = matches
    matched_pick_indices = {
        pick_index for pick_index, _ in best_matches.values()
    }
    return best_matches, len(picks) - len(matched_pick_indices)


def tracking_statistics(
    modes: np.ndarray,
    picks: list[list[float]],
) -> dict:
    detected = {mode: 0 for mode in range(3)}
    possible = {mode: 0 for mode in range(3)}
    errors = {mode: [] for mode in range(3)}
    false_picks = 0
    total_picks = 0
    for frequency_index, frequency in enumerate(FREQUENCIES):
        valid_modes = [
            (mode, value)
            for mode, value in enumerate(modes[frequency_index])
            if math.isfinite(value)
        ]
        current_picks = picks[frequency_index]
        if len(current_picks) == 3 and any(
            not math.isfinite(value) for value in current_picks
        ):
            for original_mode, mode_value in valid_modes:
                possible[original_mode] += 1
                picked = current_picks[original_mode]
                if not math.isfinite(picked):
                    continue
                error = abs(picked - mode_value)
                if error <= tolerance(frequency, mode_value):
                    detected[original_mode] += 1
                    errors[original_mode].append(error)
                else:
                    false_picks += 1
                total_picks += 1
        else:
            mode_values = [value for _, value in valid_modes]
            matches, false_count = optimal_matching(
                mode_values,
                current_picks,
                frequency,
            )
            false_picks += false_count
            total_picks += len(current_picks)
            for match_mode_index, mode_value in enumerate(mode_values):
                original_mode = valid_modes[match_mode_index][0]
                possible[original_mode] += 1
                if match_mode_index in matches:
                    _, error = matches[match_mode_index]
                    detected[original_mode] += 1
                    errors[original_mode].append(error)
    return {
        "mode_detection_percent": {
            f"mode{mode}": 100.0
            * detected[mode]
            / max(possible[mode], 1)
            for mode in range(3)
        },
        "mode_detected": {f"mode{mode}": detected[mode] for mode in range(3)},
        "mode_possible": {f"mode{mode}": possible[mode] for mode in range(3)},
        "median_pick_error_mps": {
            f"mode{mode}": (
                float(np.median(errors[mode]))
                if errors[mode]
                else None
            )
            for mode in range(3)
        },
        "false_pick_percent": 100.0 * false_picks / max(total_picks, 1),
    }


def model_prediction(h: float, beta: float):
    return true_modes(h, beta)


def robust_objective(
    h: float,
    beta: float,
    observations: list[tuple[float, int, float]],
) -> float:
    if not (0.8 <= h <= 5.5 and 800.0 <= beta <= 1450.0):
        return 1.0e6
    predicted = model_prediction(h, beta)
    residuals = []
    for frequency, mode, velocity in observations:
        frequency_index = int(round((frequency - 100.0) / 20.0))
        value = predicted[frequency_index, mode]
        if not math.isfinite(value):
            residuals.append(5.0)
            continue
        sigma = max(0.02 * velocity, 20.0)
        residuals.append((value - velocity) / sigma)
    if not residuals:
        return 1.0e6
    residuals = np.asarray(residuals)
    robust = 0.5 * (NU + 1.0) * np.log1p(residuals**2 / NU)
    prior = ((h - 2.5) / 1.2) ** 2 + ((beta - 1100.0) / 180.0) ** 2
    return float(np.mean(robust) + 0.02 * prior)


def observations_from_tracks(
    picks: list[list[float]],
) -> list[tuple[float, int, float]]:
    observations = []
    for frequency_index, values in enumerate(picks):
        for mode, velocity in enumerate(values[:3]):
            if not math.isfinite(velocity):
                continue
            observations.append(
                (float(FREQUENCIES[frequency_index]), mode, velocity)
            )
    return observations


def invert(observations: list[tuple[float, int, float]]) -> tuple[float, float]:
    best = (math.inf, 2.5, 1100.0)
    for h in np.arange(1.0, 5.51, 0.45):
        for beta in np.arange(850.0, 1401.0, 55.0):
            score = robust_objective(float(h), float(beta), observations)
            if score < best[0]:
                best = (score, float(h), float(beta))

    _, h, beta = best
    h_step = 0.25
    beta_step = 30.0
    score = best[0]
    for _ in range(50):
        improved = False
        for candidate_h, candidate_beta in (
            (h + h_step, beta),
            (h - h_step, beta),
            (h, beta + beta_step),
            (h, beta - beta_step),
        ):
            candidate_score = robust_objective(
                candidate_h,
                candidate_beta,
                observations,
            )
            if candidate_score < score:
                h, beta, score = (
                    candidate_h,
                    candidate_beta,
                    candidate_score,
                )
                improved = True
        if not improved:
            h_step *= 0.55
            beta_step *= 0.55
            if h_step < 0.002 and beta_step < 0.3:
                break
    return h, beta


def main() -> None:
    archive = np.load(ROOT / "e22_phase_shift_images.npz")
    velocities = archive["velocities"]
    rows = []
    for case_index, case in enumerate(CASES):
        modes = true_modes(case["h"], case["beta"])
        images = {
            "analytical_image": analytical_image(
                modes,
                velocities,
                20261105 + case_index,
            ),
            "finite_difference_image": archive[case["name"]],
        }
        for image_name, image in images.items():
            candidates = candidate_rows(image, velocities)
            tracked = {
                "independent_peaks": independent_peaks(candidates),
                "greedy_tracking": greedy_picks(candidates),
                "global_tracking": global_tracking(candidates),
            }
            for method, picks in tracked.items():
                stats = tracking_statistics(modes, picks)
                observations = observations_from_tracks(picks)
                if method == "global_tracking":
                    h_est, beta_est = invert(observations)
                else:
                    h_est, beta_est = math.nan, math.nan
                rows.append(
                    {
                        "case": case["name"],
                        "h_true": case["h"],
                        "beta_true": case["beta"],
                        "image": image_name,
                        "method": method,
                        "observation_count": len(observations),
                        "h_est": h_est,
                        "beta_est": beta_est,
                        "h_error": h_est - case["h"],
                        "beta_error": beta_est - case["beta"],
                        **stats,
                    }
                )
        print(f"processed {case['name']}", flush=True)

    with (ROOT / "e24_wavefield_tracking_inversion.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        flat_rows = []
        for row in rows:
            flat = {
                key: value
                for key, value in row.items()
                if not isinstance(value, dict)
            }
            for mode in range(3):
                for key in (
                    "mode_detection_percent",
                    "mode_detected",
                    "mode_possible",
                    "median_pick_error_mps",
                ):
                    flat[f"{key}_mode{mode}"] = row[key][f"mode{mode}"]
            flat_rows.append(flat)
        writer = csv.DictWriter(handle, fieldnames=list(flat_rows[0]))
        writer.writeheader()
        writer.writerows(flat_rows)

    summary = {}
    for image_name in ("analytical_image", "finite_difference_image"):
        summary[image_name] = {}
        for method in (
            "independent_peaks",
            "greedy_tracking",
            "global_tracking",
        ):
            subset = [
                row
                for row in rows
                if row["image"] == image_name
                and row["method"] == method
            ]
            inversion_subset = [
                row for row in subset if math.isfinite(row["h_est"])
            ]
            inversion_summary = (
                {
                    "h_mae": float(
                        np.mean(
                            [abs(row["h_error"]) for row in inversion_subset]
                        )
                    ),
                    "beta_mae": float(
                        np.mean(
                            [
                                abs(row["beta_error"])
                                for row in inversion_subset
                            ]
                        )
                    ),
                }
                if inversion_subset
                else {"h_mae": None, "beta_mae": None}
            )
            summary[image_name][method] = {
                **inversion_summary,
                "mode_detection_percent": {
                    f"mode{mode}": float(
                        np.mean(
                            [
                                row["mode_detection_percent"][f"mode{mode}"]
                                for row in subset
                            ]
                        )
                    )
                    for mode in range(3)
                },
                "false_pick_percent": float(
                    np.mean(
                        [row["false_pick_percent"] for row in subset]
                    )
                ),
            }
    payload = {"settings": {"frequencies_hz": FREQUENCIES.tolist()}, "cases": rows, "summary": summary}
    (ROOT / "e24_wavefield_tracking_inversion_summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
