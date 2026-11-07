"""E16: continuity-constrained tracking of multi-mode dispersion ridges."""

from __future__ import annotations

import csv
import json
import math
from itertools import combinations
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import e13_acquisition_information as e13
import e15_automatic_mode_picking as e15


SEED = 20261105
N_MODELS = 50
PEAK_THRESHOLD = 0.10
START_AMPLITUDE_THRESHOLD = 0.18
MAX_JUMP = 95.0
Z90 = e13.Z90


def detect_candidate_peaks(image: np.ndarray):
    candidates = []
    for frequency_index in range(image.shape[0]):
        row = image[frequency_index]
        frequency_candidates = []
        for index in range(1, len(row) - 1):
            if row[index] >= row[index - 1] and row[index] >= row[index + 1]:
                if row[index] >= PEAK_THRESHOLD:
                    frequency_candidates.append(
                        (float(e15.VELOCITY_GRID[index]), float(row[index]))
                    )
        candidates.append(frequency_candidates)
    return candidates


def track_modes(candidate_rows):
    tracks = [np.nan, np.nan, np.nan]
    picks = [[ ] for _ in range(len(candidate_rows))]
    start_frequency = None
    for frequency_index in range(len(candidate_rows) - 1, -1, -1):
        strong = [
            item
            for item in candidate_rows[frequency_index]
            if item[1] >= START_AMPLITUDE_THRESHOLD
        ]
        if strong:
            start_frequency = frequency_index
            strong.sort(key=lambda item: item[1], reverse=True)
            seed = sorted(
                [item[0] for item in strong[:3]],
                reverse=True,
            )
            for mode_index, velocity in enumerate(seed):
                tracks[mode_index] = velocity
                picks[frequency_index].append(velocity)
            break
    if start_frequency is None:
        return picks

    for frequency_index in range(start_frequency - 1, -1, -1):
        candidates = candidate_rows[frequency_index]
        used = set()
        selected = []
        for mode_index, previous in enumerate(tracks[: len(picks[start_frequency])]):
            if not np.isfinite(previous):
                continue
            viable = []
            for candidate_index, (velocity, amplitude) in enumerate(candidates):
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
        selected.sort(reverse=True)
        picks[frequency_index] = selected
    return picks


def candidate_states(candidate_row):
    strongest = sorted(
        candidate_row,
        key=lambda item: item[1],
        reverse=True,
    )[:6]
    states = []
    for mode_count in range(0, min(3, len(strongest)) + 1):
        for selected in combinations(strongest, mode_count):
            ordered_velocities = sorted(
                (item[0] for item in selected),
            )
            for active_modes in combinations(range(3), mode_count):
                slots = [None, None, None]
                for mode, velocity in zip(
                    active_modes,
                    ordered_velocities,
                    strict=True,
                ):
                    slots[mode] = velocity
                amplitude_sum = sum(item[1] for item in selected)
                states.append(
                    {
                        "velocities": tuple(slots),
                        "emission": amplitude_sum - 0.08 * mode_count,
                    }
                )
    # Deduplicate combinations that quantise to the same velocity tuple.
    unique = {}
    for state in states:
        key = state["velocities"]
        if key not in unique or state["emission"] > unique[key]["emission"]:
            unique[key] = state
    return list(unique.values())


def viterbi_track(candidate_rows):
    states_by_frequency = [
        candidate_states(row)
        for row in candidate_rows
    ]
    scores = {
        state["velocities"]: state["emission"]
        for state in states_by_frequency[-1]
    }
    predecessors = [None] * len(candidate_rows)
    predecessors[-1] = {}

    for frequency_index in range(len(candidate_rows) - 2, -1, -1):
        current_scores = {}
        current_predecessors = {}
        for current in states_by_frequency[frequency_index]:
            best_score = -math.inf
            best_previous = ()
            for previous in states_by_frequency[frequency_index + 1]:
                previous_key = previous["velocities"]
                current_key = current["velocities"]
                smoothness = sum(
                    ((current_key[index] - previous_key[index]) / 100.0) ** 2
                    for index in range(3)
                    if current_key[index] is not None
                    and previous_key[index] is not None
                )
                appearance_penalty = (
                    sum(
                        current_key[index] is not None
                        and previous_key[index] is None
                        for index in range(3)
                    )
                    * 0.5
                )
                candidate_score = (
                    scores[previous_key]
                    + current["emission"]
                    - 0.20 * smoothness
                    - appearance_penalty
                )
                if candidate_score > best_score:
                    best_score = candidate_score
                    best_previous = previous_key
            current_scores[current["velocities"]] = best_score
            current_predecessors[current["velocities"]] = best_previous
        scores = current_scores
        predecessors[frequency_index] = current_predecessors

    state = max(scores, key=scores.get)
    tracked = []
    for frequency_index in range(len(candidate_rows)):
        tracked.append(
            [
                float(value) if value is not None else math.nan
                for value in state
            ]
        )
        if frequency_index < len(candidate_rows) - 1:
            state = predecessors[frequency_index][state]
    return tracked


def picking_statistics(true_modes, matches):
    stats = {
        mode: {"possible": 0, "detected": 0, "errors": []}
        for mode in range(3)
    }
    false_count = 0
    for match in matches:
        if not match.get("correct", False):
            false_count += 1
    for frequency_index in range(len(e13.FREQUENCIES)):
        frequency_matches = [
            match
            for match in matches
            if match["frequency_index"] == frequency_index
        ]
        for mode_index in range(3):
            velocity = true_modes[frequency_index, mode_index]
            if not np.isfinite(velocity):
                continue
            stats[mode_index]["possible"] += 1
            candidate = min(
                frequency_matches,
                key=lambda item: abs(item["picked_velocity"] - velocity),
                default=None,
            )
            if candidate is not None and abs(
                candidate["picked_velocity"] - velocity
            ) <= e15.PICK_TOLERANCE:
                stats[mode_index]["detected"] += 1
                stats[mode_index]["errors"].append(
                    candidate["picked_velocity"] - velocity
                )
    summary = {}
    for mode_index, values in stats.items():
        summary[f"mode{mode_index}"] = {
            "possible": values["possible"],
            "detected": values["detected"],
            "detection_percent": values["detected"]
            / max(values["possible"], 1)
            * 100.0,
            "mean_absolute_error": float(np.mean(np.abs(values["errors"])))
            if values["errors"]
            else None,
        }
    summary["false_pick_percent"] = false_count / max(len(matches), 1) * 100.0
    return summary


def aggregate_method(rows, method):
    subset = [row for row in rows if row["method"] == method]
    combined_stats = {
        f"mode{mode}": {"possible": 0, "detected": 0, "errors": []}
        for mode in range(3)
    }
    total_false = 0
    total_matches = 0
    for row in subset:
        stats = row["picking"]
        for mode in range(3):
            key = f"mode{mode}"
            combined_stats[key]["possible"] += stats[key]["possible"]
            combined_stats[key]["detected"] += stats[key]["detected"]
            if stats[key]["mean_absolute_error"] is not None:
                combined_stats[key]["errors"].append(
                    stats[key]["mean_absolute_error"]
                )
        total_false += row["false_picks"]
        total_matches += row["total_picks"]
    mode_summary = {}
    for mode in range(3):
        key = f"mode{mode}"
        values = combined_stats[key]
        mode_summary[key] = {
            "possible": values["possible"],
            "detected": values["detected"],
            "detection_percent": values["detected"]
            / max(values["possible"], 1)
            * 100.0,
            "mean_absolute_pick_error_mps": float(np.mean(values["errors"]))
            if values["errors"]
            else None,
        }
    return {
        "mode_statistics": mode_summary,
        "false_pick_percent": total_false / max(total_matches, 1) * 100.0,
        "median_h_width_90": float(
            np.median([row["h_width_90"] for row in subset])
        ),
        "p90_h_width_90": float(
            np.quantile([row["h_width_90"] for row in subset], 0.90)
        ),
        "median_beta_width_90": float(
            np.median([row["beta_width_90"] for row in subset])
        ),
    }


def main() -> None:
    rng = np.random.default_rng(SEED)
    prior_precision = np.diag(1.0 / np.asarray([0.60, 100.0, 0.65]) ** 2)
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
        baseline_picks = e15.pick_ridges(image)
        candidates = detect_candidate_peaks(image)
        tracked_picks = track_modes(candidates)
        viterbi_picks = viterbi_track(candidates)
        methods = {
            "independent_peaks": baseline_picks,
            "continuity_tracking": tracked_picks,
            "viterbi_order_tracking": viterbi_picks,
        }
        for method, picks in methods.items():
            matches = e15.match_picks(center["phase"], picks)
            information, _ = e15.information_from_picks(
                center,
                plus_responses,
                minus_responses,
                matches,
            )
            information = np.nan_to_num(
                information,
                nan=0.0,
                posinf=0.0,
                neginf=0.0,
            )
            covariance = np.linalg.pinv(information + prior_precision)
            std = np.sqrt(np.maximum(np.diag(covariance), 0.0))
            rows.append(
                {
                    "model": model_index,
                    "method": method,
                    "h_true": h,
                    "total_picks": len(matches),
                    "false_picks": sum(
                        1 for match in matches if not match.get("correct", False)
                    ),
                    "h_width_90": 2.0 * Z90 * std[0],
                    "beta_width_90": 2.0 * Z90 * std[1],
                    "picking": picking_statistics(center["phase"], matches),
                }
            )
        if (model_index + 1) % 10 == 0:
            print(f"processed {model_index + 1}/{N_MODELS}")

    summary = {
        method: aggregate_method(rows, method)
        for method in (
            "independent_peaks",
            "continuity_tracking",
            "viterbi_order_tracking",
        )
    }
    with (ROOT / "e16_robust_tracking_runs.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            (
                "model",
                "method",
                "h_true",
                "total_picks",
                "false_picks",
                "h_width_90",
                "beta_width_90",
                "mode0_detection",
                "mode1_detection",
                "mode2_detection",
            )
        )
        for row in rows:
            writer.writerow(
                (
                    row["model"],
                    row["method"],
                    row["h_true"],
                    row["total_picks"],
                    row["false_picks"],
                    row["h_width_90"],
                    row["beta_width_90"],
                    row["picking"]["mode0"]["detection_percent"],
                    row["picking"]["mode1"]["detection_percent"],
                    row["picking"]["mode2"]["detection_percent"],
                )
            )
    (ROOT / "e16_robust_tracking_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
