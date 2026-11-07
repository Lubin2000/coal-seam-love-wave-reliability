"""E26: add a connected-component dispersion-image segmentation baseline."""

from __future__ import annotations

from collections import deque
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


SEED = e16.SEED
N_MODELS = e16.N_MODELS
THRESHOLD = 0.12
MIN_COMPONENT_PIXELS = 8
MIN_COMPONENT_FREQUENCIES = 4


def connected_components(mask: np.ndarray):
    visited = np.zeros_like(mask, dtype=bool)
    components = []
    row_count, column_count = mask.shape
    for row in range(row_count):
        for column in range(column_count):
            if not mask[row, column] or visited[row, column]:
                continue
            queue = deque([(row, column)])
            visited[row, column] = True
            pixels = []
            while queue:
                current_row, current_column = queue.popleft()
                pixels.append((current_row, current_column))
                for row_step in (-1, 0, 1):
                    for column_step in (-1, 0, 1):
                        if row_step == 0 and column_step == 0:
                            continue
                        next_row = current_row + row_step
                        next_column = current_column + column_step
                        if (
                            0 <= next_row < row_count
                            and 0 <= next_column < column_count
                            and mask[next_row, next_column]
                            and not visited[next_row, next_column]
                        ):
                            visited[next_row, next_column] = True
                            queue.append((next_row, next_column))
            components.append(pixels)
    return components


def segmentation_picks(image: np.ndarray) -> list[list[float]]:
    """Pick ridge peaks from sufficiently persistent bright components."""
    mask = image >= THRESHOLD
    picks = [[] for _ in range(len(image))]
    for component in connected_components(mask):
        rows = {row for row, _ in component}
        if (
            len(component) < MIN_COMPONENT_PIXELS
            or len(rows) < MIN_COMPONENT_FREQUENCIES
        ):
            continue
        pixels_by_row = {}
        for row, column in component:
            pixels_by_row.setdefault(row, []).append(column)
        for row, columns in pixels_by_row.items():
            columns = sorted(columns)
            candidates = []
            for column in columns:
                left = image[row, column - 1] if column > 0 else -1.0
                right = (
                    image[row, column + 1]
                    if column + 1 < image.shape[1]
                    else -1.0
                )
                if image[row, column] >= left and image[row, column] >= right:
                    candidates.append(
                        (
                            float(e15.VELOCITY_GRID[column]),
                            float(image[row, column]),
                        )
                    )
            candidates.sort(key=lambda item: item[1], reverse=True)
            for velocity, _ in candidates:
                if all(
                    abs(velocity - other) > 45.0
                    for other in picks[row]
                ):
                    picks[row].append(velocity)
                if len(picks[row]) == 3:
                    break
    return [sorted(values, reverse=True) for values in picks]


def information_width(center, plus, minus, picks):
    matches = e15.match_picks(center["phase"], picks)
    information, count = e15.information_from_picks(
        center,
        plus,
        minus,
        matches,
    )
    information = np.nan_to_num(
        information,
        nan=0.0,
        posinf=0.0,
        neginf=0.0,
    )
    prior_precision = np.diag(1.0 / np.asarray([0.60, 100.0, 0.65]) ** 2)
    covariance = np.linalg.pinv(information + prior_precision)
    std = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    return matches, float(2.0 * e13.Z90 * std[0])


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

        plus = []
        minus = []
        for parameter_index, step in enumerate((0.004, 1.0, 0.008)):
            trial_plus = parameters.copy()
            trial_minus = parameters.copy()
            trial_plus[parameter_index] += step
            trial_minus[parameter_index] -= step
            plus.append(e13.response_for_parameters(trial_plus))
            minus.append(e13.response_for_parameters(trial_minus))

        image = e15.make_dispersion_image(center["phase"], rng)
        candidates = e16.detect_candidate_peaks(image)
        methods = {
            "independent_peaks": e15.pick_ridges(image),
            "greedy_tracking": e16.track_modes(candidates),
            "viterbi_order_tracking": e16.viterbi_track(candidates),
            "component_segmentation": segmentation_picks(image),
        }
        for method, picks in methods.items():
            matches, h_width = information_width(
                center,
                plus,
                minus,
                picks,
            )
            statistics = e16.picking_statistics(center["phase"], matches)
            rows.append(
                {
                    "model": model_index,
                    "method": method,
                    "h_true": h,
                    "total_picks": len(matches),
                    "false_picks": sum(
                        not match.get("correct", False)
                        for match in matches
                    ),
                    "h_width_90": h_width,
                    **statistics,
                }
            )
        if (model_index + 1) % 10 == 0:
            print(f"processed {model_index + 1}/{N_MODELS}", flush=True)

    flattened = []
    for row in rows:
        flat = {
            key: value
            for key, value in row.items()
            if not isinstance(value, dict)
        }
        for mode in range(3):
            flat[f"mode{mode}_detection"] = row[f"mode{mode}"][
                "detection_percent"
            ]
            flat[f"mode{mode}_possible"] = row[f"mode{mode}"]["possible"]
            flat[f"mode{mode}_detected"] = row[f"mode{mode}"]["detected"]
        flattened.append(flat)
    with (ROOT / "e26_additional_picker_runs.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flattened[0]))
        writer.writeheader()
        writer.writerows(flattened)

    summary = {}
    for method in methods:
        subset = [row for row in rows if row["method"] == method]
        summary[method] = {
            "mode_detection_percent": {
                f"mode{mode}": float(
                    100.0
                    * sum(row[f"mode{mode}"]["detected"] for row in subset)
                    / max(
                        sum(
                            row[f"mode{mode}"]["possible"]
                            for row in subset
                        ),
                        1,
                    )
                )
                for mode in range(3)
            },
            "false_pick_percent": float(
                100.0
                * sum(row["false_picks"] for row in subset)
                / max(
                    sum(row["total_picks"] for row in subset),
                    1,
                )
            ),
            "median_h_width_90": float(
                np.median([row["h_width_90"] for row in subset])
            ),
        }
    (ROOT / "e26_additional_picker_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
