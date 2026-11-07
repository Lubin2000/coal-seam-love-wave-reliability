"""E22: independent Love-wave forward validation on Windows.

The experiment compares three calculations:
1. the analytical dispersion equation used by the paper;
2. an independently coded transfer-matrix secular equation;
3. dispersion extracted from a 2D SH finite-difference wavefield.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys
from xml.sax.saxutils import escape

import numpy as np


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import e13_acquisition_information as e13
import multilayer_model as ml
import wavefield_model as wm


VAL_FREQUENCIES = np.arange(100.0, 561.0, 20.0)
SOURCE_FREQUENCIES = (240.0, 400.0)
CASES = (
    {"name": "thin_1p8m", "h": 1.8},
    {"name": "medium_2p8m", "h": 2.8},
    {"name": "thick_3p8m", "h": 3.8},
)
BETA_COAL = 1120.0
RHO_COAL = 1380.0
Q_COAL = 300.0
BETA_ROOF = 2400.0
RHO_ROOF = 2450.0
BETA_FLOOR = 2650.0
RHO_FLOOR = 2600.0
SEAM_TOP = 24.0
DX = 0.5
RECEIVER_SPACING = 2.0
RECEIVER_COUNT = 100
RECEIVER_START_OFFSET = 20.0
TOTAL_TIME = 0.32
X_MAX = 300.0
SOURCE_X = 30.0

MODES = (0, 1, 2)
MODE_NAMES = {0: "fundamental", 1: "first_higher", 2: "second_higher"}
MODE_COLORS = {0: "#187477", 1: "#c27a32", 2: "#75549b"}


def transfer_matrix_roots(h: float, frequency: float) -> list[float]:
    """Find trapped Love roots using the general transfer-matrix equation."""
    layers = [ml.Layer(h, BETA_COAL, RHO_COAL, Q_COAL)]
    upper = ml.Halfspace(BETA_ROOF, RHO_ROOF, 100.0)
    lower = ml.Halfspace(BETA_FLOOR, RHO_FLOOR, 120.0)
    low = BETA_COAL + 0.5
    high = min(BETA_ROOF, BETA_FLOOR) - 0.5
    linear = np.linspace(low, high, 2400)
    near_high = high - np.geomspace(1.0e-5, 8.0, 180)
    grid = np.unique(np.concatenate((linear, near_high)))
    values = np.real(
        ml.secular_values_general(layers, upper, lower, frequency, grid)
    )

    roots: list[float] = []
    for index in range(len(grid) - 1):
        left = float(values[index])
        right = float(values[index + 1])
        if not math.isfinite(left) or not math.isfinite(right):
            continue
        if left == 0.0:
            roots.append(float(grid[index]))
            continue
        if left * right >= 0.0:
            continue

        a = float(grid[index])
        b = float(grid[index + 1])
        fa = left
        for _ in range(52):
            mid = 0.5 * (a + b)
            fm = float(
                np.real(
                    ml.secular_value_general(
                        layers,
                        upper,
                        lower,
                        frequency,
                        mid,
                    )
                )
            )
            if fa * fm <= 0.0:
                b = mid
            else:
                a = mid
                fa = fm
        root = 0.5 * (a + b)
        if not roots or abs(root - roots[-1]) > 2.0:
            roots.append(root)

    # Match the fixed physical ordering used by the analytical solver:
    # increasing mode order is stored as increasing root index.
    roots.sort()
    return roots[:3]


def phase_shift_image(
    times: np.ndarray,
    offsets: np.ndarray,
    data: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute a normalised phase-shift coherence image."""
    sample_dt = float(np.median(np.diff(times)))
    trace_count, sample_count = data.shape
    nfft = 1 << int(math.ceil(math.log2(max(2048, sample_count * 8))))
    taper_time = np.hanning(sample_count)[:, None]
    spectrum = np.fft.rfft(data.T * taper_time, n=nfft, axis=0)
    frequencies = np.fft.rfftfreq(nfft, sample_dt)

    low_velocity = BETA_COAL + 3.0
    high_velocity = min(BETA_ROOF, BETA_FLOOR) - 3.0
    velocities = np.linspace(low_velocity, high_velocity, 420)
    requested_indices = [
        int(np.argmin(np.abs(frequencies - frequency)))
        for frequency in VAL_FREQUENCIES
    ]

    offset_window = np.hanning(trace_count)
    offset_window /= max(float(np.mean(offset_window)), 1e-12)
    image = np.zeros((len(VAL_FREQUENCIES), len(velocities)))

    for frequency_index, spectrum_index in enumerate(requested_indices):
        frequency = float(frequencies[spectrum_index])
        trace_spectrum = spectrum[spectrum_index]
        amplitude = np.abs(trace_spectrum)
        weights = amplitude * offset_window
        active = weights > 0.03 * max(float(np.max(weights)), 1e-12)
        if not np.any(active):
            continue
        phase = np.angle(trace_spectrum[active])
        weighted = weights[active]
        x = offsets[active]
        for velocity_index, velocity in enumerate(velocities):
            phasor = np.exp(1j * 2.0 * math.pi * frequency * x / velocity)
            image[frequency_index, velocity_index] = abs(
                np.sum(weighted * phasor)
            ) / max(float(np.sum(weighted)), 1e-12)
    return VAL_FREQUENCIES, velocities, image


def local_peaks(velocities: np.ndarray, energy: np.ndarray) -> list[tuple[float, float]]:
    # Low-frequency guided modes can be coherent but weak in absolute
    # normalised energy. Use a low floor plus a relative-to-maximum condition
    # and retain the peak value for a separate contrast test.
    threshold = max(0.006, 0.03 * float(np.max(energy)))
    candidates = []
    for index in range(1, len(energy) - 1):
        if (
            energy[index] >= threshold
            and energy[index] >= energy[index - 1]
            and energy[index] > energy[index + 1]
        ):
            candidates.append(
                (float(velocities[index]), float(energy[index]))
            )
    candidates.sort(key=lambda item: item[1], reverse=True)

    peaks: list[tuple[float, float]] = []
    for velocity, value in candidates:
        if all(abs(velocity - kept[0]) > 25.0 for kept in peaks):
            peaks.append((velocity, value))
    peaks.sort(reverse=True)
    return peaks


def root_sets_match(
    analytical: list[float],
    transfer: list[float],
) -> tuple[bool, float]:
    if len(analytical) != len(transfer):
        return False, math.nan
    if not analytical:
        return True, 0.0
    differences = [
        abs(a - b) for a, b in zip(analytical, transfer, strict=True)
    ]
    return max(differences) < 2.0, max(differences)


def validate_case(
    case: dict,
) -> tuple[list[dict], dict, np.ndarray, np.ndarray]:
    h = float(case["h"])
    layers = [ml.Layer(h, BETA_COAL, RHO_COAL, Q_COAL)]
    upper = ml.Halfspace(BETA_ROOF, RHO_ROOF, 100.0)
    lower = ml.Halfspace(BETA_FLOOR, RHO_FLOOR, 120.0)

    images = []
    for source_frequency in SOURCE_FREQUENCIES:
        print(
            f"{case['name']}: SH FD source {source_frequency:.0f} Hz",
            flush=True,
        )
        times, offsets, data = wm.simulate_shot(
            layers,
            upper,
            lower,
            seam_top=SEAM_TOP,
            dx=DX,
            receiver_spacing=RECEIVER_SPACING,
            receiver_count=RECEIVER_COUNT,
            total_time=TOTAL_TIME,
            source_frequency=source_frequency,
            noise_level=0.0,
            x_max=X_MAX,
            source_x=SOURCE_X,
            receiver_start_offset=RECEIVER_START_OFFSET,
        )
        frequencies, velocities, image = phase_shift_image(
            times,
            offsets,
            data,
        )
        images.append(image)

    combined = np.maximum.reduce(images)
    rows: list[dict] = []
    root_matches = 0
    root_counts = 0

    for frequency_index, frequency in enumerate(VAL_FREQUENCIES):
        analytical = e13.all_phase_velocities(
            float(frequency),
            h,
            BETA_COAL,
        )
        transfer = transfer_matrix_roots(h, float(frequency))
        matched, root_error = root_sets_match(analytical, transfer)
        root_counts += 1
        root_matches += int(matched)

        peaks = local_peaks(velocities, combined[frequency_index])
        for mode_index, theoretical_velocity in enumerate(analytical[:3]):
            resolution = 1.5 * (
                theoretical_velocity**2
                / max(
                    float(frequency)
                    * (RECEIVER_COUNT - 1)
                    * RECEIVER_SPACING,
                    1e-9,
                )
            )
            tolerance = max(50.0, resolution)
            if peaks:
                nearest_velocity, nearest_energy = min(
                    peaks,
                    key=lambda item: abs(item[0] - theoretical_velocity),
                )
                pick_error = abs(nearest_velocity - theoretical_velocity)
                background = float(np.median(combined[frequency_index]))
                contrast = nearest_energy / max(background, 0.002)
                detected = (
                    pick_error <= tolerance
                    and nearest_energy >= 0.006
                    and contrast >= 1.5
                )
            else:
                nearest_velocity = math.nan
                nearest_energy = 0.0
                pick_error = math.nan
                contrast = 0.0
                detected = False

            transfer_velocity = (
                transfer[mode_index]
                if mode_index < len(transfer)
                else math.nan
            )
            rows.append(
                {
                    "case": case["name"],
                    "h_true_m": h,
                    "frequency_hz": float(frequency),
                    "mode": mode_index,
                    "mode_name": MODE_NAMES[mode_index],
                    "analytical_velocity": float(theoretical_velocity),
                    "transfer_velocity": transfer_velocity,
                    "transfer_error": (
                        abs(transfer_velocity - theoretical_velocity)
                        if math.isfinite(transfer_velocity)
                        else math.nan
                    ),
                    "root_count_match": matched,
                    "root_set_max_error": root_error,
                    "observed_peak_velocity": nearest_velocity,
                    "observed_peak_error": pick_error,
                    "pick_tolerance": tolerance,
                    "pick_energy": nearest_energy,
                    "pick_contrast": contrast,
                    "detected": detected,
                }
            )

    detected_by_mode = {
        MODE_NAMES[mode]: int(
            sum(
                row["detected"]
                for row in rows
                if row["mode"] == mode
            )
        )
        for mode in MODES
    }
    predicted_by_mode = {
        MODE_NAMES[mode]: int(sum(row["mode"] == mode for row in rows))
        for mode in MODES
    }
    summary = {
        "case": case["name"],
        "h_true_m": h,
        "root_validation": {
            "frequency_points": root_counts,
            "mode_count_matches": root_matches,
            "mode_count_match_percent": 100.0 * root_matches / max(root_counts, 1),
        },
        "wavefield_validation": {
            "predicted_by_mode": predicted_by_mode,
            "detected_by_mode": detected_by_mode,
            "detection_percent_by_mode": {
                name: 100.0
                * detected_by_mode[name]
                / max(predicted_by_mode[name], 1)
                for name in detected_by_mode
            },
            "median_pick_error_by_mode": {
                name: float(
                    np.nanmedian(
                        [
                            row["observed_peak_error"]
                            for row in rows
                            if row["mode_name"] == name
                        ]
                    )
                )
                for name in detected_by_mode
            },
        },
    }
    return rows, summary, combined, velocities


def make_svg(
    rows: list[dict],
    summaries: list[dict],
    images: dict[str, np.ndarray],
    velocities: np.ndarray,
) -> str:
    width, height = 1100, 440
    panel_width = 330
    panel_height = 310
    left_margin = 25
    top = 70
    body = [
        '<rect width="100%" height="100%" fill="#f8f9f7"/>',
        '<style>text{font-family:Arial,"Liberation Sans",sans-serif;}'
        ".axis{stroke:#b8c2c4;fill:none}.grid{stroke:#e3e8e7;stroke-width:1}</style>",
        '<text x="25" y="32" font-size="20" font-weight="600" fill="#18222c">'
        "Love-wave finite-difference validation</text>",
        '<text x="25" y="53" font-size="12" fill="#66727d">'
        "Background: phase-shift coherence; lines: analytical roots; circles: matched independent-wavefield peaks</text>",
    ]
    fmin, fmax = 90.0, 570.0
    cmin, cmax = 1050.0, 2470.0

    def x_map(value):
        return left_margin + (value - fmin) / (fmax - fmin) * panel_width

    def y_map(value):
        return top + panel_height - (value - cmin) / (cmax - cmin) * panel_height

    for panel_index, case in enumerate(CASES):
        x_offset = 25 + panel_index * 355
        name = case["name"]
        summary = next(item for item in summaries if item["case"] == name)
        body.append(
            f'<text x="{x_offset}" y="82" font-size="14" font-weight="600" fill="#18222c">'
            f"{escape(name)} (h={case['h']:.1f} m)</text>"
        )
        body.append(
            f'<rect x="{x_offset}" y="{top}" width="{panel_width}" height="{panel_height}" '
            'class="axis"/>'
        )

        image = images[name]
        decimation = 6
        usable_columns = (image.shape[1] // decimation) * decimation
        display_image = image[:, :usable_columns].reshape(
            image.shape[0],
            usable_columns // decimation,
            decimation,
        ).mean(axis=2)
        display_velocities = velocities[:usable_columns].reshape(
            usable_columns // decimation,
            decimation,
        ).mean(axis=1)
        color_limit = max(float(np.quantile(display_image, 0.995)), 0.01)
        row_count, column_count = display_image.shape
        cell_width = panel_width / row_count
        cell_height = panel_height / column_count
        for frequency_index in range(row_count):
            for velocity_index in range(column_count):
                intensity = min(
                    1.0,
                    math.sqrt(
                        max(float(display_image[frequency_index, velocity_index]), 0.0)
                        / color_limit
                    ),
                )
                if intensity < 0.025:
                    continue
                red = int(248 - 224 * intensity)
                green = int(249 - 133 * intensity)
                blue = int(247 - 128 * intensity)
                velocity = float(display_velocities[velocity_index]) - 0.5 * (
                    display_velocities[1] - display_velocities[0]
                )
                cell_y = y_map(
                    velocity
                    + (display_velocities[1] - display_velocities[0])
                )
                body.append(
                    f'<rect x="{x_offset + frequency_index * cell_width:.1f}" '
                    f'y="{cell_y:.1f}" width="{cell_width + 0.2:.1f}" '
                    f'height="{cell_height + 0.2:.1f}" '
                    f'fill="rgb({red},{green},{blue})" stroke="none"/>'
                )
        for velocity in (1200, 1600, 2000, 2400):
            y = y_map(velocity)
            body.append(
                f'<path d="M{x_offset} {y:.1f} H{x_offset + panel_width}" class="grid"/>'
            )
            body.append(
                f'<text x="{x_offset - 6}" y="{y + 4:.1f}" text-anchor="end" '
                f'font-size="10" fill="#66727d">{velocity}</text>'
            )
        for frequency in (100, 200, 300, 400, 500):
            x = x_offset + (frequency - fmin) / (fmax - fmin) * panel_width
            body.append(
                f'<path d="M{x:.1f} {top} V{top + panel_height}" class="grid"/>'
            )
            body.append(
                f'<text x="{x:.1f}" y="{top + panel_height + 16}" '
                f'text-anchor="middle" font-size="10" fill="#66727d">{frequency}</text>'
            )

        case_rows = [row for row in rows if row["case"] == name]
        for mode in MODES:
            mode_rows = [row for row in case_rows if row["mode"] == mode]
            points = [
                (
                    x_offset
                    + (row["frequency_hz"] - fmin) / (fmax - fmin) * panel_width,
                    y_map(row["analytical_velocity"]),
                )
                for row in mode_rows
            ]
            if len(points) > 1:
                path = " ".join(
                    f"{'M' if index == 0 else 'L'}{x:.1f},{y:.1f}"
                    for index, (x, y) in enumerate(points)
                )
                body.append(
                    f'<path d="{path}" fill="none" stroke="{MODE_COLORS[mode]}" '
                    'stroke-width="1.8"/>'
                )
            for row in mode_rows:
                if not row["detected"]:
                    continue
                x = (
                    x_offset
                    + (row["frequency_hz"] - fmin) / (fmax - fmin) * panel_width
                )
                y = y_map(row["observed_peak_velocity"])
                body.append(
                    f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.8" '
                    f'fill="{MODE_COLORS[mode]}"/>'
                )

        rates = summary["wavefield_validation"]["detection_percent_by_mode"]
        body.append(
            f'<text x="{x_offset + 8}" y="{top + panel_height - 10}" font-size="10" fill="#66727d">'
            f"detected: M0 {rates['fundamental']:.0f}% | "
            f"M1 {rates['first_higher']:.0f}% | M2 {rates['second_higher']:.0f}%</text>"
        )

    body.extend(
        [
            f'<text x="{left_margin + panel_width / 2:.1f}" y="415" text-anchor="middle" font-size="11" fill="#66727d">Frequency (Hz)</text>',
            f'<text x="14" y="{top + panel_height / 2:.1f}" transform="rotate(-90 14 {top + panel_height / 2:.1f})" text-anchor="middle" font-size="11" fill="#66727d">Phase velocity (m/s)</text>',
            f'<text x="{25 + 2 * 355 + 8}" y="{top - 8}" font-size="10" fill="#187477">mode 0</text>',
            f'<text x="{25 + 2 * 355 + 65}" y="{top - 8}" font-size="10" fill="#c27a32">mode 1</text>',
            f'<text x="{25 + 2 * 355 + 122}" y="{top - 8}" font-size="10" fill="#75549b">mode 2</text>',
        ]
    )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        + "".join(body)
        + "</svg>"
    )


def main() -> None:
    all_rows: list[dict] = []
    summaries = []
    images = {}
    phase_velocities = None
    for case in CASES:
        rows, summary, image, velocities = validate_case(case)
        all_rows.extend(rows)
        summaries.append(summary)
        images[case["name"]] = image
        phase_velocities = velocities
        print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)

    fields = list(all_rows[0])
    with (ROOT / "e22_love_forward_validation.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_rows)

    def finite_values(values):
        return [value for value in values if math.isfinite(value)]

    aggregate = {}
    for mode in MODES:
        mode_rows = [row for row in all_rows if row["mode"] == mode]
        detected_rows = [row for row in mode_rows if row["detected"]]
        errors = finite_values(
            [row["observed_peak_error"] for row in detected_rows]
        )
        aggregate[MODE_NAMES[mode]] = {
            "predicted_frequency_points": len(mode_rows),
            "detected_frequency_points": len(detected_rows),
            "detection_percent": 100.0
            * len(detected_rows)
            / max(len(mode_rows), 1),
            "median_detected_velocity_error_mps": (
                float(np.median(errors)) if errors else math.nan
            ),
            "max_detected_velocity_error_mps": (
                float(np.max(errors)) if errors else math.nan
            ),
        }

    banded = {}
    for band_name, lower, upper in (
        ("100-240_hz", 100.0, 240.0),
        ("260-400_hz", 260.0, 400.0),
        ("420-560_hz", 420.0, 560.0),
    ):
        banded[band_name] = {}
        for mode in MODES:
            mode_rows = [
                row
                for row in all_rows
                if row["mode"] == mode
                and lower <= row["frequency_hz"] <= upper
            ]
            banded[band_name][MODE_NAMES[mode]] = {
                "predicted": len(mode_rows),
                "detected": int(sum(row["detected"] for row in mode_rows)),
            }

    summary = {
        "method": {
            "wavefield": "2D SH heterogeneous finite difference",
            "dispersion_transform": "normalized phase-shift coherence",
            "analytical": "three-layer Love secular equation",
            "transfer_matrix": "independent Haskell-Thomson state transfer",
            "source_frequencies_hz": list(SOURCE_FREQUENCIES),
            "receiver_count": RECEIVER_COUNT,
            "receiver_spacing_m": RECEIVER_SPACING,
            "receiver_start_offset_m": RECEIVER_START_OFFSET,
            "grid_spacing_m": DX,
            "record_length_s": TOTAL_TIME,
            "model_length_m": X_MAX,
            "validation_frequencies_hz": VAL_FREQUENCIES.tolist(),
        },
        "aggregate_wavefield_validation": aggregate,
        "banded_wavefield_validation": banded,
        "cases": summaries,
    }
    (ROOT / "e22_love_forward_validation_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    np.savez_compressed(
        ROOT / "e22_phase_shift_images.npz",
        frequencies=VAL_FREQUENCIES,
        velocities=phase_velocities,
        **images,
    )
    svg = make_svg(
        all_rows,
        summaries,
        images,
        phase_velocities,
    )
    (ROOT / "e22_love_forward_validation.svg").write_text(
        svg,
        encoding="utf-8",
    )
    paper_figure_directory = ROOT.parent / "paper" / "figures"
    paper_figure_directory.mkdir(parents=True, exist_ok=True)
    (paper_figure_directory / "fig6_wavefield_validation.svg").write_text(
        svg,
        encoding="utf-8",
    )
    print("E22 validation complete.", flush=True)


if __name__ == "__main__":
    main()
