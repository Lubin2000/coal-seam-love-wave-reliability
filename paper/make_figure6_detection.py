"""Create a non-misleading finite-difference detection figure for paper 1."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys
from xml.sax.saxutils import escape


REPOSITORY = Path(__file__).resolve().parents[1]
ROOT = REPOSITORY / "analysis"
SOURCE_CSV = ROOT / "e22_love_forward_validation.csv"
SOURCE_JSON = ROOT / "e22_love_forward_validation_summary.json"
TARGET = REPOSITORY / "paper" / "figures" / "fig6_wavefield_validation.svg"

CASES = (
    ("thin_1p8m", "Thin seam (h = 1.8 m)"),
    ("medium_2p8m", "Medium seam (h = 2.8 m)"),
    ("thick_3p8m", "Thick seam (h = 3.8 m)"),
)
COLORS = {0: "#187477", 1: "#c27a32", 2: "#75549b"}
INKS = {
    "ink": "#18222c",
    "muted": "#66727d",
    "grid": "#dfe4e3",
    "paper": "#f8f9f7",
}


def text(x, y, value, size=12, weight=400, color=None, anchor="start"):
    fill = color or INKS["ink"]
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">'
        f"{escape(str(value))}</text>"
    )


def load_rows():
    rows = []
    numeric = {
        "h_true_m",
        "frequency_hz",
        "mode",
        "analytical_velocity",
        "observed_peak_velocity",
        "observed_peak_error",
        "pick_tolerance",
    }
    with SOURCE_CSV.open(encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row = {}
            for key, value in raw.items():
                if key in numeric and value not in ("", "nan", "NaN"):
                    row[key] = float(value)
                elif key == "detected":
                    row[key] = value == "True"
                else:
                    row[key] = value
            rows.append(row)
    return rows


def main():
    rows = load_rows()
    summary = json.loads(SOURCE_JSON.read_text(encoding="utf-8"))
    width, height = 1100, 475
    panel_width, panel_height = 330, 300
    top, left_margin = 80, 35
    fmin, fmax = 90.0, 570.0
    cmin, cmax = 1050.0, 2470.0

    def x_map(value):
        return left_margin + (value - fmin) / (fmax - fmin) * panel_width

    def y_map(value):
        return top + panel_height - (value - cmin) / (cmax - cmin) * panel_height

    body = [
        '<rect width="100%" height="100%" fill="#f8f9f7"/>',
        '<style>text{font-family:Arial,"Liberation Sans",sans-serif;}</style>',
        text(
            28,
            27,
            "Lines: analytical branches; circles: accepted coherent peaks",
            16,
            400,
            "#48545d",
        ),
    ]

    for panel_index, (case, case_label) in enumerate(CASES):
        x_offset = 35 + panel_index * 355
        case_rows = [row for row in rows if row["case"] == case]
        case_summary = next(
            item for item in summary["cases"] if item["case"] == case
        )
        body.append(
            text(
                x_offset,
                65,
                case_label,
                18,
                600,
            )
        )
        body.append(
            f'<rect x="{x_offset}" y="{top}" width="{panel_width}" '
            f'height="{panel_height}" fill="#ffffff" stroke="#b8c2c4"/>'
        )
        for velocity in (1200, 1600, 2000, 2400):
            y = y_map(velocity)
            body.append(
                f'<path d="M{x_offset} {y:.1f} H{x_offset + panel_width}" '
                f'stroke="{INKS["grid"]}"/>'
            )
            body.append(
                text(
                    x_offset - 6,
                    y + 4,
                    velocity,
                    16,
                    400,
                    "#48545d",
                    "end",
                )
            )
        for frequency in (100, 200, 300, 400, 500):
            x = x_offset + (frequency - fmin) / (fmax - fmin) * panel_width
            body.append(
                f'<path d="M{x:.1f} {top} V{top + panel_height}" '
                f'stroke="{INKS["grid"]}"/>'
            )
            body.append(
                text(
                    x,
                    top + panel_height + 17,
                    frequency,
                    16,
                    400,
                    "#48545d",
                    "middle",
                )
            )

        for mode in range(3):
            mode_rows = [row for row in case_rows if int(row["mode"]) == mode]
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
                    f'<path d="{path}" fill="none" stroke="{COLORS[mode]}" '
                    'stroke-width="2"/>'
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
                    f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.2" '
                    f'fill="{COLORS[mode]}"/>'
                )
        rates = case_summary["wavefield_validation"]["detection_percent_by_mode"]
        body.append(
            text(
                x_offset + 8,
                top + panel_height - 10,
                f"detected: M0 {rates['fundamental']:.0f}% | "
                f"M1 {rates['first_higher']:.0f}% | "
                f"M2 {rates['second_higher']:.0f}%",
                16,
                500,
                "#48545d",
            )
        )

    for mode, label in enumerate(
        ("mode 0", "mode 1", "mode 2")
    ):
        x = 670 + mode * 125
        body.append(f'<circle cx="{x}" cy="54" r="5" fill="{COLORS[mode]}"/>')
        body.append(text(x + 10, 32, label, 16))

    body.extend(
        [
            text(
                left_margin + panel_width / 2,
                414,
                "Frequency (Hz)",
                16,
                500,
                anchor="middle",
            ),
            text(
                15,
                top + panel_height / 2,
                "Phase velocity (m/s)",
                16,
                500,
                "#48545d",
            ),
            text(
                28,
                451,
                "The raw coherence image is shown separately as Supplementary Figure S1 because it also contains strong nonmodal energy.",
                16,
                400,
                "#48545d",
            ),
        ]
    )

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}" role="img">'
        + "".join(body)
        + "</svg>"
    )
    TARGET.write_text(svg, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
