"""E25: compare analytical coal-seam Love roots with QEDispInv."""

from __future__ import annotations

from collections import defaultdict
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


VALIDATION_DIR = ROOT.parent / "paper" / "qedispinv_reference"
QEDISP_PATH = VALIDATION_DIR / "qedisp_love_all.txt"
BETA_COAL = 1120.0
H = 2.4
MODES = (0, 1, 2)
COLORS = {0: "#187477", 1: "#c27a32", 2: "#75549b"}


def load_qedisp():
    values: dict[float, list[tuple[float, int]]] = defaultdict(list)
    for line in QEDISP_PATH.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        frequency, velocity_kmps, mode = line.split()
        values[round(float(frequency), 6)].append(
            (float(velocity_kmps) * 1000.0, int(mode))
        )
    return values


def make_svg(rows: list[dict]) -> str:
    width, height = 980, 460
    left, top = 95, 65
    plot_width, plot_height = 810, 310
    fmin, fmax = 90.0, 570.0
    cmin, cmax = 1100.0, 2450.0

    def x_map(value):
        return left + (value - fmin) / (fmax - fmin) * plot_width

    def y_map(value):
        return top + plot_height - (value - cmin) / (cmax - cmin) * plot_height

    body = [
        '<rect width="100%" height="100%" fill="#f8f9f7"/>',
        '<style>text{font-family:Arial,"Liberation Sans",sans-serif;}'
        '.grid{stroke:#e1e6e5;stroke-width:1}.axis{stroke:#aeb9bb;fill:none}</style>',
        '<text x="28" y="32" font-size="20" font-weight="600" fill="#18222c">'
        "QEDispInv cross-check for coal-seam Love modes</text>",
        '<text x="28" y="52" font-size="12" fill="#66727d">'
        "Lines: analytical roots; circles: independently computed QEDispInv roots</text>",
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" class="axis"/>',
    ]
    for velocity in (1200, 1600, 2000, 2400):
        y = y_map(velocity)
        body.append(f'<path d="M{left} {y:.1f} H{left + plot_width}" class="grid"/>')
        body.append(
            f'<text x="{left - 10}" y="{y + 4:.1f}" text-anchor="end" '
            f'font-size="10" fill="#66727d">{velocity}</text>'
        )
    for frequency in (100, 200, 300, 400, 500):
        x = x_map(frequency)
        body.append(f'<path d="M{x:.1f} {top} V{top + plot_height}" class="grid"/>')
        body.append(
            f'<text x="{x:.1f}" y="{top + plot_height + 17}" '
            f'text-anchor="middle" font-size="10" fill="#66727d">{frequency}</text>'
        )

    for mode in MODES:
        mode_rows = [row for row in rows if row["mode"] == mode]
        if not mode_rows:
            continue
        points = [
            (
                x_map(row["frequency_hz"]),
                y_map(row["analytical_velocity_mps"]),
            )
            for row in mode_rows
        ]
        path = " ".join(
            f"{'M' if index == 0 else 'L'}{x:.1f},{y:.1f}"
            for index, (x, y) in enumerate(points)
        )
        body.append(
            f'<path d="{path}" fill="none" stroke="{COLORS[mode]}" stroke-width="2"/>'
        )
        for row in mode_rows:
            x = x_map(row["frequency_hz"])
            y = y_map(row["qedisp_velocity_mps"])
            body.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3" '
                f'fill="none" stroke="{COLORS[mode]}" stroke-width="1.5"/>'
            )

    for index, mode in enumerate(MODES):
        x = left + 580 + index * 70
        body.append(f'<circle cx="{x}" cy="398" r="4" fill="{COLORS[mode]}"/>')
        body.append(
            f'<text x="{x + 9}" y="402" font-size="11" fill="#18222c">mode {mode}</text>'
        )
    body.extend(
        [
            f'<text x="{left + plot_width / 2:.1f}" y="427" text-anchor="middle" '
            'font-size="11" fill="#66727d">Frequency (Hz)</text>',
            f'<text x="16" y="{top + plot_height / 2:.1f}" '
            f'transform="rotate(-90 16 {top + plot_height / 2:.1f})" '
            'text-anchor="middle" font-size="11" fill="#66727d">Phase velocity (m/s)</text>',
        ]
    )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        + "".join(body)
        + "</svg>"
    )


def main() -> None:
    qedisp = load_qedisp()
    rows = []
    for frequency in range(100, 561, 20):
        analytical = e13.all_phase_velocities(
            float(frequency),
            H,
            BETA_COAL,
        )
        qedisp_roots = qedisp[float(frequency)]
        for mode, analytical_velocity in enumerate(analytical[:3]):
            qedisp_velocity, qedisp_mode = min(
                qedisp_roots,
                key=lambda item: abs(item[0] - analytical_velocity),
            )
            rows.append(
                {
                    "frequency_hz": float(frequency),
                    "mode": mode,
                    "analytical_velocity_mps": float(analytical_velocity),
                    "qedisp_velocity_mps": float(qedisp_velocity),
                    "qedisp_native_mode": qedisp_mode,
                    "absolute_error_mps": abs(
                        float(qedisp_velocity - analytical_velocity)
                    ),
                }
            )

    errors = np.asarray([row["absolute_error_mps"] for row in rows])
    summary = {
        "software": "QEDispInv",
        "upstream_commit": "ed1b5dd7",
        "compiler": "GCC 16.2.0 UCRT64",
        "build_options": {"USE_OPENMP": False},
        "model": {
            "h_m": H,
            "beta_coal_mps": BETA_COAL,
            "roof_thickness_m": 100.0,
            "note": "The finite 100-m roof approximates an upper half-space. QEDispInv also returns free-surface roof modes; each analytical coal mode is matched to its nearest independently computed root.",
        },
        "frequency_range_hz": [100.0, 560.0],
        "modal_observations": len(rows),
        "median_absolute_error_mps": float(np.median(errors)),
        "mean_absolute_error_mps": float(np.mean(errors)),
        "p95_absolute_error_mps": float(np.quantile(errors, 0.95)),
        "max_absolute_error_mps": float(np.max(errors)),
        "within_1_mps": int(np.sum(errors <= 1.0)),
        "within_10_mps": int(np.sum(errors <= 10.0)),
    }
    with (ROOT / "e25_qedispinv_validation.csv").open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (ROOT / "e25_qedispinv_validation_summary.json").write_text(
        json.dumps({"summary": summary, "rows": rows}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    svg = make_svg(rows)
    (ROOT / "e25_qedispinv_validation.svg").write_text(
        svg,
        encoding="utf-8",
    )
    (ROOT.parent / "paper" / "figures" / "fig7_qedispinv_validation.svg").write_text(
        svg,
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
