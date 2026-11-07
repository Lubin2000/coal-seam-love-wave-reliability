"""Generate publication-oriented SVG figures for paper 1."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import re
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parents[1] / "analysis"
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(parents=True, exist_ok=True)


COLORS = {
    "ink": "#18222c",
    "muted": "#66727d",
    "grid": "#dfe4e3",
    "teal": "#187477",
    "gold": "#c27a32",
    "gray": "#8a949c",
    "paper": "#f8f9f7",
}


def text(x, y, value, size=12, weight=400, color=None, anchor="start"):
    fill = color or COLORS["ink"]
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">'
        f"{escape(str(value))}</text>"
    )


def save(name, width, height, body):
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}" '
        'role="img">'
        f'<rect width="100%" height="100%" fill="{COLORS["paper"]}"/>'
        '<style>text{font-family:"Microsoft YaHei","Noto Sans SC",Arial,sans-serif;}</style>'
        f"{body}</svg>"
    )
    replacements = {
        "图1  多模态槽波信息保留与稳健拾取流程": "Figure 1. Multimodal information retention and robust picking workflow",
        "三层煤层模型": "Three-layer model",
        "厚度 / 速度 / Q": "Thickness / Vs / Q",
        "多模态正演": "Multimodal forward",
        "频散 / 衰减 / 群速度": "Dispersion / atten. / group",
        "合成能量图": "Synthetic image",
        "噪声 / 漏检 / 伪脊": "Noise / gaps / false ridges",
        "全局模态跟踪": "Global tracking",
        "顺序 / 连续性 / 消失": "Order / continuity / gaps",
        "稳健局部反演": "Robust inversion",
        "信息与可信度": "Information quality",
        "宽度 / 覆盖率 / 误差": "Width / coverage / error",
        "合成基准用于识别信息瓶颈；真实工区验证仍是该类方法外推的前提。": "The synthetic benchmark identifies information bottlenecks; field validation remains necessary.",
        "图2  不同观测组合的厚度信息预算": "Figure 2. Thickness information budget for different observation designs",
        "理论 90% 厚度区间中位数 / m": "Median theoretical 90% thickness interval / m",
        "基阶 + 1阶 + 2阶": "Fund. + 1st + 2nd higher",
        "基阶 + 1阶 + 衰减": "Fund. + 1st higher + Q",
        "多阶 + 群速度 + 衰减": "Higher modes + group + Q",
        "基阶 + 1阶": "Fund. + 1st higher",
        "基阶 + 衰减": "Fundamental + Q",
        "基阶": "Fundamental",
        "Fisher 信息下界；相同正演用于数据生成与反演，属于乐观估计。": "Fisher-information lower bound; the same forward model is used for generation and inversion, so this is optimistic.",
        "图3  基阶与高阶模态自动检出率": "Figure 3. Automatic detection of fundamental and higher-mode branches",
        "独立峰": "Independent peaks",
        "贪心跟踪": "Greedy tracking",
        "全局跟踪": "Global tracking",
        "分割基线": "Segmentation",
        "第一高阶": "First higher",
        "第二高阶": "Second higher",
        "图4  自动模态拾取下的厚度稳健性": "Figure 4. Thickness robustness under automatic modal picking",
        "额外离群点比例": "Added outlier fraction",
        "厚度 RMSE / m": "Thickness RMSE / m",
    }
    english = svg
    for source, target in replacements.items():
        english = english.replace(source, target)
    english = re.sub(r"伪拾取 ([0-9.]+)%", r"False picks \1%", english)
    figure_title = {
        "fig1_workflow": "Figure 1. Multimodal information retention and robust picking workflow",
        "fig2_information_budget": "Figure 2. Thickness information budget for different observation designs",
        "fig3_mode_picking": "Figure 3. Automatic detection of fundamental and higher-mode branches",
        "fig4_robust_likelihood": "Figure 4. Thickness robustness under automatic modal picking",
        "fig5_calibration": "Figure 5. Simulated interval coverage",
        "fig8_sls_attenuation": "Figure 8. Viscoelastic Love-wave attenuation check",
    }[Path(name).stem]
    english = re.sub(
        r'(<text x="34\.0" y="34\.0"[^>]*>).*?(</text>)',
        lambda match: match.group(1) + figure_title + match.group(2),
        english,
        count=1,
    )
    english = english.replace(
        '"Microsoft YaHei","Noto Sans SC",Arial,sans-serif',
        'Arial,"Liberation Sans",sans-serif',
    )
    stem = Path(name).stem
    english_name = f"{stem}_en.svg"
    (OUT / english_name).write_text(english, encoding="utf-8")


def figure_workflow():
    width, height = 1100, 300
    labels = [
        ("三层煤层模型", "厚度 / 速度 / Q"),
        ("多模态正演", "频散 / 衰减 / 群速度"),
        ("合成能量图", "噪声 / 漏检 / 伪脊"),
        ("全局模态跟踪", "顺序 / 连续性 / 消失"),
        ("稳健局部反演", "Huber / Student-t"),
        ("信息与可信度", "宽度 / 覆盖率 / 误差"),
    ]
    body = [text(34, 38, "图1  多模态槽波信息保留与稳健拾取流程", 19, 600)]
    box_width = 145
    gap = 31
    for index, (title, subtitle) in enumerate(labels):
        x = 32 + index * (box_width + gap)
        body.append(
            f'<rect x="{x}" y="92" width="{box_width}" height="116" '
            'rx="12" fill="#ffffff" stroke="#cdd5d4"/>'
        )
        body.append(text(x + box_width / 2, 132, title, 15, 600, anchor="middle"))
        body.append(text(x + box_width / 2, 160, subtitle, 11, 400, COLORS["muted"], "middle"))
        if index < len(labels) - 1:
            body.append(
                f'<path d="M{x + box_width + 7} 150 L{x + box_width + gap - 7} 150" '
                'stroke="#8fa7a5" stroke-width="2" marker-end="url(#arrow)"/>'
            )
    body.insert(
        1,
        '<defs><marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5" '
        'markerWidth="5" markerHeight="5" orient="auto-start-reverse">'
        '<path d="M0 0 L10 5 L0 10 z" fill="#8fa7a5"/></marker></defs>',
    )
    body.append(
        text(
            550,
            263,
            "合成基准用于识别信息瓶颈；真实工区验证仍是该类方法外推的前提。",
            12,
            400,
            COLORS["muted"],
            "middle",
        )
    )
    save("fig1_workflow.svg", width, height, "".join(body))


def figure_information():
    data = json.loads((ROOT / "e13_acquisition_information_summary.json").read_text(encoding="utf-8"))
    labels = {
        "c0": "基阶",
        "c0_c1": "基阶 + 1阶",
        "c0_c1_c2": "基阶 + 1阶 + 2阶",
        "c0_q0": "基阶 + 衰减",
        "c0_c1_q": "基阶 + 1阶 + 衰减",
        "c0_c1_group_q": "多阶 + 群速度 + 衰减",
    }
    width, height = 1000, 430
    left, top = 190, 60
    plot_width, row_height = 720, 45
    maximum = max(row["median_h_width_90"] for row in data) * 1.15
    body = [text(34, 34, "\u56fe2  \u4e0d\u540c\u89c2\u6d4b\u7ec4\u5408\u7684\u539a\u5ea6\u4fe1\u606f\u9884\u7b97", 18, 600)]
    body.append(text(left, 78, "理论 90% 厚度区间中位数 / m", 11, 400, COLORS["muted"]))
    for index, row in enumerate(data):
        y = top + 40 + index * row_height
        value = row["median_h_width_90"]
        bar_width = value / maximum * plot_width
        body.append(text(left - 15, y + 5, labels[row["design"]], 12, 400, anchor="end"))
        body.append(
            f'<rect x="{left}" y="{y - 13}" width="{bar_width:.1f}" height="20" '
            f'rx="4" fill="{COLORS["teal"]}" opacity="{0.55 + index * 0.06:.2f}"/>'
        )
        body.append(text(left + bar_width + 8, y + 3, f"{value:.4f}", 11, 500))
    body.append(text(550, 404, "Fisher 信息下界；相同正演用于数据生成与反演，属于乐观估计。", 11, 400, COLORS["muted"], "middle"))
    save("fig2_information_budget.svg", width, height, "".join(body))


def figure_picking():
    data = json.loads((ROOT / "e26_additional_picker_summary.json").read_text(encoding="utf-8"))
    methods = [
        ("independent_peaks", "\u72ec\u7acb\u5cf0"),
        ("greedy_tracking", "\u8d2a\u5fc3\u8ddf\u8e2a"),
        ("viterbi_order_tracking", "\u5168\u5c40\u8ddf\u8e2a"),
        ("component_segmentation", "\u5206\u5272\u57fa\u7ebf"),
    ]
    width, height = 1000, 390
    left, top = 110, 70
    plot_width, plot_height = 760, 230
    colors = [COLORS["gray"], COLORS["gold"], COLORS["teal"], "#75549b"]
    body = [text(34, 34, "图3  基阶与高阶模态自动检出率", 18, 600)]
    for y in (0.0, 0.25, 0.5, 0.75, 1.0):
        yy = top + plot_height * (1 - y)
        body.append(f'<path d="M{left} {yy} H{left + plot_width}" stroke="{COLORS["grid"]}"/>')
        body.append(text(left - 12, yy + 4, f"{int(y*100)}%", 11, 400, COLORS["muted"], "end"))
    group_width = plot_width / len(methods)
    bar_width = 38
    for method_index, (key, label) in enumerate(methods):
        stats = data[key]["mode_detection_percent"]
        center = left + group_width * (method_index + 0.5)
        for mode in range(3):
            detection = stats[f"mode{mode}"]
            value = detection / 100.0
            x = center + (mode - 1) * (bar_width + 12) - bar_width / 2
            bar_height = value * plot_height
            body.append(
                f'<rect x="{x:.1f}" y="{top + plot_height - bar_height:.1f}" '
                f'width="{bar_width}" height="{bar_height:.1f}" rx="4" '
                f'fill="{colors[method_index]}" opacity="{0.62 + mode * 0.13:.2f}"/>'
            )
            body.append(
                text(x + bar_width / 2, top + plot_height - bar_height - 6, f"{value*100:.0f}", 10, 500, anchor="middle")
            )
        body.append(text(center, top + plot_height + 28, label, 12, 500, anchor="middle"))
        false_rate = data[key]["false_pick_percent"]
        body.append(
            text(
                center,
                top + plot_height + 50,
                f"\u4f2a\u62fe\u53d6 {false_rate:.1f}%",
                10,
                400,
                COLORS["muted"],
                "middle",
            )
        )
    for mode, label in enumerate(("\u57fa\u9636", "\u7b2c\u4e00\u9ad8\u9636", "\u7b2c\u4e8c\u9ad8\u9636")):
        x = left + plot_width - 290 + mode * 100
        body.append(f'<rect x="{x}" y="316" width="12" height="12" fill="{COLORS["teal"]}" opacity="{0.62 + mode * 0.13:.2f}"/>')
        body.append(text(x + 17, 327, label, 11))
    save("fig3_mode_picking.svg", width, height, "".join(body))


def figure_robustness():
    data = json.loads((ROOT / "e17_robust_likelihood_summary.json").read_text(encoding="utf-8"))
    width, height = 980, 430
    left, top = 90, 65
    plot_width, plot_height = 820, 280
    likelihoods = ("gaussian", "huber", "student_t")
    colors = {
        "gaussian": COLORS["gray"],
        "huber": COLORS["gold"],
        "student_t": COLORS["teal"],
    }
    labels = {"gaussian": "Gaussian", "huber": "Huber", "student_t": "Student-t"}
    series = {
        likelihood: [
            row
            for row in data
            if row["include_attenuation"] and row["likelihood"] == likelihood
        ]
        for likelihood in likelihoods
    }
    maximum = max(
        row["h_rmse"] for rows in series.values() for row in rows
    ) * 1.15
    body = [text(34, 34, "图4  自动模态拾取下的厚度稳健性", 18, 600)]
    for y in (0.0, 0.25, 0.5, 0.75, 1.0):
        yy = top + plot_height * (1 - y)
        body.append(f'<path d="M{left} {yy} H{left + plot_width}" stroke="{COLORS["grid"]}"/>')
        body.append(text(left - 12, yy + 4, f"{y * maximum:.2f}", 10, 400, COLORS["muted"], "end"))
    fractions = (0.0, 0.10, 0.20, 0.30)
    for x_index, fraction in enumerate(fractions):
        x = left + x_index / (len(fractions) - 1) * plot_width
        body.append(text(x, top + plot_height + 25, f"{fraction*100:.0f}%", 11, 400, anchor="middle"))
    for likelihood, rows in series.items():
        points = []
        for x_index, row in enumerate(rows):
            x = left + x_index / (len(fractions) - 1) * plot_width
            y = top + plot_height - row["h_rmse"] / maximum * plot_height
            points.append((x, y))
        body.append(
            '<path d="' + " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(points)) + f'" fill="none" stroke="{colors[likelihood]}" stroke-width="2.5"/>'
        )
        for x, y in points:
            body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{colors[likelihood]}"/>')
        legend_x = left + 20 + list(likelihoods).index(likelihood) * 130
        body.append(f'<circle cx="{legend_x}" cy="370" r="5" fill="{colors[likelihood]}"/>')
        body.append(text(legend_x + 12, 374, labels[likelihood], 11))
    body.append(text(left + plot_width / 2, top + plot_height + 55, "额外离群点比例", 12, 500, anchor="middle"))
    body.append(text(20, top + plot_height / 2, "厚度 RMSE / m", 12, 500, COLORS["muted"]))
    save("fig4_robust_likelihood.svg", width, height, "".join(body))


def figure_calibration():
    data = json.loads(
        (
            ROOT / "e18_robust_posterior_calibration_summary.json"
        ).read_text(encoding="utf-8")
    )
    raw = data["overall"]["raw"]["coverage_percent"]
    conformal = data["overall"]["global_conformal"]["coverage_percent"]
    width, height = 900, 430
    left, top = 90, 65
    plot_width, plot_height = 760, 285
    body = [
        text(34, 34, "Figure 5. Simulated interval coverage", 19, 600),
        text(
            left + plot_width / 2,
            402,
            "Parameter",
            12,
            500,
            anchor="middle",
        ),
        text(
            20,
            top + plot_height / 2,
            "Empirical coverage (%)",
            12,
            500,
            COLORS["muted"],
        ),
    ]
    for value in (0, 25, 50, 75, 90, 100):
        y = top + plot_height * (1 - value / 100.0)
        body.append(
            f'<path d="M{left} {y:.1f} H{left + plot_width}" '
            f'stroke="{COLORS["grid"]}" '
            f'{"stroke-dasharray=\"7 5\"" if value == 90 else ""}/>'
        )
        body.append(
            text(
                left - 12,
                y + 4,
                str(value),
                10,
                400,
                COLORS["muted"],
                "end",
            )
        )

    labels = ("Thickness", "Velocity", "Q")
    group_width = plot_width / 3
    bar_width = 55
    for index, label in enumerate(labels):
        center = left + group_width * (index + 0.5)
        for series_index, (values, fill) in enumerate(
            (
                (raw, COLORS["gray"]),
                (conformal, COLORS["teal"]),
            )
        ):
            value = values[index]
            bar_height = value / 100.0 * plot_height
            x = center + (series_index - 0.5) * (bar_width + 18) - bar_width / 2
            body.append(
                f'<rect x="{x:.1f}" y="{top + plot_height - bar_height:.1f}" '
                f'width="{bar_width}" height="{bar_height:.1f}" rx="5" '
                f'fill="{fill}" opacity="{0.72 + series_index * 0.16:.2f}"/>'
            )
            body.append(
                text(
                    x + bar_width / 2,
                    top + plot_height - bar_height - 7,
                    f"{value:.1f}",
                    11,
                    500,
                    anchor="middle",
                )
            )
        body.append(text(center, 378, label, 12, 500, anchor="middle"))

    for index, (label, fill) in enumerate(
        (("Raw local", COLORS["gray"]), ("Global conformal", COLORS["teal"]))
    ):
        x = left + 235 + index * 190
        body.append(f'<rect x="{x}" y="43" width="14" height="14" fill="{fill}"/>')
        body.append(text(x + 20, 55, label, 11))
    body.append(
        text(
            left + 610,
            55,
            "Nominal 90%",
            11,
            500,
            COLORS["gold"],
        )
    )
    save("fig5_calibration.svg", width, height, "".join(body))


def figure_sls_attenuation():
    summary = json.loads(
        (
            ROOT / "e28_sls_attenuation_validation_summary.json"
        ).read_text(encoding="utf-8")
    )
    rows = list(
        csv.DictReader(
            (
                ROOT / "e28_sls_attenuation_validation.csv"
            ).open(encoding="utf-8")
        )
    )
    centers = [float(row["offset_m"]) for row in rows]
    log_ratio = [float(row["log_ratio"]) for row in rows]
    fit_slope, fit_intercept = _linear_fit(centers, log_ratio)
    expected_slope = float(summary["expected_log_slope"])
    expected_intercept = sum(
        value - expected_slope * center
        for center, value in zip(centers, log_ratio, strict=True)
    ) / len(centers)

    width, height = 1040, 450
    left, top = 80, 72
    left_width, plot_height = 590, 300
    right_left = 735
    right_width = 245
    x_min, x_max = min(centers) - 4.0, max(centers) + 4.0
    y_min = min(log_ratio) - 0.12
    y_max = max(log_ratio) + 0.12

    def x_map(value):
        return left + (value - x_min) / (x_max - x_min) * left_width

    def y_map(value):
        return top + plot_height - (value - y_min) / (y_max - y_min) * plot_height

    body = [
        text(32, 34, "Figure 8. Viscoelastic Love-wave attenuation check", 19, 600),
        text(
            left + left_width / 2,
            415,
            "Source-receiver offset (m)",
            12,
            500,
            anchor="middle",
        ),
        text(
            18,
            top + plot_height / 2,
            "ln(A_low / A_high)",
            12,
            500,
            COLORS["muted"],
        ),
    ]
    body.append(
        f'<rect x="{left}" y="{top}" width="{left_width}" height="{plot_height}" '
        'fill="#ffffff" stroke="#b8c2c4"/>'
    )
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        value = y_min + fraction * (y_max - y_min)
        y = y_map(value)
        body.append(
            f'<path d="M{left} {y:.1f} H{left + left_width}" '
            f'stroke="{COLORS["grid"]}"/>'
        )
        body.append(
            text(
                left - 10,
                y + 4,
                f"{value:.2f}",
                10,
                400,
                COLORS["muted"],
                "end",
            )
        )
    for value in (60, 80, 100, 120, 140, 160):
        x = x_map(value)
        body.append(
            f'<path d="M{x:.1f} {top} V{top + plot_height}" '
            f'stroke="{COLORS["grid"]}"/>'
        )
        body.append(
            text(
                x,
                top + plot_height + 18,
                f"{value:.0f}",
                10,
                400,
                COLORS["muted"],
                "middle",
            )
        )
    body.append(
        '<path d="'
        + " ".join(
            f"{'M' if index == 0 else 'L'}"
            f"{x_map(center):.1f},{y_map(fit_intercept + fit_slope * center):.1f}"
            for index, center in enumerate((x_min, x_max))
        )
        + f'" fill="none" stroke="{COLORS["teal"]}" stroke-width="2.5"/>'
    )
    body.append(
        '<path d="'
        + " ".join(
            f"{'M' if index == 0 else 'L'}"
            f"{x_map(center):.1f},"
            f"{y_map(expected_intercept + expected_slope * center):.1f}"
            for index, center in enumerate((x_min, x_max))
        )
        + f'" fill="none" stroke="{COLORS["gold"]}" stroke-width="2.5" '
        'stroke-dasharray="8 5"/>'
    )
    for center, value in zip(centers, log_ratio, strict=True):
        body.append(
            f'<circle cx="{x_map(center):.1f}" cy="{y_map(value):.1f}" '
            'r="4.2" fill="#18222c"/>'
        )
    body.append(
        f'<circle cx="{left + 24}" cy="391" r="5" fill="{COLORS["teal"]}"/>'
    )
    body.append(text(left + 38, 395, "Observed linear fit", 11))
    body.append(
        f'<circle cx="{left + 205}" cy="391" r="5" fill="{COLORS["gold"]}"/>'
    )
    body.append(text(left + 219, 395, "First-order modal-Q theory", 11))

    body.append(
        text(
            right_left + right_width / 2,
            52,
            "SLS quality-factor fit",
            13,
            600,
            anchor="middle",
        )
    )
    q40 = summary["sls_fit"]["q40_min_max_mean_error_percent"]
    q300 = summary["sls_fit"]["q300_min_max_mean_error_percent"]
    groups = (
        ("Qc = 40", q40, 40.0),
        ("Qc = 300", q300, 300.0),
    )
    for index, (label, values, target) in enumerate(groups):
        center = right_left + right_width * (index + 0.5) / 2
        base_y = top + plot_height
        target_height = target / 350.0 * plot_height
        mean_height = values[2] / 350.0 * plot_height
        minimum_height = values[0] / 350.0 * plot_height
        maximum_height = values[1] / 350.0 * plot_height
        body.append(
            f'<path d="M{center - 30:.1f} {base_y} V{base_y - maximum_height:.1f}" '
            'stroke="#9ca7aa" stroke-width="2"/>'
        )
        body.append(
            f'<path d="M{center - 38:.1f} {base_y - minimum_height:.1f} '
            f'H{center - 22:.1f} M{center - 38:.1f} '
            f'{base_y - maximum_height:.1f} H{center - 22:.1f}" '
            'stroke="#9ca7aa" stroke-width="2"/>'
        )
        body.append(
            f'<rect x="{center - 12:.1f}" y="{base_y - mean_height:.1f}" '
            f'width="24" height="{mean_height:.1f}" rx="3" '
            f'fill="{COLORS["teal"]}" opacity="0.82"/>'
        )
        body.append(
            f'<path d="M{center - 30:.1f} {base_y - target_height:.1f} '
            f'H{center + 30:.1f}" stroke="{COLORS["gold"]}" '
            'stroke-width="2" stroke-dasharray="6 4"/>'
        )
        body.append(
            text(
                center,
                base_y + 20,
                label,
                11,
                500,
                anchor="middle",
            )
        )
        body.append(
            text(
                center,
                base_y - mean_height - 10,
                f"{values[2]:.0f}",
                10,
                500,
                anchor="middle",
            )
        )
    body.append(
        text(
            right_left + right_width / 2,
            391,
            "Bars: fitted mean Q; whiskers: 80-550 Hz range",
            9,
            400,
            COLORS["muted"],
            "middle",
        )
    )
    body.append(
        text(
            right_left + right_width / 2,
            408,
            "Dashed lines: input Q",
            9,
            400,
            COLORS["muted"],
            "middle",
        )
    )
    save("fig8_sls_attenuation.svg", width, height, "".join(body))


def _linear_fit(x_values, y_values):
    count = len(x_values)
    mean_x = sum(x_values) / count
    mean_y = sum(y_values) / count
    denominator = sum((value - mean_x) ** 2 for value in x_values)
    slope = sum(
        (x_value - mean_x) * (y_value - mean_y)
        for x_value, y_value in zip(x_values, y_values, strict=True)
    ) / denominator
    return slope, mean_y - slope * mean_x


def main():
    figure_workflow()
    figure_information()
    figure_picking()
    figure_robustness()
    figure_calibration()
    figure_sls_attenuation()
    print(f"figures written to {OUT}")


if __name__ == "__main__":
    main()
