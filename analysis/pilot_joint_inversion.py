"""Pilot experiment: dispersion-only vs dispersion-attenuation inversion.

This is a deliberately reduced-order study for coal-seam Love waves. It uses
the exact elastic dispersion relation for a low-velocity layer between two
half-spaces and first-order perturbation theory for mode attenuation.

The goal is not to simulate a complete mine seismic survey. The goal is to test
whether attenuation supplies identifiable information beyond dispersion before
investing in a larger neural inverse model.
"""

from __future__ import annotations

import csv
import json
import math
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
RNG_SEED = 20261105
N_REPLICATES = 36

TRUE_MODEL = {
    "h": 2.40,
    "beta_coal": 1120.0,
    "q_coal": 30.0,
    "beta_roof": 2400.0,
    "beta_floor": 2650.0,
    "rho_coal": 1380.0,
    "rho_roof": 2450.0,
    "rho_floor": 2600.0,
    "q_roof": 100.0,
    "q_floor": 120.0,
}

START_MODEL = {
    "h": 2.80,
    "beta_coal": 1050.0,
    "q_coal": 45.0,
}

LOWER = np.array([1.20, 850.0, math.log(5.0)])
UPPER = np.array([4.20, 1400.0, math.log(100.0)])
FREQUENCIES = np.linspace(100.0, 850.0, 48)


def shear_modulus(density: float, velocity: float) -> float:
    return density * velocity**2


def secular_value(
    phase_velocity: float,
    frequency: float,
    h: float,
    beta1: float,
) -> float:
    """Fundamental-mode Love-wave secular function in angle form."""
    omega = 2.0 * math.pi * frequency
    beta0 = TRUE_MODEL["beta_roof"]
    beta2 = TRUE_MODEL["beta_floor"]
    mu0 = shear_modulus(TRUE_MODEL["rho_roof"], beta0)
    mu1 = shear_modulus(TRUE_MODEL["rho_coal"], beta1)
    mu2 = shear_modulus(TRUE_MODEL["rho_floor"], beta2)

    if not beta1 < phase_velocity < min(beta0, beta2):
        return math.nan

    k = omega / phase_velocity
    nu = math.sqrt(max(omega**2 / beta1**2 - k**2, 0.0))
    gamma0 = math.sqrt(max(k**2 - omega**2 / beta0**2, 0.0))
    gamma2 = math.sqrt(max(k**2 - omega**2 / beta2**2, 0.0))
    if nu <= 1e-10:
        return -math.pi

    r = mu2 * gamma2 / (mu1 * nu)
    s = mu0 * gamma0 / (mu1 * nu)
    return nu * h - (math.atan(r) + math.atan(s))


def fundamental_phase_velocity(
    frequency: float,
    h: float,
    beta1: float,
) -> float:
    """Find the fundamental Love-wave root by stable bisection."""
    low = beta1 + 1e-4
    high = min(TRUE_MODEL["beta_roof"], TRUE_MODEL["beta_floor"]) - 1e-4

    flow = secular_value(low, frequency, h, beta1)
    fhigh = secular_value(high, frequency, h, beta1)
    if not math.isfinite(flow) or not math.isfinite(fhigh) or flow * fhigh > 0:
        # The fundamental root may approach the upper bound at very low
        # frequencies. Move inward while preserving the physical velocity band.
        for shift in (1.0, 10.0, 50.0, 100.0):
            high_test = high - shift
            fhigh = secular_value(high_test, frequency, h, beta1)
            if math.isfinite(fhigh) and flow * fhigh <= 0:
                high = high_test
                break
        else:
            raise RuntimeError(
                f"No bracketed fundamental Love root at f={frequency}, h={h}, beta={beta1}"
            )

    for _ in range(70):
        mid = 0.5 * (low + high)
        fm = secular_value(mid, frequency, h, beta1)
        if flow * fm <= 0:
            high = mid
            fhigh = fm
        else:
            low = mid
            flow = fm
    return 0.5 * (low + high)


def modal_energy_fractions(
    frequency: float,
    phase_velocity: float,
    h: float,
    beta1: float,
) -> tuple[float, float, float]:
    """Energy fractions in roof, coal and floor from the Love eigenfunction."""
    omega = 2.0 * math.pi * frequency
    beta0 = TRUE_MODEL["beta_roof"]
    beta2 = TRUE_MODEL["beta_floor"]
    rho0 = TRUE_MODEL["rho_roof"]
    rho1 = TRUE_MODEL["rho_coal"]
    rho2 = TRUE_MODEL["rho_floor"]
    mu0 = shear_modulus(rho0, beta0)
    mu1 = shear_modulus(rho1, beta1)
    mu2 = shear_modulus(rho2, beta2)

    k = omega / phase_velocity
    nu = math.sqrt(max(omega**2 / beta1**2 - k**2, 0.0))
    gamma0 = math.sqrt(max(k**2 - omega**2 / beta0**2, 0.0))
    gamma2 = math.sqrt(max(k**2 - omega**2 / beta2**2, 0.0))
    r = mu2 * gamma2 / (mu1 * nu) if nu > 1e-12 else 1.0
    s = mu0 * gamma0 / (mu1 * nu) if nu > 1e-12 else 1.0

    z_layer = np.linspace(0.0, h, 801)
    v_layer = np.cos(nu * z_layer) + r * np.sin(nu * z_layer)
    dv_layer = -nu * np.sin(nu * z_layer) + r * nu * np.cos(nu * z_layer)
    energy_density_layer = mu1 * (dv_layer**2 + k**2 * v_layer**2)
    energy_coal = float(np.trapezoid(energy_density_layer, z_layer))

    # Integrate the evanescent tails analytically from the boundary eigenfunction.
    v_roof_interface = math.cos(nu * h) + r * math.sin(nu * h)
    v_floor_interface = 1.0
    tail_roof = (gamma0**2 + k**2) * v_roof_interface**2 / (2.0 * gamma0)
    tail_floor = (gamma2**2 + k**2) * v_floor_interface**2 / (2.0 * gamma2)
    energy_roof = mu0 * tail_roof
    energy_floor = mu2 * tail_floor

    total = energy_roof + energy_coal + energy_floor
    return energy_roof / total, energy_coal / total, energy_floor / total


def forward_response(
    frequencies: np.ndarray,
    h: float,
    beta1: float,
    q_coal: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    phase = []
    q_inverse = []
    coal_fraction = []
    q0 = 1.0 / TRUE_MODEL["q_roof"]
    q1 = 1.0 / q_coal
    q2 = 1.0 / TRUE_MODEL["q_floor"]

    for frequency in frequencies:
        try:
            c = fundamental_phase_velocity(float(frequency), h, beta1)
        except RuntimeError:
            # Some trial parameter combinations do not support a trapped
            # fundamental mode. Return an invalid response and let the
            # optimiser reject that trial rather than terminating the run.
            phase.append(math.nan)
            q_inverse.append(math.nan)
            coal_fraction.append(math.nan)
            continue
        e0, e1, e2 = modal_energy_fractions(float(frequency), c, h, beta1)
        qeff_inv = e0 * q0 + e1 * q1 + e2 * q2
        phase.append(c)
        q_inverse.append(qeff_inv)
        coal_fraction.append(e1)

    return np.asarray(phase), np.asarray(q_inverse), np.asarray(coal_fraction)


def make_observations(rng: np.random.Generator) -> dict[str, np.ndarray]:
    phase, q_inverse, _ = forward_response(
        FREQUENCIES,
        TRUE_MODEL["h"],
        TRUE_MODEL["beta_coal"],
        TRUE_MODEL["q_coal"],
    )
    phase_sigma = 0.006 * phase
    q_sigma = 0.06 * q_inverse
    return {
        "phase": phase + rng.normal(0.0, phase_sigma),
        "phase_sigma": phase_sigma,
        "q_inverse": q_inverse + rng.normal(0.0, q_sigma),
        "q_sigma": q_sigma,
    }


def unpack_joint(x: np.ndarray) -> tuple[float, float, float]:
    return float(x[0]), float(x[1]), float(math.exp(x[2]))


def joint_residual(x: np.ndarray, observations: dict[str, np.ndarray]) -> np.ndarray:
    h, beta, q = unpack_joint(x)
    phase_pred, q_pred, _ = forward_response(FREQUENCIES, h, beta, q)
    residual = np.concatenate(
        (
            (phase_pred - observations["phase"]) / observations["phase_sigma"],
            (q_pred - observations["q_inverse"]) / observations["q_sigma"],
        )
    )
    return np.nan_to_num(residual, nan=1e6, posinf=1e6, neginf=-1e6)


def dispersion_residual(x: np.ndarray, observations: dict[str, np.ndarray]) -> np.ndarray:
    h, beta = float(x[0]), float(x[1])
    phase_pred, _, _ = forward_response(
        FREQUENCIES,
        h,
        beta,
        START_MODEL["q_coal"],
    )
    residual = (phase_pred - observations["phase"]) / observations["phase_sigma"]
    return np.nan_to_num(residual, nan=1e6, posinf=1e6, neginf=-1e6)


def finite_difference_jacobian(
    function,
    x: np.ndarray,
    observations: dict[str, np.ndarray],
    steps: np.ndarray,
) -> np.ndarray:
    columns = []
    for index, step in enumerate(steps):
        xp = x.copy()
        xm = x.copy()
        xp[index] = min(xp[index] + step, UPPER[index] if index < len(UPPER) else xp[index] + step)
        xm[index] = max(xm[index] - step, LOWER[index] if index < len(LOWER) else xm[index] - step)
        actual = xp[index] - xm[index]
        columns.append(
            (function(xp, observations) - function(xm, observations)) / actual
        )
    return np.column_stack(columns)


def optimize(
    residual_function,
    initial: np.ndarray,
    observations: dict[str, np.ndarray],
    parameter_steps: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    x = initial.copy()
    lower = LOWER[: len(x)]
    upper = UPPER[: len(x)]
    residual = residual_function(x, observations)
    objective = float(residual @ residual)

    for _ in range(30):
        jacobian = finite_difference_jacobian(
            residual_function,
            x,
            observations,
            parameter_steps,
        )
        normal = jacobian.T @ jacobian
        gradient = jacobian.T @ residual
        accepted = False

        for damping in (1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0):
            system = normal + damping * np.diag(
                np.maximum(np.diag(normal), 1e-8)
            )
            try:
                delta = np.linalg.solve(system, -gradient)
            except np.linalg.LinAlgError:
                continue

            # Bound-aware line search.
            max_scale = 1.0
            for idx in range(len(x)):
                if delta[idx] > 0:
                    max_scale = min(
                        max_scale,
                        max((upper[idx] - x[idx]) / delta[idx], 0.0),
                    )
                elif delta[idx] < 0:
                    max_scale = min(
                        max_scale,
                        max((lower[idx] - x[idx]) / delta[idx], 0.0),
                    )
            delta *= max_scale

            for scale in (1.0, 0.5, 0.25, 0.1):
                trial = np.clip(x + scale * delta, lower, upper)
                trial_residual = residual_function(trial, observations)
                trial_objective = float(trial_residual @ trial_residual)
                if trial_objective < objective:
                    x = trial
                    residual = trial_residual
                    objective = trial_objective
                    accepted = True
                    break
            if accepted:
                break

        if not accepted or np.linalg.norm(delta) < 1e-7:
            break

    jacobian = finite_difference_jacobian(
        residual_function,
        x,
        observations,
        parameter_steps,
    )
    return x, jacobian, objective


def covariance_from_jacobian(jacobian: np.ndarray) -> np.ndarray:
    normal = jacobian.T @ jacobian
    return np.linalg.pinv(normal)


def normalized_condition(jacobian: np.ndarray) -> float:
    norms = np.linalg.norm(jacobian, axis=0)
    normalized = jacobian / np.maximum(norms, 1e-12)
    singular_values = np.linalg.svd(normalized, compute_uv=False)
    positive = singular_values[singular_values > 1e-10]
    if positive.size < 2:
        return math.inf
    return float(positive[0] / positive[-1])


def run_replicates() -> list[dict[str, float | int | str | bool]]:
    rng = np.random.default_rng(RNG_SEED)
    rows: list[dict[str, float | int | str | bool]] = []
    initial_dispersion = np.array(
        [START_MODEL["h"], START_MODEL["beta_coal"]], dtype=float
    )
    initial_joint = np.array(
        [START_MODEL["h"], START_MODEL["beta_coal"], math.log(START_MODEL["q_coal"])]
    )
    start = time.time()

    for replicate in range(N_REPLICATES):
        observations = make_observations(rng)

        x_disp, jac_disp, obj_disp = optimize(
            dispersion_residual,
            initial_dispersion,
            observations,
            np.array([0.005, 2.0]),
        )
        cov_disp = covariance_from_jacobian(jac_disp)
        std_disp = np.sqrt(np.maximum(np.diag(cov_disp), 0.0))

        x_joint, jac_joint, obj_joint = optimize(
            joint_residual,
            initial_joint,
            observations,
            np.array([0.005, 2.0, 0.01]),
        )
        h_joint, beta_joint, q_joint = unpack_joint(x_joint)
        cov_joint = covariance_from_jacobian(jac_joint)
        std_joint = np.sqrt(np.maximum(np.diag(cov_joint), 0.0))
        q_joint_std = q_joint * std_joint[2]

        z90 = 1.6448536269514722
        row = {
            "replicate": replicate + 1,
            "method": "dispersion",
            "h_true": TRUE_MODEL["h"],
            "h_est": float(x_disp[0]),
            "h_std": float(std_disp[0]),
            "h_covered90": bool(
                abs(x_disp[0] - TRUE_MODEL["h"]) <= z90 * std_disp[0]
            ),
            "beta_true": TRUE_MODEL["beta_coal"],
            "beta_est": float(x_disp[1]),
            "beta_std": float(std_disp[1]),
            "beta_covered90": bool(
                abs(x_disp[1] - TRUE_MODEL["beta_coal"]) <= z90 * std_disp[1]
            ),
            "q_true": TRUE_MODEL["q_coal"],
            "q_est": math.nan,
            "q_std": math.nan,
            "q_covered90": False,
            "condition": normalized_condition(jac_disp),
            "objective": obj_disp,
        }
        rows.append(row)

        rows.append(
            {
                "replicate": replicate + 1,
                "method": "joint",
                "h_true": TRUE_MODEL["h"],
                "h_est": h_joint,
                "h_std": float(std_joint[0]),
                "h_covered90": bool(
                    abs(h_joint - TRUE_MODEL["h"]) <= z90 * std_joint[0]
                ),
                "beta_true": TRUE_MODEL["beta_coal"],
                "beta_est": beta_joint,
                "beta_std": float(std_joint[1]),
                "beta_covered90": bool(
                    abs(beta_joint - TRUE_MODEL["beta_coal"]) <= z90 * std_joint[1]
                ),
                "q_true": TRUE_MODEL["q_coal"],
                "q_est": q_joint,
                "q_std": q_joint_std,
                "q_covered90": bool(
                    abs(q_joint - TRUE_MODEL["q_coal"]) <= z90 * q_joint_std
                ),
                "condition": normalized_condition(jac_joint),
                "objective": obj_joint,
            }
        )

        if (replicate + 1) % 6 == 0:
            print(
                f"replicate {replicate + 1}/{N_REPLICATES} "
                f"elapsed={time.time() - start:.1f}s"
            )
    return rows


def summarize(rows: list[dict]) -> dict:
    result: dict[str, dict] = {
        "settings": {
            "seed": RNG_SEED,
            "replicates": N_REPLICATES,
            "frequencies_hz": [float(FREQUENCIES[0]), float(FREQUENCIES[-1])],
            "frequency_points": int(FREQUENCIES.size),
            "phase_noise_relative_std": 0.006,
            "attenuation_noise_relative_std": 0.06,
            "true_model": TRUE_MODEL,
            "start_model": START_MODEL,
        }
    }
    for method in ("dispersion", "joint"):
        subset = [row for row in rows if row["method"] == method]
        parameters = ("h", "beta") if method == "dispersion" else ("h", "beta", "q")
        metrics = {}
        for parameter in parameters:
            true_value = float(subset[0][f"{parameter}_true"])
            estimates = np.array([row[f"{parameter}_est"] for row in subset])
            stds = np.array([row[f"{parameter}_std"] for row in subset])
            errors = estimates - true_value
            metrics[parameter] = {
                "bias": float(np.mean(errors)),
                "rmse": float(np.sqrt(np.mean(errors**2))),
                "mae_relative_percent": float(
                    np.mean(np.abs(errors) / true_value) * 100.0
                ),
                "mean_estimated_std": float(np.mean(stds)),
                "coverage90_percent": float(
                    np.mean([row[f"{parameter}_covered90"] for row in subset]) * 100.0
                ),
            }
        metrics["median_normalized_condition"] = float(
            np.median([row["condition"] for row in subset])
        )
        metrics["mean_objective"] = float(
            np.mean([row["objective"] for row in subset])
        )
        result[method] = metrics
    return result


def write_csv(rows: list[dict]) -> None:
    path = ROOT / "pilot_results.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def svg_line(points: list[tuple[float, float]]) -> str:
    return " ".join(
        f"{'M' if index == 0 else 'L'}{x:.2f},{y:.2f}"
        for index, (x, y) in enumerate(points)
    )


def write_figure(rows: list[dict]) -> None:
    phase, q_inverse, coal_fraction = forward_response(
        FREQUENCIES,
        TRUE_MODEL["h"],
        TRUE_MODEL["beta_coal"],
        TRUE_MODEL["q_coal"],
    )
    width, height = 1160, 760
    left, top = 70, 60
    panel_w, panel_h = 490, 270

    def map_values(
        xs: np.ndarray,
        ys: np.ndarray,
        panel_x: float,
        xlim: tuple[float, float],
        ylim: tuple[float, float],
    ) -> list[tuple[float, float]]:
        return [
            (
                panel_x + (x - xlim[0]) / (xlim[1] - xlim[0]) * panel_w,
                top + panel_h - (y - ylim[0]) / (ylim[1] - ylim[0]) * panel_h,
            )
            for x, y in zip(xs, ys)
        ]

    phase_points = map_values(
        FREQUENCIES, phase, left, (100, 850), (1050, 2450)
    )
    q_points = map_values(
        FREQUENCIES,
        q_inverse * 1000.0,
        left + 580,
        (100, 850),
        (0, max(q_inverse * 1000.0) * 1.15),
    )

    def dispersion_scatter(method: str, panel_x: float) -> str:
        subset = [row for row in rows if row["method"] == method]
        return "".join(
            (
                f'<circle cx="{panel_x + 10 + (float(row["replicate"]) - 1) / (N_REPLICATES - 1) * 470:.1f}" '
                f'cy="{535 + (float(row["h_est"]) - 1.9) * 180:.1f}" r="3.2" '
                f'fill="{("#196f72" if method == "joint" else "#5f6b7a")}" opacity="0.8"/>'
            )
            for row in subset
        )

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">
<title id="title">槽波频散—衰减联合反演初步合成实验</title>
<desc id="desc">左侧为基阶 Love 波相速度曲线，右侧为模式衰减曲线，下方比较仅频散与联合反演得到的厚度。</desc>
<rect width="100%" height="100%" fill="#f7f8f6"/>
<style>
  text {{ font-family: "Microsoft YaHei", "PingFang SC", Arial, sans-serif; fill: #18222c; }}
  .small {{ font-size: 12px; fill: #6a7480; }}
  .axis {{ stroke: #c8ced3; stroke-width: 1; }}
  .grid {{ stroke: #e3e7e6; stroke-width: 1; }}
  .line {{ fill: none; stroke: #196f72; stroke-width: 2.6; }}
  .qline {{ fill: none; stroke: #c27a32; stroke-width: 2.6; }}
</style>
<text x="70" y="32" font-size="21" font-weight="600">物理基础验证：频散与衰减提供了不同但相关的地下信息</text>
<text x="70" y="52" class="small">三层 Love 槽波模型 · h = 2.40 m · βcoal = 1,120 m/s · Qcoal = 30 · 48 个频率点</text>

<text x="{left}" y="96" font-size="15" font-weight="600">A. 基阶频散曲线</text>
<rect x="{left}" y="{top + 45}" width="{panel_w}" height="{panel_h}" fill="none" class="axis"/>
<g transform="translate(0,45)">
  <path d="{svg_line([(left + panel_w * i / 5, top) for i in range(6)])}" class="grid"/>
</g>
<g transform="translate(0,45)">
  <path d="{svg_line(phase_points)}" class="line"/>
</g>
<text x="{left}" y="{top + panel_h + 70}" class="small">100 Hz</text>
<text x="{left + panel_w - 38}" y="{top + panel_h + 70}" class="small">850 Hz</text>
<text x="{left + panel_w / 2 - 40}" y="{top + panel_h + 92}" class="small">频率 / Hz</text>
<text x="{left - 53}" y="{top + 55}" class="small">2,400</text>
<text x="{left - 53}" y="{top + panel_h + 48}" class="small">1,120</text>
<text x="{left - 22}" y="{top + panel_h + 72}" class="small">m/s</text>

<text x="{left + 580}" y="96" font-size="15" font-weight="600">B. 基阶模式衰减</text>
<rect x="{left + 580}" y="{top + 45}" width="{panel_w}" height="{panel_h}" fill="none" class="axis"/>
<g transform="translate(0,45)">
  <path d="{svg_line(q_points)}" class="qline"/>
</g>
<text x="{left + 580}" y="{top + panel_h + 70}" class="small">100 Hz</text>
<text x="{left + 580 + panel_w - 38}" y="{top + panel_h + 70}" class="small">850 Hz</text>
<text x="{left + 580 + panel_w / 2 - 40}" y="{top + panel_h + 92}" class="small">频率 / Hz</text>
<text x="{left + 555}" y="{top + 55}" class="small">×10⁻³</text>
<text x="{left + 580 - 27}" y="{top + panel_h + 72}" class="small">Q⁻¹</text>

<text x="{left}" y="496" font-size="15" font-weight="600">C. 36 次含噪反演得到的煤层厚度</text>
<line x1="{left}" y1="{535 + (TRUE_MODEL['h'] - 1.9) * 180:.1f}" x2="{left + 480}" y2="{535 + (TRUE_MODEL['h'] - 1.9) * 180:.1f}" stroke="#c9a24a" stroke-width="1.5" stroke-dasharray="5 5"/>
<text x="{left + 486}" y="{535 + (TRUE_MODEL['h'] - 1.9) * 180 + 4:.1f}" class="small">真实 2.40 m</text>
{dispersion_scatter("dispersion", left)}
{dispersion_scatter("joint", left + 580)}
<text x="{left + 176}" y="725" class="small" text-anchor="middle">● 仅频散反演</text>
<text x="{left + 580 + 176}" y="725" class="small" text-anchor="middle" fill="#196f72">● 频散—衰减联合反演</text>
<text x="{left + 580 - 2}" y="496" font-size="15" font-weight="600">D. 联合反演的厚度估计</text>
<text x="{left + 580}" y="{top + panel_h + 92}" class="small">频率 / Hz</text>
<text x="70" y="750" class="small">合成基线仅用于可辨识性判断；后续需加入复杂 Q、各向异性、道间相干性与现场验证。</text>
</svg>"""
    (ROOT / "pilot_results.svg").write_text(svg, encoding="utf-8")


def main() -> None:
    rows = run_replicates()
    summary = summarize(rows)
    write_csv(rows)
    write_figure(rows)
    (ROOT / "pilot_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
