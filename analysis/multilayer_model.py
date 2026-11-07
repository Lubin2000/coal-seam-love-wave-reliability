"""Elastic Love-wave forward model for piecewise layers between half-spaces.

The mode solver is used only to generate higher-fidelity synthetic truth in the
E1 model-error experiment. The inversion deliberately remains a simpler
three-layer model.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass
class Layer:
    thickness: float
    beta: float
    density: float
    q: float


@dataclass
class Halfspace:
    beta: float
    density: float
    q: float

    @property
    def mu(self) -> float:
        return self.density * self.beta**2


def shear_modulus(density: float, beta: float) -> float:
    return density * beta**2


def propagate_layer(
    layer: Layer,
    frequency: float,
    phase_velocity: float,
) -> tuple[complex, complex]:
    """Return the state transfer matrix for one finite layer."""
    omega = 2.0 * math.pi * frequency
    k = omega / phase_velocity
    mu = shear_modulus(layer.density, layer.beta)
    nu = np.sqrt(complex(omega**2 / layer.beta**2 - k**2, 0.0))
    phase = nu * layer.thickness

    if abs(nu) < 1e-10:
        return complex(1.0, 0.0), complex(layer.thickness, 0.0)

    c = np.cos(phase)
    s = np.sin(phase)
    return (
        c,
        s / (mu * nu),
    )


def propagate_state(
    layer: Layer,
    state: tuple[complex, complex],
    frequency: float,
    phase_velocity: float,
) -> tuple[complex, complex]:
    v0, t0 = state
    c, inv_impedance = propagate_layer(layer, frequency, phase_velocity)
    mu = shear_modulus(layer.density, layer.beta)
    omega = 2.0 * math.pi * frequency
    k = omega / phase_velocity
    nu = np.sqrt(complex(omega**2 / layer.beta**2 - k**2, 0.0))
    impedance = mu * nu
    v1 = c * v0 + inv_impedance * t0
    t1 = -impedance * np.sin(nu * layer.thickness) * v0 + c * t0
    return complex(v1), complex(t1)


def secular_value_general(
    layers: list[Layer],
    upper: Halfspace,
    lower: Halfspace,
    frequency: float,
    phase_velocity: float,
) -> complex:
    omega = 2.0 * math.pi * frequency
    k = omega / phase_velocity
    gamma0 = np.sqrt(complex(k**2 - omega**2 / upper.beta**2, 0.0))
    gamma2 = np.sqrt(complex(k**2 - omega**2 / lower.beta**2, 0.0))

    # At z=0 the upper half-space decays away from the layer; its traction is
    # mu * gamma * v. Set v=1 and propagate to the bottom interface.
    state = (complex(1.0, 0.0), complex(upper.mu * gamma0, 0.0))
    for layer in layers:
        state = propagate_state(layer, state, frequency, phase_velocity)
    v_bottom, t_bottom = state

    # The lower half-space decays for z > bottom; its traction is
    # -mu * gamma * v. The secular equation is t_bottom + mu2*gamma2*v=0.
    return t_bottom + lower.mu * gamma2 * v_bottom


def secular_values_general(
    layers: list[Layer],
    upper: Halfspace,
    lower: Halfspace,
    frequency: float,
    phase_velocities: np.ndarray,
) -> np.ndarray:
    """Vectorised secular function over an array of trial phase velocities."""
    omega = 2.0 * math.pi * frequency
    k = omega / phase_velocities
    gamma0 = np.sqrt(k**2 - omega**2 / upper.beta**2 + 0j)
    gamma2 = np.sqrt(k**2 - omega**2 / lower.beta**2 + 0j)
    values = np.ones_like(k, dtype=complex)
    tractions = upper.mu * gamma0

    for layer in layers:
        mu = shear_modulus(layer.density, layer.beta)
        nu = np.sqrt(omega**2 / layer.beta**2 - k**2 + 0j)
        phase = nu * layer.thickness
        cos_phase = np.cos(phase)
        sin_phase = np.sin(phase)
        inverse_impedance = np.divide(
            sin_phase,
            mu * nu,
            out=np.zeros_like(phase),
            where=np.abs(nu) > 1e-12,
        )
        new_values = cos_phase * values + inverse_impedance * tractions
        new_tractions = (
            -mu * nu * sin_phase * values + cos_phase * tractions
        )
        values = new_values
        tractions = new_tractions

    return tractions + lower.mu * gamma2 * values


def secular_values_multi(
    layers: list[Layer],
    upper: Halfspace,
    lower: Halfspace,
    frequencies: np.ndarray,
    phase_velocities: np.ndarray,
) -> np.ndarray:
    """Vectorised secular values with one frequency-velocity pair per point."""
    omega = 2.0 * math.pi * frequencies
    k = omega / phase_velocities
    gamma0 = np.sqrt(k**2 - omega**2 / upper.beta**2 + 0j)
    gamma2 = np.sqrt(k**2 - omega**2 / lower.beta**2 + 0j)
    values = np.ones_like(k, dtype=complex)
    tractions = upper.mu * gamma0

    for layer in layers:
        mu = shear_modulus(layer.density, layer.beta)
        nu = np.sqrt(omega**2 / layer.beta**2 - k**2 + 0j)
        phase = nu * layer.thickness
        cos_phase = np.cos(phase)
        sin_phase = np.sin(phase)
        inverse_impedance = np.divide(
            sin_phase,
            mu * nu,
            out=np.zeros_like(phase),
            where=np.abs(nu) > 1e-12,
        )
        new_values = cos_phase * values + inverse_impedance * tractions
        new_tractions = (
            -mu * nu * sin_phase * values + cos_phase * tractions
        )
        values = new_values
        tractions = new_tractions

    return tractions + lower.mu * gamma2 * values


def fundamental_phase_velocities_general(
    layers: list[Layer],
    upper: Halfspace,
    lower: Halfspace,
    frequencies: np.ndarray,
    grid_points: int = 520,
) -> np.ndarray:
    frequencies = np.asarray(frequencies, dtype=float)
    beta_min = min(layer.beta for layer in layers)
    beta_max = min(upper.beta, lower.beta)
    deltas = np.geomspace(
        1e-5,
        beta_max - beta_min - 1e-4,
        grid_points,
    )
    grid = beta_max - deltas

    lows = np.empty_like(frequencies)
    highs = np.empty_like(frequencies)
    low_values = np.empty_like(frequencies)
    high_values = np.empty_like(frequencies)

    for frequency_index, frequency in enumerate(frequencies):
        values = np.real(
            secular_values_general(
                layers,
                upper,
                lower,
                float(frequency),
                grid,
            )
        )
        found = False
        for index in range(len(grid) - 1, 0, -1):
            c_high = float(grid[index])
            c_low = float(grid[index - 1])
            f_high = float(values[index])
            f_low = float(values[index - 1])
            if f_high == 0.0:
                lows[frequency_index] = c_high
                highs[frequency_index] = c_high
                low_values[frequency_index] = f_high
                high_values[frequency_index] = f_high
                found = True
                break
            if f_high * f_low < 0.0:
                lows[frequency_index] = c_low
                highs[frequency_index] = c_high
                low_values[frequency_index] = f_low
                high_values[frequency_index] = f_high
                found = True
                break
        if not found:
            raise RuntimeError(
                f"No fundamental root found at {frequency} Hz; beta range "
                f"{beta_min:.2f}-{beta_max:.2f} m/s."
            )

    for _ in range(55):
        mids = 0.5 * (lows + highs)
        mid_values = np.real(
            secular_values_multi(
                layers,
                upper,
                lower,
                frequencies,
                mids,
            )
        )
        same_as_low = low_values * mid_values <= 0.0
        highs = np.where(same_as_low, mids, highs)
        high_values = np.where(same_as_low, mid_values, high_values)
        lows = np.where(same_as_low, lows, mids)
        low_values = np.where(same_as_low, low_values, mid_values)

    return 0.5 * (lows + highs)


def fundamental_phase_velocity_general(
    layers: list[Layer],
    upper: Halfspace,
    lower: Halfspace,
    frequency: float,
    grid_points: int = 520,
) -> float:
    return float(
        fundamental_phase_velocities_general(
            layers,
            upper,
            lower,
            np.asarray([frequency], dtype=float),
            grid_points=grid_points,
        )[0]
    )


def mode_energy_and_q(
    layers: list[Layer],
    upper: Halfspace,
    lower: Halfspace,
    frequency: float,
    phase_velocity: float,
) -> tuple[float, float]:
    omega = 2.0 * math.pi * frequency
    k = omega / phase_velocity
    gamma0 = np.sqrt(complex(k**2 - omega**2 / upper.beta**2, 0.0))
    gamma2 = np.sqrt(complex(k**2 - omega**2 / lower.beta**2, 0.0))

    energies: list[float] = []
    state = (complex(1.0, 0.0), complex(upper.mu * gamma0, 0.0))
    for layer in layers:
        v0, t0 = state
        mu = shear_modulus(layer.density, layer.beta)
        nu = np.sqrt(complex(omega**2 / layer.beta**2 - k**2, 0.0))
        if abs(nu) < 1e-10:
            amplitude = v0
            derivative = t0 / mu
            energy = float(
                np.real(mu * (abs(derivative) ** 2 + k**2 * abs(amplitude) ** 2))
                * layer.thickness
            )
        else:
            amplitude_0 = v0
            amplitude_1 = t0 / (mu * nu)
            z = np.linspace(0.0, layer.thickness, 301)
            amplitude = amplitude_0 * np.cos(nu * z) + amplitude_1 * np.sin(nu * z)
            derivative = (
                -amplitude_0 * nu * np.sin(nu * z)
                + amplitude_1 * nu * np.cos(nu * z)
            )
            integrand = mu * (np.abs(derivative) ** 2 + k**2 * np.abs(amplitude) ** 2)
            energy = float(np.trapezoid(np.real(integrand), z))
        energies.append(energy)
        state = propagate_state(layer, state, frequency, phase_velocity)

    v_bottom, _ = state
    upper_tail = (
        upper.mu
        * (abs(gamma0) ** 2 + k**2)
        / (2.0 * max(np.real(gamma0), 1e-12))
    )
    lower_tail = (
        lower.mu
        * (abs(gamma2) ** 2 + k**2)
        * abs(v_bottom) ** 2
        / (2.0 * max(np.real(gamma2), 1e-12))
    )

    energy_values = [float(np.real(upper_tail))] + energies + [
        float(np.real(lower_tail))
    ]
    q_values = [upper.q] + [layer.q for layer in layers] + [lower.q]
    total = sum(energy_values)
    q_inverse = sum(energy / q for energy, q in zip(energy_values, q_values)) / total
    # Return the coal-layer energy fraction specifically. For gradient models
    # this is the energy in all finite layers.
    finite_fraction = sum(energies) / total
    return float(q_inverse), float(finite_fraction)


def response_general(
    layers: list[Layer],
    upper: Halfspace,
    lower: Halfspace,
    frequencies: np.ndarray,
    grid_points: int = 520,
) -> tuple[np.ndarray, np.ndarray]:
    phase = fundamental_phase_velocities_general(
        layers,
        upper,
        lower,
        frequencies,
        grid_points=grid_points,
    )
    q_inverse = []
    for frequency, c in zip(frequencies, phase):
        qinv, _ = mode_energy_and_q(
            layers,
            upper,
            lower,
            float(frequency),
            c,
        )
        q_inverse.append(qinv)
    return np.asarray(phase), np.asarray(q_inverse)


def build_three_layer(
    h: float,
    beta_coal: float,
    q_coal: float,
    beta_roof: float = 2400.0,
    beta_floor: float = 2650.0,
    rho_coal: float = 1380.0,
    rho_roof: float = 2450.0,
    rho_floor: float = 2600.0,
    q_roof: float = 100.0,
    q_floor: float = 120.0,
) -> tuple[list[Layer], Halfspace, Halfspace]:
    return (
        [Layer(h, beta_coal, rho_coal, q_coal)],
        Halfspace(beta_roof, rho_roof, q_roof),
        Halfspace(beta_floor, rho_floor, q_floor),
    )


def smoothstep(value: float) -> float:
    x = min(max(value, 0.0), 1.0)
    return x * x * (3.0 - 2.0 * x)


def build_gradient_seam(
    h: float,
    transition_half_width: float,
    n_layers: int = 60,
    beta_coal: float = 1120.0,
    q_coal: float = 30.0,
    density_coal: float = 1380.0,
) -> tuple[list[Layer], Halfspace, Halfspace]:
    """A seam with smooth roof/seam and seam/floor transition zones."""
    upper = Halfspace(2400.0, 2450.0, 100.0)
    lower = Halfspace(2650.0, 2600.0, 120.0)
    coal_beta = beta_coal
    coal_density = density_coal
    coal_q = q_coal
    thickness = h / n_layers
    layers: list[Layer] = []

    for index in range(n_layers):
        z = (index + 0.5) * thickness
        if z < transition_half_width:
            weight = smoothstep(z / transition_half_width)
            beta = upper.beta * (1.0 - weight) + coal_beta * weight
            density = upper.density * (1.0 - weight) + coal_density * weight
            q = math.exp(
                (1.0 - weight) * math.log(upper.q) + weight * math.log(coal_q)
            )
        elif z > h - transition_half_width:
            distance = h - z
            weight = smoothstep(distance / transition_half_width)
            beta = lower.beta * (1.0 - weight) + coal_beta * weight
            density = lower.density * (1.0 - weight) + coal_density * weight
            q = math.exp(
                (1.0 - weight) * math.log(lower.q) + weight * math.log(coal_q)
            )
        else:
            beta = coal_beta
            density = coal_density
            q = coal_q
        layers.append(Layer(thickness, beta, density, q))
    return layers, upper, lower


def build_interlayer_seam(
    h: float,
    interlayer_thickness: float,
    beta_coal: float = 1120.0,
    q_coal: float = 30.0,
    density_coal: float = 1380.0,
) -> tuple[list[Layer], Halfspace, Halfspace]:
    upper = Halfspace(2400.0, 2450.0, 100.0)
    lower = Halfspace(2650.0, 2600.0, 120.0)
    coal_1 = 0.92
    coal_2 = h - coal_1 - interlayer_thickness
    layers = [
        Layer(coal_1, beta_coal, density_coal, q_coal),
        Layer(interlayer_thickness, 1750.0, 2200.0, 70.0),
        Layer(coal_2, beta_coal, density_coal, q_coal),
    ]
    return layers, upper, lower
