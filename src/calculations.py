from __future__ import annotations

from dataclasses import dataclass, asdict
from math import log, pi, sqrt
from typing import Dict, List


@dataclass
class InputData:
    feed_flow_kmol_h: float
    feed_composition_light: float
    distillate_composition_light: float
    bottoms_composition_light: float
    relative_volatility: float
    reflux_ratio_factor_to_rmin: float
    feed_thermal_condition_q: float
    vapor_density_kg_m3: float
    souders_brown_factor_m_s: float


@dataclass
class StagePoint:
    stage: int
    x: float
    y: float


@dataclass
class CalculationResult:
    inputs: Dict[str, float]
    distillate_flow_kmol_h: float
    bottoms_flow_kmol_h: float
    r_min: float
    r_work: float
    n_min_fenske: float
    n_theoretical: int
    stage_profile: List[StagePoint]
    vapor_flow_kmol_h: float
    vapor_vol_flow_m3_s: float
    column_diameter_m: float


def validate_mole_fraction(value: float, name: str) -> None:
    if not (0.0 < value < 1.0):
        raise ValueError(f"{name} должен быть в диапазоне (0, 1). Получено: {value}")


def y_equilibrium(x: float, alpha: float) -> float:
    return alpha * x / (1.0 + (alpha - 1.0) * x)


def rectifying_line(x: float, r: float, x_d: float) -> float:
    return (r / (r + 1.0)) * x + x_d / (r + 1.0)


def q_line_y(x: float, q: float, z_f: float) -> float:
    if abs(q - 1.0) < 1e-9:
        # Вертикальная линия; в y(x) форме не определена.
        raise ValueError("Для q=1 используйте отдельную обработку пересечения q-линии.")
    return (q / (q - 1.0)) * x - z_f / (q - 1.0)


def calc_intersection_feed_rectifying(z_f: float, q: float, r: float, x_d: float) -> tuple[float, float]:
    if abs(q - 1.0) < 1e-9:
        x_i = z_f
        y_i = rectifying_line(x_i, r, x_d)
        return x_i, y_i

    # q/(q-1)*x - zF/(q-1) = R/(R+1)*x + xD/(R+1)
    a1 = q / (q - 1.0)
    b1 = -z_f / (q - 1.0)
    a2 = r / (r + 1.0)
    b2 = x_d / (r + 1.0)
    x_i = (b2 - b1) / (a1 - a2)
    y_i = a2 * x_i + b2
    return x_i, y_i


def stripping_line(x: float, x_i: float, y_i: float, x_b: float) -> float:
    m = (y_i - x_b) / (x_i - x_b)
    return m * (x - x_b) + x_b


def invert_equilibrium_y_to_x(y: float, alpha: float) -> float:
    # y = alpha*x / (1 + (alpha-1)x)
    # => x = y / (alpha - y*(alpha-1))
    den = alpha - y * (alpha - 1.0)
    if den <= 0:
        raise ValueError("Некорректные параметры для обращения равновесной кривой.")
    return y / den


def calc_stages_mccabe_thiele(
    x_d: float,
    x_b: float,
    z_f: float,
    q: float,
    alpha: float,
    r: float,
    max_iter: int = 500,
) -> List[StagePoint]:
    x_i, y_i = calc_intersection_feed_rectifying(z_f, q, r, x_d)

    points: List[StagePoint] = []
    y_curr = x_d
    stage = 0

    while stage < max_iter:
        x_eq = invert_equilibrium_y_to_x(y_curr, alpha)
        stage += 1
        points.append(StagePoint(stage=stage, x=x_eq, y=y_curr))

        if x_eq <= x_b:
            break

        if x_eq >= x_i:
            y_next = rectifying_line(x_eq, r, x_d)
        else:
            y_next = stripping_line(x_eq, x_i, y_i, x_b)

        y_curr = y_next

        if y_curr < 0 or y_curr > 1:
            raise ValueError("Операционная линия вышла за физически допустимые границы (0..1).")

    return points


def calculate(data: InputData) -> CalculationResult:
    validate_mole_fraction(data.feed_composition_light, "feed_composition_light")
    validate_mole_fraction(data.distillate_composition_light, "distillate_composition_light")
    validate_mole_fraction(data.bottoms_composition_light, "bottoms_composition_light")

    if data.distillate_composition_light <= data.bottoms_composition_light:
        raise ValueError("xD должно быть больше xB для ректификации легколетучего компонента.")

    f = data.feed_flow_kmol_h
    z_f = data.feed_composition_light
    x_d = data.distillate_composition_light
    x_b = data.bottoms_composition_light
    alpha = data.relative_volatility

    # Баланс
    d = f * (z_f - x_b) / (x_d - x_b)
    b = f - d

    # Приближенный Rmin для бинарной смеси при q=1 (или близко к нему)
    y_f_star = y_equilibrium(z_f, alpha)
    r_min = max((x_d - y_f_star) / (y_f_star - z_f), 0.01)
    r_work = data.reflux_ratio_factor_to_rmin * r_min

    # Фенске
    n_min = log((x_d / (1.0 - x_d)) * ((1.0 - x_b) / x_b)) / log(alpha)

    # МакКэб-Тили
    stage_points = calc_stages_mccabe_thiele(
        x_d=x_d,
        x_b=x_b,
        z_f=z_f,
        q=data.feed_thermal_condition_q,
        alpha=alpha,
        r=r_work,
    )
    n_theoretical = len(stage_points)

    # Упрощенная гидравлика и диаметр
    vapor_flow_kmol_h = (r_work + 1.0) * d
    vapor_vol_flow_m3_s = (vapor_flow_kmol_h / 3600.0) * (1.0 / data.vapor_density_kg_m3)
    diameter = sqrt((4.0 * vapor_vol_flow_m3_s) / (pi * data.souders_brown_factor_m_s))

    return CalculationResult(
        inputs=asdict(data),
        distillate_flow_kmol_h=d,
        bottoms_flow_kmol_h=b,
        r_min=r_min,
        r_work=r_work,
        n_min_fenske=n_min,
        n_theoretical=n_theoretical,
        stage_profile=stage_points,
        vapor_flow_kmol_h=vapor_flow_kmol_h,
        vapor_vol_flow_m3_s=vapor_vol_flow_m3_s,
        column_diameter_m=diameter,
    )
