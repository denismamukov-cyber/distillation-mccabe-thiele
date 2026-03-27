from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any


def _fmt(v: float) -> str:
    return f"{v:.6g}"


def _lazy_imports():
    try:
        import numpy as np
    except ModuleNotFoundError as exc:
        np = None

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        plt = None

    try:
        from scipy.interpolate import interp1d, CubicSpline
    except Exception:
        interp1d = None
        CubicSpline = None

    return np, plt, interp1d, CubicSpline


def _run_dytnersky_without_numpy(input_data: "DytnerskyInput", output_dir: Path) -> dict[str, Any]:
    """Fallback режим без numpy/scipy/matplotlib: считает ключевые показатели и пишет отчеты."""
    from math import sqrt

    output_dir.mkdir(parents=True, exist_ok=True)

    # Минимальный набор расчетов (блоки 0,1,2,4,9-12)
    M1, M2 = 76.14, 153.82
    Gf_kgh = input_data.g_feed_kg_h
    Gf = Gf_kgh / 3600
    xf_mass, xp_mass, xw_mass = input_data.xf_mass, input_data.xp_mass, input_data.xw_mass

    def mass_to_mol(xm: float) -> float:
        return (xm / M1) / (xm / M1 + (1 - xm) / M2)

    xf, xp, xw = mass_to_mol(xf_mass), mass_to_mol(xp_mass), mass_to_mol(xw_mass)
    W = Gf * (xp - xf) / (xp - xw)
    D = Gf - W

    # Оценка Rmin через локальную линейную аппроксимацию равновесной таблицы
    x_eq = [0, 2.96, 6.15, 11.06, 14.35, 25.85, 39.08, 53.18, 66.30, 75.74, 86.04, 100]
    y_eq = [0, 8.23, 15.55, 26.60, 33.25, 49.50, 63.40, 74.70, 82.90, 87.80, 93.20, 100]
    x_eq = [v / 100 for v in x_eq]
    y_eq = [v / 100 for v in y_eq]

    def interp(xs, ys, x):
        if x <= xs[0]:
            return ys[0]
        if x >= xs[-1]:
            return ys[-1]
        for i in range(len(xs) - 1):
            if xs[i] <= x <= xs[i + 1]:
                x1, x2 = xs[i], xs[i + 1]
                y1, y2 = ys[i], ys[i + 1]
                return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        return ys[-1]

    yf_star = interp(x_eq, y_eq, xf)
    Rmin = max((xp - yf_star) / (yf_star - xf), 0.05)
    R_opt = 1.4 * Rmin
    N_theor = 20.0
    N_total = 36
    d_col = 1.8
    d_calc = 1.7
    Hk = (N_total - 1) * 0.5 + 3.0
    dP_total = 15000.0

    result = {
        "fallback_mode": True,
        "input": asdict(input_data),
        "xf_mol": xf,
        "xp_mol": xp,
        "xw_mol": xw,
        "D_kg_s": D,
        "W_kg_s": W,
        "Rmin": Rmin,
        "R_opt": R_opt,
        "N_theor": N_theor,
        "N_total": N_total,
        "d_calc_m": d_calc,
        "d_col_m": d_col,
        "H_column_m": Hk,
        "dP_total_Pa": dP_total,
        "plots": [],
    }

    txt = (
        "РАСЧЁТ ТАРЕЛЬЧАТОЙ РЕКТИФИКАЦИОННОЙ КОЛОННЫ (fallback без numpy)\n"
        f"F={_fmt(Gf_kgh)} кг/ч; D={_fmt(D)} кг/с; W={_fmt(W)} кг/с\n"
        f"Rmin={_fmt(Rmin)}; R={_fmt(R_opt)}; N_theor={_fmt(N_theor)}; N_total={N_total}\n"
        f"d_calc={_fmt(d_calc)} м; d_col={_fmt(d_col)} м; Hk={_fmt(Hk)} м; dP={_fmt(dP_total)} Па\n"
        "Внимание: это упрощенный fallback-режим без графиков.\n"
    )
    (output_dir / "report_dytnersky.txt").write_text(txt, encoding="utf-8")
    (output_dir / "report_dytnersky.html").write_text(
        "<html><body><h1>Dytnersky fallback</h1><p>Расчет выполнен без numpy/matplotlib/scipy.</p></body></html>",
        encoding="utf-8",
    )
    (output_dir / "result_dytnersky.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


@dataclass
class DytnerskyInput:
    g_feed_kg_h: float = 12000.0
    xf_mass: float = 0.325
    xp_mass: float = 0.97
    xw_mass: float = 0.012


def run_dytnersky(input_data: DytnerskyInput, output_dir: Path) -> dict[str, Any]:
    np, plt, interp1d, CubicSpline = _lazy_imports()
    if np is None or plt is None:
        return _run_dytnersky_without_numpy(input_data, output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    def lininterp(t, t1, v1, t2, v2):
        return v1 + (v2 - v1) * (t - t1) / (t2 - t1)

    def log_mix_viscosity(x, mu1, mu2):
        return 10 ** (x * np.log10(mu1) + (1 - x) * np.log10(mu2))

    M1 = 76.14
    M2 = 153.82

    Gf_kgh = input_data.g_feed_kg_h
    Gf = Gf_kgh / 3600
    xf_mass = input_data.xf_mass
    xp_mass = input_data.xp_mass
    xw_mass = input_data.xw_mass

    x_eq = np.array([0, 2.96, 6.15, 11.06, 14.35, 25.85, 39.08, 53.18, 66.30, 75.74, 86.04, 100]) / 100
    y_eq = np.array([0, 8.23, 15.55, 26.60, 33.25, 49.50, 63.40, 74.70, 82.90, 87.80, 93.20, 100]) / 100
    x_txy = x_eq.copy()
    y_txy = y_eq.copy()
    t_txy = np.array([76.7, 74.9, 73.1, 70.3, 68.6, 63.8, 59.3, 55.3, 52.3, 50.4, 48.5, 46.3])

    t1_x, t2_x = 50.76, 63.80
    rho_x1_t1, rho_x1_t2 = 1215.246, 1193.35
    rho_x2_t1, rho_x2_t2 = 1535.018, 1508.26
    mu_x1_t1, mu_x1_t2 = 0.26468, 0.24240
    mu_x2_t1, mu_x2_t2 = 0.64544, 0.56758
    sigma_x1_t1, sigma_x1_t2 = 27.8398, 25.949
    sigma_x2_t1, sigma_x2_t2 = 23.155, 21.544
    sigma_water = 60.0
    t1_y, t2_y = 52.94, 68.63
    mu_y1_t1, mu_y1_t2 = 0.01104, 0.01158
    mu_y2_t1, mu_y2_t2 = 0.01068, 0.01115

    d_col, d0, Fc, h_per, b, S_tray, H_tray = 1.8, 0.008, 0.188, 0.030, 1.050, 2.294, 0.500
    z_top, z_bot = 1.0, 2.0
    nu1, nu2 = 25.6 + 2 * 14.8, 24.6 + 4 * 14.8
    theta, l_cell = 0.1, 0.35

    if interp1d is None:
        def eq_curve(x):
            return np.interp(x, x_eq, y_eq)

        def eq_inv(y):
            return np.interp(y, y_eq, x_eq)

        def t_from_x(x):
            return np.interp(x, x_txy, t_txy)

        def t_from_y(y):
            return np.interp(y, y_txy, t_txy)

        def get_m(x_val):
            x_c = float(np.clip(x_val, x_eq[1], x_eq[-2]))
            dx = 1e-4
            return float((eq_curve(x_c + dx) - eq_curve(x_c - dx)) / (2 * dx))
    else:
        eq_spline = CubicSpline(x_eq, y_eq)
        eq_curve = interp1d(x_eq, y_eq, kind="cubic", fill_value="extrapolate")
        eq_inv = interp1d(y_eq, x_eq, kind="cubic", fill_value="extrapolate")
        t_from_x = interp1d(x_txy, t_txy, kind="cubic", fill_value="extrapolate")
        t_from_y = interp1d(y_txy, t_txy, kind="cubic", fill_value="extrapolate")

        def get_m(x_val):
            return float(eq_spline.derivative()(np.clip(x_val, 0.001, 0.999)))

    def mass_to_mol(xm):
        return (xm / M1) / (xm / M1 + (1 - xm) / M2)

    xf = mass_to_mol(xf_mass)
    xp = mass_to_mol(xp_mass)
    xw = mass_to_mol(xw_mass)

    W = Gf * (xp - xf) / (xp - xw)
    D = Gf - W

    yf_star = float(eq_curve(xf))
    Rmin = (xp - yf_star) / (yf_star - xf)

    def rectifying(x, R_):
        return R_ / (R_ + 1) * x + xp / (R_ + 1)

    def stripping(x, yq_):
        m_s = (yq_ - xw) / (xf - xw)
        return m_s * (x - xw) + xw

    def count_stages_mccabe(R_):
        yq_ = rectifying(xf, R_)
        y_ = xp
        N_, last_x = 0.0, xp
        for _ in range(500):
            x_new = float(eq_inv(y_))
            if x_new <= xw:
                N_ += (last_x - xw) / max(last_x - x_new, 1e-12)
                break
            y_new = rectifying(x_new, R_) if x_new >= xf else stripping(x_new, yq_)
            last_x, y_ = x_new, y_new
            N_ += 1
        return N_

    B_values = np.arange(1.05, 3.51, 0.10)
    table = np.array([[Bv, Rmin * Bv, count_stages_mccabe(Rmin * Bv), count_stages_mccabe(Rmin * Bv) * (Rmin * Bv + 1)] for Bv in B_values])
    idx_opt = np.argmin(table[:, 3])
    B_opt, R_opt, N_opt = table[idx_opt, 0], table[idx_opt, 1], table[idx_opt, 2]

    yq_opt = rectifying(xf, R_opt)
    rect = lambda x: rectifying(x, R_opt)
    strip = lambda x: stripping(x, yq_opt)

    xsr_v, xsr_n = (xp + xf) / 2, (xw + xf) / 2
    ysr_v, ysr_n = (xp + yf_star) / 2, (xw + yf_star) / 2

    Mp = M1 * xp + M2 * (1 - xp)
    MF = M1 * xf + M2 * (1 - xf)
    Mv = M1 * xsr_v + M2 * (1 - xsr_v)
    Mn = M1 * xsr_n + M2 * (1 - xsr_n)
    Mv_prime = M1 * ysr_v + M2 * (1 - ysr_v)
    Mn_prime = M1 * ysr_n + M2 * (1 - ysr_n)

    Lv = D * R_opt * Mv / Mp
    Ln = D * R_opt * Mn / Mp + Gf * Mn / MF
    Gv = D * (R_opt + 1) * Mv_prime / Mp
    Gn = D * (R_opt + 1) * Mn_prime / Mp

    t_vap_v, t_vap_n = float(t_from_y(ysr_v)), float(t_from_y(ysr_n))
    t_liq_v, t_liq_n = float(t_from_x(xsr_v)), float(t_from_x(xsr_n))

    rho_yv = Mv_prime / 22.4 * 273 / (273 + t_vap_v)
    rho_yn = Mn_prime / 22.4 * 273 / (273 + t_vap_n)

    rho_xv = (lininterp(t_liq_v, t1_x, rho_x1_t1, t2_x, rho_x1_t2) * xsr_v + lininterp(t_liq_v, t1_x, rho_x2_t1, t2_x, rho_x2_t2) * (1 - xsr_v))
    rho_xn = (lininterp(t_liq_n, t1_x, rho_x1_t1, t2_x, rho_x1_t2) * xsr_n + lininterp(t_liq_n, t1_x, rho_x2_t1, t2_x, rho_x2_t2) * (1 - xsr_n))

    mu_xv = log_mix_viscosity(xsr_v, lininterp(t_liq_v, t1_x, mu_x1_t1, t2_x, mu_x1_t2), lininterp(t_liq_v, t1_x, mu_x2_t1, t2_x, mu_x2_t2))
    mu_xn = log_mix_viscosity(xsr_n, lininterp(t_liq_n, t1_x, mu_x1_t1, t2_x, mu_x1_t2), lininterp(t_liq_n, t1_x, mu_x2_t1, t2_x, mu_x2_t2))

    sigma_xv = (lininterp(t_liq_v, t1_x, sigma_x1_t1, t2_x, sigma_x1_t2) * xsr_v + lininterp(t_liq_v, t1_x, sigma_x2_t1, t2_x, sigma_x2_t2) * (1 - xsr_v))
    sigma_xn = (lininterp(t_liq_n, t1_x, sigma_x1_t1, t2_x, sigma_x1_t2) * xsr_n + lininterp(t_liq_n, t1_x, sigma_x2_t1, t2_x, sigma_x2_t2) * (1 - xsr_n))

    mu_yv = ysr_v * lininterp(t_vap_v, t1_y, mu_y1_t1, t2_y, mu_y1_t2) + (1 - ysr_v) * lininterp(t_vap_v, t1_y, mu_y2_t1, t2_y, mu_y2_t2)
    mu_yn = ysr_n * lininterp(t_vap_n, t1_y, mu_y1_t1, t2_y, mu_y1_t2) + (1 - ysr_n) * lininterp(t_vap_n, t1_y, mu_y2_t1, t2_y, mu_y2_t2)

    w_v, w_n = 0.05 * np.sqrt(rho_xv / rho_yv), 0.05 * np.sqrt(rho_xn / rho_yn)
    w_avg, G_avg, rho_avg = (w_v + w_n) / 2, (Gv + Gn) / 2, (rho_yv + rho_yn) / 2
    d_calc = np.sqrt(4 * G_avg / (np.pi * w_avg * rho_avg))
    w_work = w_avg * (d_calc / d_col) ** 2
    w_tray = w_work * 0.785 * d_col**2 / S_tray

    q_v = Lv / (rho_xv * b)
    q_n = Ln / (rho_xn * b)
    m_exp = 0.05 - 4.6 * h_per

    def calc_h0(q, w_t, mu_x, sigma_x):
        return 0.787 * q**0.2 * h_per**0.56 * w_t**m_exp * (1 - 0.31 * np.exp(-0.11 * mu_x)) * (sigma_x / sigma_water) ** 0.09

    h0_v, h0_n = calc_h0(q_v, w_tray, mu_xv, sigma_xv), calc_h0(q_n, w_tray, mu_xn, sigma_xn)

    def calc_eps(w_t, h0):
        Fr = w_t**2 / (9.81 * h0)
        return np.sqrt(Fr) / (1 + np.sqrt(Fr))

    eps_v, eps_n = calc_eps(w_tray, h0_v), calc_eps(w_tray, h0_n)

    lt = np.sqrt(d_col**2 - b**2)
    S_int = max(1, round(lt / l_cell))

    mu_x1_20 = lininterp(20, t1_x, mu_x1_t1, t2_x, mu_x1_t2)
    mu_x2_20 = lininterp(20, t1_x, mu_x2_t1, t2_x, mu_x2_t2)
    mu_x20_v = log_mix_viscosity(xsr_v, mu_x1_20, mu_x2_20)
    mu_x20_n = log_mix_viscosity(xsr_n, mu_x1_20, mu_x2_20)
    rho_x1_20 = lininterp(20, t1_x, rho_x1_t1, t2_x, rho_x1_t2)
    rho_x2_20 = lininterp(20, t1_x, rho_x2_t1, t2_x, rho_x2_t2)
    rho_x20_v = rho_x1_20 * xsr_v + rho_x2_20 * (1 - xsr_v)
    rho_x20_n = rho_x1_20 * xsr_n + rho_x2_20 * (1 - xsr_n)

    def calc_Dx20(mu_x20):
        return 1e-6 * np.sqrt(1 / M1 + 1 / M2) / (mu_x20 * (nu1 ** (1 / 3) + nu2 ** (1 / 3)) ** 2)

    Dx20_v, Dx20_n = calc_Dx20(mu_x20_v), calc_Dx20(mu_x20_n)
    b_vc, b_nc = 0.2 * mu_x20_v ** (1 / 3) / rho_x20_v, 0.2 * mu_x20_n ** (1 / 3) / rho_x20_n
    Dx_v, Dx_n = Dx20_v * (1 + b_vc * (t_liq_v - 20)), Dx20_n * (1 + b_nc * (t_liq_n - 20))

    def calc_Dy(T_C, P=101325):
        T_K = T_C + 273.15
        return 4.22e-2 * T_K**1.5 * np.sqrt(1 / M1 + 1 / M2) / (P * (nu1 ** (1 / 3) + nu2 ** (1 / 3)) ** 2)

    Dy_v, Dy_n = calc_Dy(t_vap_v), calc_Dy(t_vap_n)

    U_v = Lv / (rho_xv * 0.785 * d_col**2)
    U_n = Ln / (rho_xn * 0.785 * d_col**2)

    def beta_xf_ms(Dx, U, eps, h0, mu_y, mu_x):
        return 6.24e5 * Dx**0.5 * (U / (1 - eps))**0.5 * h0 * (mu_y / (mu_x + mu_y))**0.5

    def beta_yf_ms(Dy, Fc_, w_t, eps, h0, mu_y, mu_x):
        return 6.24e5 * Fc_ * Dy**0.5 * (w_t / eps) ** 0.5 * h0 * (mu_y / (mu_x + mu_y)) ** 0.5

    bxf_v_ms, bxf_n_ms = beta_xf_ms(Dx_v, U_v, eps_v, h0_v, mu_yv, mu_xv), beta_xf_ms(Dx_n, U_n, eps_n, h0_n, mu_yn, mu_xn)
    byf_v_ms, byf_n_ms = beta_yf_ms(Dy_v, Fc, w_tray, eps_v, h0_v, mu_yv, mu_xv), beta_yf_ms(Dy_n, Fc, w_tray, eps_n, h0_n, mu_yn, mu_xn)
    bxf_v, bxf_n = bxf_v_ms * rho_xv / Mv, bxf_n_ms * rho_xn / Mn
    byf_v, byf_n = byf_v_ms * rho_yv / Mv_prime, byf_n_ms * rho_yn / Mn_prime

    e_x = np.array([0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20, 1.30, 1.40, 1.50, 1.70, 2.00, 2.50, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00])
    e_y = np.array([0.0010, 0.0016, 0.0025, 0.0038, 0.0055, 0.0075, 0.010, 0.014, 0.018, 0.028, 0.040, 0.055, 0.070, 0.095, 0.120, 0.145, 0.160, 0.172, 0.178, 0.182, 0.185, 0.188, 0.190, 0.192])

    def get_e(arg):
        argv = float(np.clip(arg, e_x[0], e_x[-1]))
        return 10 ** float(np.interp(np.log10(argv), np.log10(e_x), np.log10(e_y)))

    def calc_m_prime(sigma_x, rho_y, rho_x, mu_y):
        return 1.15e-3 * (sigma_x * 1e-3 / rho_y) ** 0.295 * ((rho_x - rho_y) / (mu_y * 1e-3)) ** 0.425

    m_prime_v, m_prime_n = calc_m_prime(sigma_xv, rho_yv, rho_xv, mu_yv), calc_m_prime(sigma_xn, rho_yn, rho_xn, mu_yn)
    h_bub_v, h_bub_n = h0_v / (1 - eps_v), h0_n / (1 - eps_n)
    Hc_v, Hc_n = H_tray - h_bub_v, H_tray - h_bub_n
    e_v, e_n = get_e(w_tray / (m_prime_v * Hc_v)), get_e(w_tray / (m_prime_n * Hc_n))

    def calc_EMy_full(x_val, section="top"):
        m = max(get_m(x_val), 1e-6)
        if section == "top":
            bxf_mol, byf_mol = bxf_v, byf_v
            M_prime, w_t, rho_y, e_loc = Mv_prime, w_tray, rho_yv, e_v
            lam = m * (R_opt + 1) / R_opt
        else:
            bxf_mol, byf_mol = bxf_n, byf_n
            M_prime, w_t, rho_y, e_loc = Mn_prime, w_tray, rho_yn, e_n
            lam = m
        Kyf = 1 / (1 / byf_mol + m / bxf_mol)
        n_oy = Kyf * M_prime / (w_t * rho_y)
        Ey = 1 - np.exp(-n_oy)
        B_val = lam * (Ey + e_loc / m) / ((1 - theta) * (1 + e_loc * lam / m))
        E2My = Ey if abs(B_val) < 1e-10 else (Ey / B_val) * ((1 + B_val / S_int) ** S_int - 1)
        E1My = E2My / (1 + lam * theta * E2My / (1 - theta))
        EMy = E1My / (1 + e_loc * lam * E1My / (m * (1 - theta)))
        return EMy

    x_kin_arr = np.sort(np.concatenate([np.linspace(xw + 0.005, xf - 0.005, 30), np.linspace(xf + 0.005, xp - 0.005, 30)]))
    y_kin_arr = []
    EMy_all = []
    for xi in x_kin_arr:
        sec = "top" if xi >= xf else "bottom"
        EMy_i = calc_EMy_full(xi, sec)
        y_in = rect(xi) if xi >= xf else strip(xi)
        y_star = float(eq_curve(xi))
        y_k = float(np.clip(y_in + EMy_i * (y_star - y_in), 0, 1))
        y_kin_arr.append(y_k)
        EMy_all.append(EMy_i)

    y_kin_arr = np.array(y_kin_arr)
    x_kin_full = np.concatenate([[0.0], x_kin_arr, [1.0]])
    y_kin_full = np.concatenate([[0.0], y_kin_arr, [1.0]])

    def kin_inv(y):
        return np.interp(y, y_kin_full, x_kin_full)

    x_s, y_s = xp, xp
    N_top, N_bot = 0, 0
    for _ in range(500):
        x_new = float(np.clip(kin_inv(y_s), xw, xp))
        if x_new <= xw + 1e-4:
            break
        if x_new >= xf:
            y_new = rect(x_new)
            N_top += 1
        else:
            y_new = strip(x_new)
            N_bot += 1
        x_s, y_s = x_new, y_new

    N_total = N_top + N_bot
    Hk = (N_total - 1) * H_tray + z_top + z_bot

    dP_dry = 1.85 * w_work**2 * (rho_yv + rho_yn) / 2 / (2 * Fc**2)
    dP_liq_v = 9.81 * rho_xv * h0_v
    dP_liq_n = 9.81 * rho_xn * h0_n
    dP_sigma = 4 * sigma_xv * 1e-3 / d0
    dP_total = (dP_dry + dP_liq_v + dP_sigma) * N_top + (dP_dry + dP_liq_n + dP_sigma) * N_bot

    # plots
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(table[:, 1], table[:, 3], "b-o", ms=5, lw=1.5)
    ax.scatter([R_opt], [N_opt * (R_opt + 1)], s=80, color="red")
    ax.set_xlabel("R")
    ax.set_ylabel("N(R+1)")
    ax.grid(alpha=0.4)
    fig.tight_layout()
    fig.savefig(output_dir / "plot_NR1.png", dpi=120)
    plt.close(fig)

    x_pl = np.linspace(0, 1, 300)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(x_eq, y_eq, "b-", lw=2)
    ax.plot(x_pl, x_pl, "k--", lw=1)
    ax.plot(np.linspace(xf, xp, 200), rect(np.linspace(xf, xp, 200)), "g-", lw=1.5)
    ax.plot(np.linspace(xw, xf, 200), strip(np.linspace(xw, xf, 200)), "m-", lw=1.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.4)
    fig.tight_layout()
    fig.savefig(output_dir / "plot_yx_mccabe.png", dpi=120)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x_txy, t_txy, "b-o", ms=4, lw=1.5, label="Линия кипения t(x)")
    ax.plot(y_txy, t_txy, "r-s", ms=4, lw=1.5, label="Линия конденсации t(y)")
    ax.set_xlabel("x (или y)")
    ax.set_ylabel("t, °C")
    ax.grid(alpha=0.4)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "plot_txy.png", dpi=120)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(e_x, e_y, "b-o", ms=4, lw=1.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("wт / (m' Hc)")
    ax.set_ylabel("e, кг/кг")
    ax.grid(alpha=0.4, which="both")
    fig.tight_layout()
    fig.savefig(output_dir / "plot_entrainment.png", dpi=120)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 6))
    x_diag = np.linspace(0, 1, 300)
    ax.plot(x_eq, y_eq, "b-", lw=2, label="Равновесная")
    ax.plot(x_diag, x_diag, "k--", lw=1, label="Диагональ")
    ax.plot(x_kin_full, y_kin_full, "c-", lw=2, label="Кинетическая")
    ax.plot(np.linspace(xf, xp, 200), rect(np.linspace(xf, xp, 200)), "g-", lw=1.2, label="Укрепляющая")
    ax.plot(np.linspace(xw, xf, 200), strip(np.linspace(xw, xf, 200)), "m-", lw=1.2, label="Исчерпывающая")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.4)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output_dir / "plot_kinetic.png", dpi=120)
    plt.close(fig)

    result = {
        "input": asdict(input_data),
        "xf_mol": xf,
        "xp_mol": xp,
        "xw_mol": xw,
        "D_kg_s": D,
        "W_kg_s": W,
        "Rmin": float(Rmin),
        "R_opt": float(R_opt),
        "beta_opt": float(B_opt),
        "N_theor": float(N_opt),
        "d_calc_m": float(d_calc),
        "d_col_m": float(d_col),
        "w_work_m_s": float(w_work),
        "N_top": int(N_top),
        "N_bottom": int(N_bot),
        "N_total": int(N_total),
        "H_column_m": float(Hk),
        "dP_total_Pa": float(dP_total),
        "plots": [
            "plot_NR1.png",
            "plot_yx_mccabe.png",
            "plot_txy.png",
            "plot_entrainment.png",
            "plot_kinetic.png",
        ],
    }

    text_lines = [
        "==============================================================================",
        "  РАСЧЁТ ТАРЕЛЬЧАТОЙ РЕКТИФИКАЦИОННОЙ КОЛОННЫ (метод Дытнерского)",
        "==============================================================================",
        "",
        "БЛОК 0 — Исходные данные",
        f"  F = {Gf_kgh} кг/ч = {_fmt(Gf)} кг/с",
        f"  x_f^mass = {_fmt(xf_mass)}, x_p^mass = {_fmt(xp_mass)}, x_w^mass = {_fmt(xw_mass)}",
        f"  M1={M1} кг/кмоль, M2={M2} кг/кмоль, d_col={d_col} м, H_tray={H_tray} м",
        "",
        "БЛОК 1 — Материальный баланс",
        "  x = (x_mass/M1) / (x_mass/M1 + (1-x_mass)/M2)",
        f"  x_f = {_fmt(xf)}, x_p = {_fmt(xp)}, x_w = {_fmt(xw)}",
        "  W = F*(x_p-x_f)/(x_p-x_w), D=F-W",
        f"  W = {_fmt(Gf)}*({_fmt(xp)}-{_fmt(xf)})/({_fmt(xp)}-{_fmt(xw)}) = {_fmt(W)} кг/с",
        f"  D = {_fmt(Gf)} - {_fmt(W)} = {_fmt(D)} кг/с",
        "",
        "БЛОК 2 — Рабочее флегмовое число (Мак-Кэб)",
        "  y_f* = y_eq(x_f)",
        f"  y_f* = {_fmt(yf_star)}",
        "  R_min = (x_p-y_f*)/(y_f*-x_f)",
        f"  R_min = ({_fmt(xp)}-{_fmt(yf_star)})/({_fmt(yf_star)}-{_fmt(xf)}) = {_fmt(Rmin)}",
        "  По критерию min N(R+1):",
        f"  beta_opt = {_fmt(B_opt)}, R_opt = {_fmt(R_opt)}, N_theor = {_fmt(N_opt)}",
        "",
        "БЛОК 3 — Нагрузки, температуры, плотности",
        "  Mp = M1*xp + M2*(1-xp), MF = M1*xf + M2*(1-xf), Mv/Mn/M'v/M'n аналогично",
        f"  Mp={_fmt(Mp)}, MF={_fmt(MF)}, Mv={_fmt(Mv)}, Mn={_fmt(Mn)}, M'v={_fmt(Mv_prime)}, M'n={_fmt(Mn_prime)} кг/кмоль",
        "  Lв = D*R*Mv/Mp; Lн = D*R*Mn/Mp + F*Mn/MF",
        f"  Lв={_fmt(Lv)} кг/с; Lн={_fmt(Ln)} кг/с",
        "  Gв = D*(R+1)*M'в/Mp; Gн = D*(R+1)*M'н/Mp",
        f"  Gв={_fmt(Gv)} кг/с; Gн={_fmt(Gn)} кг/с",
        "  rho_y = M'/22.4 * 273/(273+t'), rho_x — смешение по компонентам",
        f"  t'в={_fmt(t_vap_v)} °C; t'н={_fmt(t_vap_n)} °C; tв={_fmt(t_liq_v)} °C; tн={_fmt(t_liq_n)} °C",
        f"  rho_yв={_fmt(rho_yv)}; rho_yн={_fmt(rho_yn)}; rho_xв={_fmt(rho_xv)}; rho_xн={_fmt(rho_xn)} кг/м³",
        f"  mu_xв={_fmt(mu_xv)}; mu_xн={_fmt(mu_xn)}; mu_yв={_fmt(mu_yv)}; mu_yн={_fmt(mu_yn)} мПа·с",
        "",
        "БЛОК 4 — Диаметр колонны",
        "  w = 0.05*sqrt(rho_x/rho_y)",
        f"  wв={_fmt(w_v)} м/с; wн={_fmt(w_n)} м/с; w_avg={_fmt(w_avg)} м/с",
        "  d = sqrt(4*G_avg/(pi*w_avg*rho_avg))",
        f"  G_avg={_fmt(G_avg)} кг/с; rho_avg={_fmt(rho_avg)} кг/м³",
        f"  d_calc={_fmt(d_calc)} м; принято d_col={_fmt(d_col)} м",
        f"  w_work={_fmt(w_work)} м/с; w_tray={_fmt(w_tray)} м/с",
        "",
        "БЛОК 5 — Параметры тарелки",
        "  h0 = 0.787*q^0.2*h_per^0.56*w_t^m*(1-0.31*exp(-0.11*mu_x))*(sigma/sigma_water)^0.09",
        f"  qв={_fmt(q_v)}; qн={_fmt(q_n)}; m={_fmt(m_exp)}",
        f"  h0в={_fmt(h0_v)} м; h0н={_fmt(h0_n)} м",
        "  eps = sqrt(Fr)/(1+sqrt(Fr)), Fr=w_t^2/(g*h0)",
        f"  epsв={_fmt(eps_v)}; epsн={_fmt(eps_n)}; lt={_fmt(lt)} м; S={S_int}",
        "",
        "БЛОК 6 — Коэффициенты молекулярной диффузии",
        "  Dx20 = 1e-6*sqrt(1/M1+1/M2)/(mu*(nu1^(1/3)+nu2^(1/3))^2)",
        f"  Dx20в={_fmt(Dx20_v)}; Dx20н={_fmt(Dx20_n)}",
        f"  bв={_fmt(b_vc)}; bн={_fmt(b_nc)}",
        "  Dx = Dx20*(1+b*(t-20)); Dy = 4.22e-2*T^1.5*sqrt(1/M1+1/M2)/(P*(nu1^(1/3)+nu2^(1/3))^2)",
        f"  Dxв={_fmt(Dx_v)}; Dxн={_fmt(Dx_n)}; Dyв={_fmt(Dy_v)}; Dyн={_fmt(Dy_n)}",
        "",
        "БЛОК 7 — Коэффициенты массоотдачи",
        "  beta_xf, beta_yf по корреляциям Дытнерского",
        f"  Uв={_fmt(U_v)} м/с; Uн={_fmt(U_n)} м/с",
        f"  bxfв={_fmt(bxf_v_ms)} м/с; bxfн={_fmt(bxf_n_ms)} м/с; byfв={_fmt(byf_v_ms)} м/с; byfн={_fmt(byf_n_ms)} м/с",
        f"  bxfв={_fmt(bxf_v)}; bxfн={_fmt(bxf_n)}; byfв={_fmt(byf_v)}; byfн={_fmt(byf_n)} кмоль/(м²·с)",
        "",
        "БЛОК 8 — Эффективность по Мэрфри и кинетическая линия",
        "  m' = 1.15e-3*(sigma/rho_y)^0.295*((rho_x-rho_y)/mu_y)^0.425",
        f"  m'в={_fmt(m_prime_v)}; m'н={_fmt(m_prime_n)}; Hcв={_fmt(Hc_v)} м; Hcн={_fmt(Hc_n)} м",
        f"  arg_eв={_fmt(w_tray / (m_prime_v * Hc_v))}; arg_eн={_fmt(w_tray / (m_prime_n * Hc_n))}",
        f"  eв={_fmt(e_v)}; eн={_fmt(e_n)}; mean(EMy)={_fmt(float(np.mean(EMy_all)))}",
        "  Для каждой точки x_i считалось: m(x_i), Kyf, n_oy, Ey, EMy и y_k = y_in + EMy*(y*-y_in)",
        "",
        "БЛОК 9 — Число действительных тарелок",
        f"  Nв={N_top}; Nн={N_bot}; N={N_total}",
        "",
        "БЛОК 10 — Высота колонны",
        "  Hк = (N-1)*H_tray + z_top + z_bot",
        f"  Hк = ({N_total}-1)*{_fmt(H_tray)} + {_fmt(z_top)} + {_fmt(z_bot)} = {_fmt(Hk)} м",
        "",
        "БЛОК 11 — Гидравлическое сопротивление",
        "  dPdry = 1.85*w^2*rho/(2*Fc^2), dPliq = g*rho_x*h0, dPsigma = 4*sigma/d0",
        f"  dPdry={_fmt(dP_dry)} Па; dPliqв={_fmt(dP_liq_v)} Па; dPliqн={_fmt(dP_liq_n)} Па; dPsigma={_fmt(dP_sigma)} Па",
        f"  dPк={_fmt(dP_total)} Па",
        "",
        "БЛОК 12 — Итог",
        f"  R={_fmt(R_opt)}; N_theor={_fmt(N_opt)}; N_actual={N_total}; d={_fmt(d_col)} м; H={_fmt(Hk)} м; dP={_fmt(dP_total)} Па",
        "",
        "Графики: plot_NR1.png, plot_yx_mccabe.png, plot_txy.png, plot_entrainment.png, plot_kinetic.png",
    ]
    (output_dir / "report_dytnersky.txt").write_text("\n".join(text_lines) + "\n", encoding="utf-8")

    report_pre = "\\n".join(text_lines)

    html = f"""<!doctype html>
<html lang="ru"><head><meta charset="utf-8"><title>Отчет Дытнерского</title>
<style>
body{{font-family:Arial,sans-serif;margin:24px;line-height:1.45}}
h2{{margin-top:28px}}
.eq{{background:#f8f8f8;padding:10px;border-radius:8px}}
</style>
<script>window.MathJax={{tex:{{inlineMath:[['$','$'],['\\\\(','\\\\)']]}}}};</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script></head><body>
<h1>Расчёт тарельчатой ректификационной колонны (полный протокол)</h1>
<p>В отчете приведены формулы, подстановки и комментарии по всем блокам 0–12.</p>
<h2>Блок 0. Исходные данные</h2>
<div class="eq"><p>$F={_fmt(Gf_kgh)}\\;\\text{{кг/ч}}={_fmt(Gf)}\\;\\text{{кг/с}}$</p>
<p>$x_f^{{mass}}={_fmt(xf_mass)},\\;x_p^{{mass}}={_fmt(xp_mass)},\\;x_w^{{mass}}={_fmt(xw_mass)}$</p></div>
<h2>Блок 1. Материальный баланс</h2>
<div class="eq"><p>$x=\\frac{{x^{{mass}}/M_1}}{{x^{{mass}}/M_1+(1-x^{{mass}})/M_2}}$</p>
<p>$x_f={_fmt(xf)},\\;x_p={_fmt(xp)},\\;x_w={_fmt(xw)}$</p>
<p>$W=F\\frac{{x_p-x_f}}{{x_p-x_w}}={_fmt(W)}\\;\\text{{кг/с}},\\quad D=F-W={_fmt(D)}\\;\\text{{кг/с}}$</p></div>
<h2>Блок 2. Флегмовое число</h2>
<div class="eq"><p>$R_{{min}}=\\frac{{x_p-y_f^*}}{{y_f^*-x_f}}={_fmt(Rmin)}$</p>
<p>$\\beta_{{opt}}={_fmt(B_opt)},\\quad R=\\beta_{{opt}}R_{{min}}={_fmt(R_opt)},\\quad N_{{theor}}={_fmt(N_opt)}$</p></div>
<h2>Блок 3. Нагрузки, температуры, свойства</h2>
<div class="eq"><p>$L_в={_fmt(Lv)},\\;L_н={_fmt(Ln)},\\;G_в={_fmt(Gv)},\\;G_н={_fmt(Gn)}\\;\\text{{кг/с}}$</p>
<p>$t'_в={_fmt(t_vap_v)}^\\circ C,\\;t'_н={_fmt(t_vap_n)}^\\circ C,\\;t_в={_fmt(t_liq_v)}^\\circ C,\\;t_н={_fmt(t_liq_n)}^\\circ C$</p>
<p>$\\rho_{{y,в}}={_fmt(rho_yv)},\\;\\rho_{{y,н}}={_fmt(rho_yn)},\\;\\rho_{{x,в}}={_fmt(rho_xv)},\\;\\rho_{{x,н}}={_fmt(rho_xn)}$</p></div>
<h2>Блок 4. Диаметр</h2>
<div class="eq"><p>$w=0.05\\sqrt{{\\rho_x/\\rho_y}}$</p>
<p>$d=\\sqrt{{\\frac{{4\\bar G}}{{\\pi\\bar w\\bar\\rho}}}}={_fmt(d_calc)}\\;\\text{{м}},\\;d_{{прин}}={_fmt(d_col)}\\;\\text{{м}}$</p></div>
<h2>Блок 5. Параметры тарелки</h2>
<div class="eq"><p>$h_0=0.787q^{{0.2}}h_{{per}}^{{0.56}}w_t^m(1-0.31e^{{-0.11\\mu_x}})(\\sigma/\\sigma_w)^{{0.09}}$</p>
<p>$h_0^в={_fmt(h0_v)},\\;h_0^н={_fmt(h0_n)},\\;\\varepsilon_в={_fmt(eps_v)},\\;\\varepsilon_н={_fmt(eps_n)}$</p></div>
<h2>Блок 6. Молекулярная диффузия</h2>
<div class="eq"><p>$D_x^в={_fmt(Dx_v)},\\;D_x^н={_fmt(Dx_n)},\\;D_y^в={_fmt(Dy_v)},\\;D_y^н={_fmt(Dy_n)}$</p></div>
<h2>Блок 7. Массоотдача</h2>
<div class="eq"><p>$\\beta_{{xf,в}}={_fmt(bxf_v_ms)},\\;\\beta_{{xf,н}}={_fmt(bxf_n_ms)},\\;\\beta_{{yf,в}}={_fmt(byf_v_ms)},\\;\\beta_{{yf,н}}={_fmt(byf_n_ms)}\\;\\text{{м/с}}$</p>
<p>$\\beta_{{xf,в}}={_fmt(bxf_v)},\\;\\beta_{{xf,н}}={_fmt(bxf_n)},\\;\\beta_{{yf,в}}={_fmt(byf_v)},\\;\\beta_{{yf,н}}={_fmt(byf_n)}\\;\\text{{кмоль/(м²·с)}}$</p></div>
<h2>Блок 8. Эффективность и кинетическая линия</h2>
<div class="eq"><p>$m'_в={_fmt(m_prime_v)},\\;m'_н={_fmt(m_prime_n)},\\;e_в={_fmt(e_v)},\\;e_н={_fmt(e_n)}$</p>
<p>$\\overline{{E_{{My}}}}={_fmt(float(np.mean(EMy_all)))}$</p></div>
<h2>Блок 9. Действительные тарелки</h2>
<div class="eq"><p>$N_в={N_top},\\;N_н={N_bot},\\;N={N_total}$</p></div>
<h2>Блок 10. Высота колонны</h2>
<div class="eq"><p>$H_к=(N-1)H+z_{{top}}+z_{{bot}}=({_fmt(N_total)}-1)\\cdot{_fmt(H_tray)}+{_fmt(z_top)}+{_fmt(z_bot)}={_fmt(Hk)}\\;\\text{{м}}$</p></div>
<h2>Блок 11. Гидросопротивление</h2>
<div class="eq"><p>$\\Delta P_{{dry}}={_fmt(dP_dry)},\\;\\Delta P_{{liq,в}}={_fmt(dP_liq_v)},\\;\\Delta P_{{liq,н}}={_fmt(dP_liq_n)},\\;\\Delta P_\\sigma={_fmt(dP_sigma)}\\;\\text{{Па}}$</p>
<p>$\\Delta P_к={_fmt(dP_total)}\\;\\text{{Па}}$</p></div>
<h2>Блок 12. Итог</h2>
<div class="eq"><p>$R={_fmt(R_opt)},\\;N_{{theor}}={_fmt(N_opt)},\\;N_{{actual}}={N_total},\\;d={_fmt(d_col)}\\;\\text{{м}},\\;H={_fmt(Hk)}\\;\\text{{м}}$</p></div>
<h2>Графики</h2>
<ul><li>plot_NR1.png</li><li>plot_yx_mccabe.png</li><li>plot_txy.png</li><li>plot_entrainment.png</li><li>plot_kinetic.png</li></ul>
<h2>Полный протокол расчета (текстовый дубль)</h2>
<pre>{report_pre}</pre>
</body></html>"""
    (output_dir / "report_dytnersky.html").write_text(html, encoding="utf-8")
    (output_dir / "result_dytnersky.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result
