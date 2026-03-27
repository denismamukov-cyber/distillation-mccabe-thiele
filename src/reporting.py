from __future__ import annotations

from pathlib import Path
from typing import Optional

from src.calculations import CalculationResult


def fmt(value: float, digits: int = 6) -> str:
    return f"{value:.{digits}g}"


def build_text_report(result: CalculationResult) -> str:
    s = []
    s.append("РАСЧЕТ РЕКТИФИКАЦИОННОЙ КОЛОННЫ (ПРОТОТИП)\n")

    s.append("1. Исходные данные")
    for k, v in result.inputs.items():
        s.append(f"   - {k} = {fmt(v)}")

    x_d = result.inputs["distillate_composition_light"]
    x_b = result.inputs["bottoms_composition_light"]
    z_f = result.inputs["feed_composition_light"]
    f = result.inputs["feed_flow_kmol_h"]

    s.append("\n2. Материальный баланс")
    s.append("   Формулы:")
    s.append("   F = D + B")
    s.append("   F*zF = D*xD + B*xB")
    s.append("   D = F*(zF - xB)/(xD - xB)")
    s.append("   B = F - D")
    s.append(
        f"   Подстановка: D = {fmt(f)}*({fmt(z_f)} - {fmt(x_b)})/({fmt(x_d)} - {fmt(x_b)}) = {fmt(result.distillate_flow_kmol_h)} кмоль/ч"
    )
    s.append(f"   B = {fmt(f)} - {fmt(result.distillate_flow_kmol_h)} = {fmt(result.bottoms_flow_kmol_h)} кмоль/ч")

    s.append("\n3. Флегмовое число")
    s.append(f"   R_min = {fmt(result.r_min)}")
    s.append(
        f"   R_work = factor*R_min = {fmt(result.inputs['reflux_ratio_factor_to_rmin'])}*{fmt(result.r_min)} = {fmt(result.r_work)}"
    )

    s.append("\n4. Минимальное число тарелок (Фенске)")
    s.append("   N_min = ln((xD/(1-xD))*((1-xB)/xB))/ln(alpha)")
    s.append(f"   N_min = {fmt(result.n_min_fenske)}")

    s.append("\n5. Число теоретических тарелок (МакКэб–Тили)")
    s.append(f"   N_theoretical = {result.n_theoretical}")
    s.append("   Профиль по ступеням (stage: x, y):")
    for p in result.stage_profile:
        s.append(f"   {p.stage:>3}: x={fmt(p.x)}, y={fmt(p.y)}")

    s.append("\n6. Гидравлика и диаметр")
    s.append(f"   V = (R_work + 1)*D = {fmt(result.vapor_flow_kmol_h)} кмоль/ч")
    s.append(f"   Qv = {fmt(result.vapor_vol_flow_m3_s)} м³/с")
    s.append(f"   D_col = sqrt(4*Qv/(pi*u_allow)) = {fmt(result.column_diameter_m)} м")

    s.append(
        "\nПримечание: формулы гидравлики и R_min здесь даны как рабочий шаблон. "
        "Для строгого совпадения с Mathcad-проектом перенесите ваши исходные корреляции и коэффициенты."
    )

    return "\n".join(s) + "\n"


def build_html_report(result: CalculationResult) -> str:
    rows = "\n".join(
        f"<tr><td>{p.stage}</td><td>{fmt(p.x)}</td><td>{fmt(p.y)}</td></tr>" for p in result.stage_profile
    )
    input_rows = "\n".join(
        f"<li><code>{k}</code> = {fmt(v)}</li>" for k, v in result.inputs.items()
    )

    return f"""<!doctype html>
<html lang=\"ru\">
<head>
  <meta charset=\"utf-8\" />
  <title>Отчет расчета ректификационной колонны</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 28px; line-height: 1.4; }}
    h1, h2 {{ margin-bottom: 0.25em; }}
    code {{ background: #f2f2f2; padding: 2px 6px; border-radius: 4px; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 8px; }}
    th, td {{ border: 1px solid #ddd; padding: 6px; text-align: right; }}
    th:first-child, td:first-child {{ text-align: center; }}
  </style>
  <script>
    window.MathJax = {{ tex: {{ inlineMath: [['$', '$'], ['\\\\(', '\\\\)']] }} }};
  </script>
  <script src=\"https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js\"></script>
</head>
<body>
  <h1>Расчет ректификационной колонны</h1>
  <p><strong>Режим:</strong> автоматизированный шаблон, совместимый с методикой МакКэба–Тили.</p>

  <h2>1. Исходные данные</h2>
  <ul>{input_rows}</ul>

  <h2>2. Материальный баланс</h2>
  <p>$F = D + B$</p>
  <p>$F z_F = D x_D + B x_B$</p>
  <p>$D = \\frac{{F(z_F - x_B)}}{{x_D - x_B}}, B=F-D$</p>
  <p>Результат: $D = {fmt(result.distillate_flow_kmol_h)}$ кмоль/ч, $B = {fmt(result.bottoms_flow_kmol_h)}$ кмоль/ч.</p>

  <h2>3. Флегмовое число</h2>
  <p>$R_{{min}} = {fmt(result.r_min)}$</p>
  <p>$R = {fmt(result.r_work)}$</p>

  <h2>4. Тарелки</h2>
  <p>$N_{{min}}$ (Фенске): {fmt(result.n_min_fenske)}</p>
  <p>$N_{{theor}}$ (МакКэб–Тили): {result.n_theoretical}</p>

  <table>
    <thead><tr><th>Ступень</th><th>x</th><th>y</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>

  <h2>5. Гидравлический расчет (упрощенный)</h2>
  <p>$V = (R+1)D = {fmt(result.vapor_flow_kmol_h)}$ кмоль/ч</p>
  <p>$Q_v = {fmt(result.vapor_vol_flow_m3_s)}$ м³/с</p>
  <p>$D_{{col}} = \\sqrt{{\\frac{{4Q_v}}{{\\pi u_{{allow}}}}}} = {fmt(result.column_diameter_m)}$ м</p>

  <hr>
  <p><em>Для строгого совпадения с Mathcad Prime 10 перенесите оригинальные корреляции и коэффициенты в модуль расчета.</em></p>
</body></html>
"""


def save_reports(result: CalculationResult, output_dir: Path, template_dir: Path | None = None) -> dict[str, Optional[Path]]:
    _ = template_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    txt_path = output_dir / "report.txt"
    txt_path.write_text(build_text_report(result), encoding="utf-8")

    html_path = output_dir / "report.html"
    html_content = build_html_report(result)
    html_path.write_text(html_content, encoding="utf-8")

    pdf_path: Optional[Path] = None
    try:
        from weasyprint import HTML

        pdf_path = output_dir / "report.pdf"
        HTML(string=html_content).write_pdf(str(pdf_path))
    except Exception:
        pdf_path = None

    return {"txt": txt_path, "html": html_path, "pdf": pdf_path}
