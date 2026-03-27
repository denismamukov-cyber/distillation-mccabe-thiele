# Distillation McCabe–Thiele Calculator

Прототип программы для автоматизации технологического расчета ректификационной колонны
(баланс, число теоретических тарелок по МакКэбу–Тили, упрощенный гидравлический расчет, подбор диаметра)
с генерацией подробного отчета в TXT/HTML/PDF.

## Что умеет

- Принимает исходные данные из JSON.
- Выполняет:
  - материальный баланс (F, D, B);
  - оценку `R_min` (упрощенно), выбор рабочего флегмового числа;
  - расчет `N_min` по Фенске;
  - пошаговый расчет тарелок методом МакКэба–Тили (численно);
  - упрощенный расчет диаметра колонны по допустимой скорости пара.
- Формирует отчет:
  - `report.txt` (подробный, с подстановкой чисел);
  - `report.html` (с LaTeX-формулами для MathJax);
  - `report.pdf` (если установлен `weasyprint`).

> ⚠️ Важно: чтобы формулы совпали **строго** с вашим Mathcad Prime 10, замените/уточните формулы в `src/calculations.py` согласно вашему методическому расчету.

## Быстрый старт

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python cli.py --input examples/input.json --output-dir out
```

## Формат входных данных

См. `examples/input.json`.

Ключевые параметры:

- `feed_flow_kmol_h` — расход питания F, кмоль/ч
- `feed_composition_light` — состав легколетучего в питании zF
- `distillate_composition_light` — состав легколетучего в дистилляте xD
- `bottoms_composition_light` — состав легколетучего в кубовом остатке xB
- `relative_volatility` — относительная летучесть α
- `reflux_ratio_factor_to_rmin` — множитель к Rmin (например, 1.3…1.8)
- `feed_thermal_condition_q` — тепловое состояние питания q
- `molar_mass_vapor_kg_kmol`, `vapor_density_kg_m3` — параметры для гидравлики
- `souders_brown_factor_m_s` — допустимая скорость пара

## Примечание по PDF

`report.pdf` создается, если установлен пакет `weasyprint` и его системные зависимости.
Если библиотека недоступна, программа завершится успешно и оставит TXT/HTML.
