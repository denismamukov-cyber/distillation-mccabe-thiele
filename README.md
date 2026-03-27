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

## Как запустить локально

```bash
python -m venv .venv
source .venv/bin/activate
python cli.py --input examples/input.json --output-dir out
```

После выполнения появятся файлы:
- `out/report.txt`
- `out/report.html`
- `out/report.pdf` (только если доступен `weasyprint`)

## Как запустить код на GitHub

### Вариант 1: GitHub Actions (автоматически)

В репозитории добавлен workflow `.github/workflows/ci.yml`.
Он запускается на каждый `push` и `pull_request`, выполняет:

```bash
python cli.py --input examples/input.json --output-dir out
```

И загружает результаты из папки `out/` как артефакт `distillation-report`.

Как посмотреть результат:
1. Откройте вкладку **Actions** в вашем репозитории.
2. Выберите последний запуск workflow **ci**.
3. Скачайте артефакт **distillation-report**.

### Вариант 2: Через Codespaces в браузере

1. На странице репозитория нажмите **Code → Codespaces → Create codespace**.
2. В терминале Codespaces выполните:

```bash
python cli.py --input examples/input.json --output-dir out
```

3. Откройте папку `out` и скачайте отчеты.

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
- `vapor_density_kg_m3` — плотность пара для гидравлики
- `souders_brown_factor_m_s` — допустимая скорость пара

## Примечание по PDF

`report.pdf` создается, если установлен пакет `weasyprint` и его системные зависимости.
Если библиотека недоступна, программа завершится успешно и оставит TXT/HTML.
