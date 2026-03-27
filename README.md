# Distillation McCabe–Thiele Calculator

Теперь проект поддерживает **два режима расчета**:

1. `dytnersky` (по умолчанию) — расширенный порядок расчета по структуре блоков Дытнерского (на основе вашей вариации).
2. `prototype` — упрощенный шаблонный расчет (старый режим).

## Запуск локально

### Режим Дытнерского (рекомендуется)

```bash
python cli.py --method dytnersky --input examples/input.json --output-dir out
```

Результаты:
- `out/report_dytnersky.txt`
- `out/report_dytnersky.html`
- `out/result_dytnersky.json`
- графики `plot_NR1.png`, `plot_yx_mccabe.png`, `plot_txy.png`, `plot_entrainment.png`

### Упрощенный режим

```bash
python cli.py --method prototype --input examples/input_prototype.json --output-dir out
```

## Входные данные

### `dytnersky`

```json
{
  "g_feed_kg_h": 12000.0,
  "xf_mass": 0.325,
  "xp_mass": 0.97,
  "xw_mass": 0.012
}
```

### `prototype`

Старый формат параметров сохранен; пример в `examples/input_prototype.json`.

## Запуск на GitHub

В workflow `.github/workflows/ci.yml` добавлен запуск режима `dytnersky`.
Артефакт `distillation-report` содержит файлы из директории `out/`.

## Что важно

- В коде перенесены порядок и ключевые формулы из вашего скрипта (матбаланс, подбор `R`, нагрузки, диаметр, гидравлика, эффективность и расчет действительных тарелок).
- Для интерполяции используется `scipy` при наличии; если `scipy` недоступен — включается fallback на линейную интерполяцию `numpy`.
