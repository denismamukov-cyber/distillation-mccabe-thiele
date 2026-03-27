from __future__ import annotations

import argparse
import json
from pathlib import Path

from src.calculations import InputData, calculate
from src.reporting import save_reports


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Расчет ректификационной колонны + отчет")
    parser.add_argument("--input", required=True, help="Путь к JSON с исходными данными")
    parser.add_argument("--output-dir", default="out", help="Директория для отчетов")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)

    data = json.loads(input_path.read_text(encoding="utf-8"))
    model = InputData(**data)

    result = calculate(model)
    paths = save_reports(result, output_dir=output_dir, template_dir=Path("templates"))

    print("Расчет завершен.")
    print(f"TXT:  {paths['txt']}")
    print(f"HTML: {paths['html']}")
    print(f"PDF:  {paths['pdf'] if paths['pdf'] else 'не создан (weasyprint недоступен)'}")


if __name__ == "__main__":
    main()
