from __future__ import annotations

import argparse
import json
from pathlib import Path

from src.calculations import InputData, calculate
from src.dytnersky_method import DytnerskyInput, run_dytnersky
from src.reporting import save_reports


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Расчет ректификационной колонны + отчет")
    parser.add_argument("--input", required=True, help="Путь к JSON с исходными данными")
    parser.add_argument("--output-dir", default="out", help="Директория для отчетов")
    parser.add_argument(
        "--method",
        choices=["prototype", "dytnersky"],
        default="dytnersky",
        help="Метод расчета: упрощенный prototype или расширенный dytnersky",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)

    data = json.loads(input_path.read_text(encoding="utf-8"))

    if args.method == "prototype":
        model = InputData(**data)
        result = calculate(model)
        paths = save_reports(result, output_dir=output_dir, template_dir=Path("templates"))

        print("Расчет завершен (prototype).")
        print(f"TXT:  {paths['txt']}")
        print(f"HTML: {paths['html']}")
        print(f"PDF:  {paths['pdf'] if paths['pdf'] else 'не создан (weasyprint недоступен)'}")
    else:
        model = DytnerskyInput(**data)
        result = run_dytnersky(model, output_dir=output_dir)
        print("Расчет завершен (dytnersky).")
        print(f"TXT:  {output_dir / 'report_dytnersky.txt'}")
        print(f"HTML: {output_dir / 'report_dytnersky.html'}")
        print(f"JSON: {output_dir / 'result_dytnersky.json'}")
        print(f"N_total = {result['N_total']}, H = {result['H_column_m']:.3f} м")


if __name__ == "__main__":
    main()
