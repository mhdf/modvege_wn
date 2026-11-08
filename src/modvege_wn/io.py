"""CSV input and output using only the Python standard library."""

import csv
from dataclasses import asdict
from datetime import date
from pathlib import Path

from .model import SimulationResult, WeatherDay


def read_weather_csv(path: str | Path) -> list[WeatherDay]:
    """Read ISO dates and daily weather with the column names in README.md."""
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        needed = {"date", "max_temp_c", "min_temp_c", "reference_et_mm",
                  "solar_radiation_mj_m2", "precipitation_mm"}
        if not needed.issubset(reader.fieldnames or []):
            raise ValueError(f"weather CSV needs columns: {', '.join(sorted(needed))}")
        return [WeatherDay(
            date=date.fromisoformat(row["date"]),
            max_temp_c=float(row["max_temp_c"]), min_temp_c=float(row["min_temp_c"]),
            reference_et_mm=float(row["reference_et_mm"]),
            solar_radiation_mj_m2=float(row["solar_radiation_mj_m2"]),
            precipitation_mm=float(row["precipitation_mm"]),
        ) for row in reader]


def write_results(result: SimulationResult, output_dir: str | Path) -> None:
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    for name, rows, total_field in (
        ("daily.csv", result.daily, "standing_biomass_kg_ha"),
        ("harvests.csv", result.harvests, "yield_kg_ha"),
    ):
        if not rows:
            continue
        records = [{**asdict(row), total_field: getattr(row, total_field)} for row in rows]
        with (directory/name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=records[0].keys())
            writer.writeheader()
            writer.writerows(records)
