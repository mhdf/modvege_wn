"""Minimal command line run for a single site and fixed nitrogen index."""

import argparse
from datetime import date

from .io import read_weather_csv, write_results
from .model import simulate
from .parameters import ModelConfig, Soil


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the ModVege-WN daily pasture model")
    parser.add_argument("weather_csv", help="daily weather CSV; see README.md")
    parser.add_argument("--cut", action="append", default=[], help="cut date YYYY-MM-DD; repeat as needed")
    parser.add_argument("--soil", nargs=4, type=float, metavar=("CLAY_PCT", "SAND_PCT", "OM_G_KG", "ROOT_M"), required=True)
    parser.add_argument("--groups", nargs=4, type=float, metavar=("A", "B", "C", "D"), required=True)
    parser.add_argument("--ni", type=float, default=.85, help="fixed ModVege nutrition index")
    parser.add_argument("--groundwater-floor", type=float, default=0.0)
    parser.add_argument("--water-response", choices=("fao56", "calanca"), default="fao56")
    parser.add_argument("--environment-response", choices=("minimum", "product"), default="minimum")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    soil = Soil(*args.soil)
    config = ModelConfig(nitrogen_index=args.ni, groundwater_floor=args.groundwater_floor,
                         water_response=args.water_response,
                         environment_response=args.environment_response)
    result = simulate(read_weather_csv(args.weather_csv),
                      [date.fromisoformat(value) for value in args.cut],
                      soil, tuple(args.groups), config)
    write_results(result, args.output_dir)
    print(f"Simulated {len(result.daily)} days and {len(result.harvests)} harvests")


if __name__ == "__main__":
    main()
