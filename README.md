# ModVege-WN

A small, installable Python implementation of a daily pasture growth model
adapted from ModVege (Jouven, Carrère and Baumont, 2006). It estimates four above-ground dry-matter
pools, leaf area, soil-water status, and biomass removed by cutting.

**Status:** experimental implementation, version 0.1.0. The package has not yet
been independently evaluated against observational data. It uses a **fixed
ModVege nutrition index (NI)** for each simulation.

## Install and run

Python 3.10 or newer is required. The simulation has no third-party runtime
dependencies.

```powershell
py -m pip install -e .
py -m unittest discover -s tests -v
py -m modvege_wn examples/weather.csv --cut 2020-01-25 --soil 20 40 20 1.0 --groups 1 0 0 0 --ni 0.85 --output-dir outputs/example
```

If `py` is unavailable, replace it with the path to your Python executable.
The example weather and soil inputs are **illustrative** and test the interface
only. `daily.csv`
reports the standing crop at the end of each day's growth; on cut dates, that
value is recorded immediately **before** biomass is removed. `harvests.csv`
contains the removed biomass by pool and the total yield. Mass is kg dry matter
per hectare. `outputs/` is ignored by Git.

## Python API

```python
from datetime import date
from modvege_wn import ModelConfig, Soil, simulate
from modvege_wn.io import read_weather_csv

result = simulate(
    weather=read_weather_csv("examples/weather.csv"),
    cuts=[date(2020, 1, 25)],
    soil=Soil(clay_pct=20, sand_pct=40, organic_matter_g_kg=20, root_depth_m=1.0),
    sward_proportions=(1.0, 0.0, 0.0, 0.0),  # groups A, B, C, D
    config=ModelConfig(nitrogen_index=.85, groundwater_floor=0.0),
)
print(result.harvests[0].yield_kg_ha)
```

## Input data

Weather CSV columns (one row per consecutive day, starting on 1 January):

| Column | Unit |
| --- | --- |
| `date` | ISO date, `YYYY-MM-DD` |
| `max_temp_c`, `min_temp_c` | degrees Celsius |
| `reference_et_mm` | reference evapotranspiration, mm/day |
| `solar_radiation_mj_m2` | incoming solar radiation, MJ/m²/day |
| `precipitation_mm` | precipitation, mm/day |

The four sward proportions must be nonnegative and sum to one. Soil texture is
percent, organic matter is g/kg, and root depth is metres. NI is restricted to
0.35–1.20. `groundwater_floor` is a fraction of available water capacity; it is
an empirical lower bound, not a groundwater-flow calculation.

The defaults use FAO-56 style water stress and the minimum of nitrogen and
water stress.
Other original options are available through `ModelConfig` and the CLI.

## Model provenance and limits

The original 2021 `ModVege.py` is preserved in [`legacy/`](legacy/ModVege.py).
[`docs/model-notes.md`](docs/model-notes.md) records implementation changes and
limitations. The original model description is Jouven et al. (2006), *Grass and Forage Science* 61,
112–124, [doi:10.1111/j.1365-2494.2006.00515.x](https://doi.org/10.1111/j.1365-2494.2006.00515.x).

The example inputs are illustrative; no observational dataset is included.
Select an open-source license with the code authors before publishing this
repository on GitHub.
