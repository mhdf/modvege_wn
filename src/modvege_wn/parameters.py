"""Site and sward parameters. Biomass uses kg dry matter per hectare."""

from dataclasses import dataclass
from math import isclose, isfinite


@dataclass(frozen=True)
class Soil:
    clay_pct: float
    sand_pct: float
    organic_matter_g_kg: float
    root_depth_m: float

    def __post_init__(self) -> None:
        if not 0 <= self.clay_pct <= 100 or not 0 <= self.sand_pct <= 100:
            raise ValueError("clay_pct and sand_pct must be percentages from 0 to 100")
        if self.clay_pct + self.sand_pct > 100:
            raise ValueError("clay_pct + sand_pct cannot exceed 100")
        if self.organic_matter_g_kg < 0 or self.root_depth_m <= 0:
            raise ValueError("organic matter must be nonnegative and root depth positive")
        if not all(isfinite(v) for v in (self.clay_pct, self.sand_pct, self.organic_matter_g_kg, self.root_depth_m)):
            raise ValueError("soil parameters must be finite")


@dataclass(frozen=True)
class ModelConfig:
    nitrogen_index: float = 0.85
    groundwater_floor: float = 0.0
    water_response: str = "fao56"
    environment_response: str = "minimum"
    min_reproductive_days: int = 10
    residual_height_m: float = 0.05

    def __post_init__(self) -> None:
        if not 0.35 <= self.nitrogen_index <= 1.20:
            raise ValueError("nitrogen_index must be from 0.35 to 1.20")
        if not 0 <= self.groundwater_floor <= 1:
            raise ValueError("groundwater_floor must be from 0 to 1")
        if self.water_response not in {"fao56", "calanca"}:
            raise ValueError("water_response must be 'fao56' or 'calanca'")
        if self.environment_response not in {"minimum", "product"}:
            raise ValueError("environment_response must be 'minimum' or 'product'")
        if self.min_reproductive_days < 0 or self.residual_height_m < 0:
            raise ValueError("reproductive days and residual height must be nonnegative")
        if not all(isfinite(v) for v in (self.nitrogen_index, self.groundwater_floor, self.residual_height_m)):
            raise ValueError("model parameters must be finite")


@dataclass(frozen=True)
class SwardTraits:
    sla_m2_g: float
    lamina_fraction: float
    st1_degree_days: float
    st2_degree_days: float
    seasonal_max: float
    seasonal_min: float
    leaf_lifespan_degree_days: float
    omd_max_gv: float
    omd_min_gv: float
    omd_max_gr: float
    omd_min_gr: float
    bd_gv_g_m3: float
    bd_dv_g_m3: float
    bd_gr_g_m3: float
    bd_dr_g_m3: float


_GROUPS = {
    "sla_m2_g": (.033, .025, .022, .019),
    "lamina_fraction": (.680, .680, .680, .680),
    "st1_degree_days": (600, 700, 850, 1000),
    "st2_degree_days": (1200, 1350, 1550, 1850),
    "seasonal_max": (1.20, 1.30, 1.40, 1.50),
    "seasonal_min": (.80, .70, .60, .50),
    "leaf_lifespan_degree_days": (500, 800, 900, 1400),
    "omd_max_gv": (.90, .90, .85, .75),
    "omd_min_gv": (.75, .60, .65, .65),
    "omd_max_gr": (.90, .90, .85, .75),
    "omd_min_gr": (.65, .45, .45, .45),
    "bd_gv_g_m3": (850, 850, 1200, 800),
    "bd_dv_g_m3": (500, 500, 1800, 2200),
    "bd_gr_g_m3": (300, 300, 200, 150),
    "bd_dr_g_m3": (150, 150, 300, 450),
}


def functional_traits(proportions: tuple[float, float, float, float]) -> SwardTraits:
    """Mix the original functional groups A, B, C and D by their proportions."""
    if len(proportions) != 4 or any(not isfinite(v) or v < 0 for v in proportions):
        raise ValueError("provide four finite, nonnegative functional-group proportions")
    if not isclose(sum(proportions), 1.0, rel_tol=0, abs_tol=1e-8):
        raise ValueError("functional-group proportions must sum to one")
    return SwardTraits(**{
        name: sum(value * weight for value, weight in zip(values, proportions))
        for name, values in _GROUPS.items()
    })
