"""Daily ModVege-style simulation with explicit cutting events and outputs.

The equations follow ``legacy/ModVege.py``. Changes needed for a reusable
simulation are recorded in README.md and docs/model-notes.md.
"""

from bisect import bisect_left
from dataclasses import dataclass
from datetime import date, timedelta
from math import exp, isfinite
from typing import Iterable

from .nitrogen import environment_multiplier
from .parameters import ModelConfig, Soil, SwardTraits, functional_traits
from .water import holding_capacity_mm, water_stress


@dataclass(frozen=True)
class WeatherDay:
    date: date
    max_temp_c: float
    min_temp_c: float
    reference_et_mm: float
    solar_radiation_mj_m2: float
    precipitation_mm: float

    def __post_init__(self) -> None:
        values = (self.max_temp_c, self.min_temp_c, self.reference_et_mm,
                  self.solar_radiation_mj_m2, self.precipitation_mm)
        if not all(isfinite(v) for v in values):
            raise ValueError(f"weather values must be finite on {self.date}")
        if self.max_temp_c < self.min_temp_c:
            raise ValueError(f"maximum temperature below minimum on {self.date}")
        if any(v < 0 for v in values[2:]):
            raise ValueError(f"ET, radiation and precipitation must be nonnegative on {self.date}")


@dataclass(frozen=True)
class DailyResult:
    date: date
    green_vegetative_kg_ha: float
    green_reproductive_kg_ha: float
    dead_vegetative_kg_ha: float
    dead_reproductive_kg_ha: float
    lai_m2_m2: float
    relative_soil_water: float
    crop_et_mm: float
    par_mj_m2: float
    absorbed_par_fraction: float
    thermal_time_degree_days: float
    seasonal_effect: float
    water_stress: float
    temperature_stress: float
    par_stress: float
    environment_multiplier: float
    growth_kg_ha: float
    senescence_kg_ha: float
    abscission_kg_ha: float
    omd_green_vegetative: float
    omd_green_reproductive: float

    @property
    def standing_biomass_kg_ha(self) -> float:
        return (self.green_vegetative_kg_ha + self.green_reproductive_kg_ha
                + self.dead_vegetative_kg_ha + self.dead_reproductive_kg_ha)


@dataclass(frozen=True)
class HarvestResult:
    date: date
    green_vegetative_kg_ha: float
    green_reproductive_kg_ha: float
    dead_vegetative_kg_ha: float
    dead_reproductive_kg_ha: float

    @property
    def yield_kg_ha(self) -> float:
        return (self.green_vegetative_kg_ha + self.green_reproductive_kg_ha
                + self.dead_vegetative_kg_ha + self.dead_reproductive_kg_ha)


@dataclass(frozen=True)
class SimulationResult:
    daily: tuple[DailyResult, ...]
    harvests: tuple[HarvestResult, ...]


@dataclass
class _State:
    gv: float
    gr: float
    dv: float
    dr: float
    age_gv: float = 400.0
    age_gr: float = 0.0
    age_dv: float = 500.0
    age_dr: float = 0.0
    water: float = 1.0


def _linear(x: float, points: tuple[tuple[float, float], ...]) -> float:
    if x <= points[0][0]:
        return points[0][1]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x <= x1:
            return y0 + (x-x0) * (y1-y0) / (x1-x0)
    return points[-1][1]


def _seasonal(thermal: float, traits: SwardTraits) -> float:
    return _linear(thermal, (
        (0, traits.seasonal_min), (200, traits.seasonal_min),
        (traits.st1_degree_days-200, traits.seasonal_max),
        (traits.st1_degree_days-100, traits.seasonal_max),
        (traits.st2_degree_days, traits.seasonal_min),
        (15*traits.st2_degree_days, traits.seasonal_min),
    ))


def _senescence_age_factor(age: float, lifespan: float) -> float:
    return _linear(age, ((0, 1), (lifespan/3, 1), (2*lifespan/3, 2),
                         (lifespan, 3), (15*lifespan, 3)))


def _abscission_age_factor(age: float, lifespan: float) -> float:
    epsilon = .01
    return _linear(age, ((0, 1), (lifespan/3-epsilon, 1),
                         (lifespan/3, 2), (2*lifespan/3-epsilon, 2),
                         (2*lifespan/3, 3), (lifespan, 3),
                         (15*lifespan, 3)))


def _average_age(old_mass: float, loss: float, addition: float, age: float, temperature: float) -> float:
    """Age retained material, then dilute its mean age with new material."""
    remaining = max(0.0, old_mass-loss)
    total = remaining + addition
    return remaining * (age + max(temperature, 0.0)) / total if total > 0 else 0.0


def _validate_dates(weather: list[WeatherDay], cuts: list[date]) -> None:
    if not weather:
        raise ValueError("at least one weather day is required")
    if weather[0].date.month != 1 or weather[0].date.day != 1:
        raise ValueError("weather must begin on 1 January so thermal time is defined")
    for previous, current in zip(weather, weather[1:]):
        if current.date != previous.date + timedelta(days=1):
            raise ValueError("weather must contain exactly one row per consecutive day")
    if len(cuts) != len(set(cuts)):
        raise ValueError("cut dates must be unique")
    if any(d < weather[0].date or d > weather[-1].date for d in cuts):
        raise ValueError("all cuts must lie within the weather period")


def simulate(
    weather: Iterable[WeatherDay],
    cuts: Iterable[date],
    soil: Soil,
    sward_proportions: tuple[float, float, float, float],
    config: ModelConfig = ModelConfig(),
) -> SimulationResult:
    """Simulate daily pasture state and cutting yield.

    Each cut is applied *after* that day's growth. Daily output on a cut date
    therefore reports standing biomass immediately before removal. A fixed
    nitrogen index applies to the whole run in this initial package release.
    """
    days = list(weather)
    cut_dates = sorted(cuts)
    _validate_dates(days, cut_dates)
    traits = functional_traits(sward_proportions)
    capacity_mm = holding_capacity_mm(soil)
    residual = {
        "gv": config.residual_height_m * 10 * traits.bd_gv_g_m3,
        "gr": config.residual_height_m * 10 * traits.bd_gr_g_m3,
        "dv": config.residual_height_m * 10 * traits.bd_dv_g_m3,
        "dr": config.residual_height_m * 10 * traits.bd_dr_g_m3,
    }
    state = _State(gv=residual["gv"], gr=0.0, dv=residual["dv"], dr=0.0)
    daily: list[DailyResult] = []
    harvests: list[HarvestResult] = []
    thermal = 0.0
    reproductive_done = False
    reproductive_this_cut = False
    current_year = days[0].date.year
    reproductive_lifespan = traits.st2_degree_days - traits.st1_degree_days

    for day in days:
        if day.date.year != current_year:
            current_year = day.date.year
            thermal = 0.0
            reproductive_done = False
            reproductive_this_cut = False

        temperature = (day.max_temp_c + day.min_temp_c) / 2
        thermal += max(temperature - 3.0, 0.0)
        seasonal = _seasonal(thermal, traits)
        temperature_stress = max(0.0, min(1.0, (temperature-3.0)/7.0))
        par = .47 * day.solar_radiation_mj_m2
        par_stress = max(0.0, min(1.0, 1-.0445*(par-5)))
        lai_before = max(0.0, traits.sla_m2_g * state.gv/10 * traits.lamina_fraction)
        absorbed_par = 1 - exp(-.6*lai_before)
        crop_et_potential = 1.15 * day.reference_et_mm
        water_factor = water_stress(state.water, crop_et_potential, config.water_response)
        potential_transpiration = crop_et_potential * absorbed_par
        crop_et = potential_transpiration * water_factor + (crop_et_potential-potential_transpiration)*state.water
        environment = environment_multiplier(config.nitrogen_index, water_factor,
                                             par_stress, temperature_stress,
                                             config.environment_response)
        growth = max(0.0, par * 3.0 * absorbed_par * 10 * environment * seasonal)

        cut_index = bisect_left(cut_dates, day.date)
        next_cut = cut_dates[cut_index] if cut_index < len(cut_dates) else None
        if not reproductive_done and traits.st1_degree_days <= thermal <= traits.st2_degree_days:
            partition = .25 + .75*(config.nitrogen_index-.35)/.65
            if next_cut is not None and (next_cut-day.date).days >= config.min_reproductive_days:
                reproductive_this_cut = True
        else:
            partition = 0.0
        growth_gr = growth * partition
        growth_gv = growth - growth_gr

        if temperature > 3:
            sen_gv = min(state.gv, .002*state.gv*temperature *
                         _senescence_age_factor(state.age_gv, traits.leaf_lifespan_degree_days))
            sen_gr = min(state.gr, .001*state.gr*temperature *
                         _senescence_age_factor(state.age_gr, reproductive_lifespan))
        elif temperature < 0:
            sen_gv = min(state.gv, .002*state.gv*abs(temperature))
            sen_gr = min(state.gr, .001*state.gr*abs(temperature))
        else:
            sen_gv = sen_gr = 0.0
        if temperature > 0:
            abs_dv = min(state.dv, .001*state.dv*temperature *
                         _abscission_age_factor(state.age_dv, traits.leaf_lifespan_degree_days))
            abs_dr = min(state.dr, .0005*state.dr*temperature *
                         _abscission_age_factor(state.age_dr, reproductive_lifespan))
        else:
            abs_dv = abs_dr = 0.0
        dead_from_gv = .6 * sen_gv
        dead_from_gr = .8 * sen_gr
        state.age_gv = _average_age(state.gv, sen_gv, growth_gv, state.age_gv, temperature)
        state.age_gr = _average_age(state.gr, sen_gr, growth_gr, state.age_gr, temperature)
        state.age_dv = _average_age(state.dv, abs_dv, dead_from_gv, state.age_dv, temperature)
        state.age_dr = _average_age(state.dr, abs_dr, dead_from_gr, state.age_dr, temperature)
        state.gv += growth_gv - sen_gv
        state.gr += growth_gr - sen_gr
        state.dv += dead_from_gv - abs_dv
        state.dr += dead_from_gr - abs_dr
        state.water = max(config.groundwater_floor,
                          min(1.0, state.water + (day.precipitation_mm-crop_et)/capacity_mm))
        lai_after = max(0.0, traits.sla_m2_g * state.gv/10 * traits.lamina_fraction)
        omd_gv = max(traits.omd_min_gv, traits.omd_max_gv -
                     state.age_gv*(traits.omd_max_gv-traits.omd_min_gv)/traits.leaf_lifespan_degree_days)
        omd_gr = max(traits.omd_min_gr, traits.omd_max_gr -
                     state.age_gr*(traits.omd_max_gr-traits.omd_min_gr)/reproductive_lifespan)
        daily.append(DailyResult(
            date=day.date, green_vegetative_kg_ha=state.gv,
            green_reproductive_kg_ha=state.gr, dead_vegetative_kg_ha=state.dv,
            dead_reproductive_kg_ha=state.dr, lai_m2_m2=lai_after,
            relative_soil_water=state.water, crop_et_mm=crop_et, par_mj_m2=par,
            absorbed_par_fraction=absorbed_par, thermal_time_degree_days=thermal,
            seasonal_effect=seasonal, water_stress=water_factor,
            temperature_stress=temperature_stress, par_stress=par_stress,
            environment_multiplier=environment, growth_kg_ha=growth,
            senescence_kg_ha=sen_gv+sen_gr, abscission_kg_ha=abs_dv+abs_dr,
            omd_green_vegetative=omd_gv, omd_green_reproductive=omd_gr,
        ))
        if next_cut == day.date:
            removed = {name: max(0.0, getattr(state, name)-limit)
                       for name, limit in residual.items()}
            harvests.append(HarvestResult(day.date, removed["gv"], removed["gr"],
                                          removed["dv"], removed["dr"]))
            for name, mass in removed.items():
                setattr(state, name, getattr(state, name)-mass)
            reproductive_done = reproductive_done or reproductive_this_cut
            reproductive_this_cut = False
    return SimulationResult(tuple(daily), tuple(harvests))
