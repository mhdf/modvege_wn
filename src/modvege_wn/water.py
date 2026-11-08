"""Simple root-zone water balance functions inherited from the research code."""

from .parameters import Soil


def holding_capacity_mm(soil: Soil) -> float:
    """Available water in the root zone from the original pedotransfer equations."""
    om_pct = soil.organic_matter_g_kg / 10
    field_capacity = .2576 - .002 * soil.sand_pct + .0036 * soil.clay_pct + .0299 * om_pct
    wilting_point = .026 + .005 * soil.clay_pct + .0158 * om_pct
    capacity = (field_capacity - wilting_point) * soil.root_depth_m * 1000
    if capacity <= 0:
        raise ValueError("soil properties produce nonpositive available water capacity")
    return capacity


def water_stress(relative_water: float, crop_et_mm: float, response: str) -> float:
    """Return the water multiplier, constrained to the physical interval [0, 1]."""
    w = max(0.0, min(1.0, relative_water))
    if response == "fao56":
        depletion_fraction = .6 + .04 * (5 - crop_et_mm)
        critical_water = 1 - depletion_fraction
        return max(0.0, min(1.0, w / critical_water))
    if response == "calanca":
        if crop_et_mm <= 3.81:
            value = 4*w if w <= .2 else (.75*w + .65 if w <= .4 else (.25*w + .85 if w <= .6 else 1))
        elif crop_et_mm <= 6.35:
            value = 2*w if w <= .2 else (1.5*w + .1 if w <= .4 else (w + .3 if w <= .6 else (.5*w + .6 if w <= .8 else 1)))
        else:
            value = w
        return max(0.0, min(1.0, value))
    raise ValueError(f"unknown water response: {response}")
