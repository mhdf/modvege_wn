"""Nitrogen is currently represented by a fixed ModVege nutrition index."""


def environment_multiplier(nitrogen_index: float, water: float, par: float, temperature: float, mode: str) -> float:
    """Combine stress factors as in the two alternatives in the research code."""
    if mode == "minimum":
        nitrogen_water = min(nitrogen_index, water)
    elif mode == "product":
        nitrogen_water = nitrogen_index * water
    else:
        raise ValueError(f"unknown environment response: {mode}")
    return nitrogen_water * par * temperature
