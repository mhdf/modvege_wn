"""Small, unit-preserving summary metrics for observed and simulated yields."""

from math import sqrt
from statistics import mean


def yield_metrics(observed_kg_ha: list[float], predicted_kg_ha: list[float]) -> dict[str, float]:
    """Return count, mean bias, mean absolute error and RMSE in kg DM ha-1."""
    if not observed_kg_ha or len(observed_kg_ha) != len(predicted_kg_ha):
        raise ValueError("observed and predicted yields need equal, nonzero lengths")
    errors = [predicted-observed for observed, predicted in zip(observed_kg_ha, predicted_kg_ha)]
    return {
        "n": len(errors),
        "bias_kg_ha": mean(errors),
        "mae_kg_ha": mean(abs(error) for error in errors),
        "rmse_kg_ha": sqrt(mean(error*error for error in errors)),
    }
