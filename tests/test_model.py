import unittest
from datetime import date, timedelta

from modvege_wn import ModelConfig, Soil, WeatherDay, functional_traits, simulate
from modvege_wn.water import holding_capacity_mm


SOIL = Soil(20, 40, 20, 1.0)
GROUPS = (1.0, 0.0, 0.0, 0.0)


def weather(days: int, maximum: float = 20, minimum: float = 10,
            radiation: float = 18, et: float = 2.0, rain: float = 2.0):
    start = date(2020, 1, 1)
    return [WeatherDay(start + timedelta(days=i), maximum, minimum, et, radiation, rain)
            for i in range(days)]


class ModelTests(unittest.TestCase):
    def test_original_group_a_traits_and_water_capacity(self):
        traits = functional_traits(GROUPS)
        self.assertEqual(traits.sla_m2_g, .033)
        self.assertEqual(traits.st1_degree_days, 600)
        self.assertAlmostEqual(holding_capacity_mm(SOIL), 151.8)

    def test_each_day_simulated_once_across_cut(self):
        days = weather(35)
        cut_date = date(2020, 1, 20)
        result = simulate(days, [cut_date], SOIL, GROUPS)
        uncut = simulate(days, [], SOIL, GROUPS)
        self.assertEqual(len(result.daily), 35)
        self.assertEqual([r.date for r in result.daily], [d.date for d in days])
        self.assertEqual(len(result.harvests), 1)
        self.assertAlmostEqual(result.daily[19].standing_biomass_kg_ha,
                               uncut.daily[19].standing_biomass_kg_ha)
        self.assertGreater(result.harvests[0].yield_kg_ha, 0)
        self.assertLess(result.daily[20].standing_biomass_kg_ha,
                        uncut.daily[20].standing_biomass_kg_ha)

    def test_cold_zero_radiation_keeps_biomass_finite_and_nonnegative(self):
        days = weather(12, maximum=0, minimum=-8, radiation=0)
        result = simulate(days, [], SOIL, GROUPS)
        for row in result.daily:
            self.assertGreaterEqual(row.standing_biomass_kg_ha, 0)
            self.assertEqual(row.growth_kg_ha, 0)
            self.assertGreaterEqual(row.relative_soil_water, 0)
            self.assertLessEqual(row.relative_soil_water, 1)

    def test_missing_weather_date_is_rejected(self):
        days = weather(3)
        with self.assertRaisesRegex(ValueError, "consecutive"):
            simulate([days[0], days[2]], [], SOIL, GROUPS)

    def test_options_change_environmental_response(self):
        days = weather(50, et=6.0, rain=0.0)
        minimum = simulate(days, [], SOIL, GROUPS,
                           ModelConfig(environment_response="minimum"))
        product = simulate(days, [], SOIL, GROUPS,
                           ModelConfig(environment_response="product"))
        self.assertNotEqual(minimum.daily[-1].standing_biomass_kg_ha,
                            product.daily[-1].standing_biomass_kg_ha)


if __name__ == "__main__":
    unittest.main()
