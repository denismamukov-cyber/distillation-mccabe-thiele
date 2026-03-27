import unittest

from src.calculations import InputData, calculate, y_equilibrium


class CalculationTests(unittest.TestCase):
    def test_equilibrium_monotonic_sample(self):
        alpha = 2.4
        y1 = y_equilibrium(0.2, alpha)
        y2 = y_equilibrium(0.6, alpha)
        self.assertGreater(y2, y1)

    def test_material_balance_closure(self):
        data = InputData(
            feed_flow_kmol_h=100.0,
            feed_composition_light=0.45,
            distillate_composition_light=0.95,
            bottoms_composition_light=0.05,
            relative_volatility=2.4,
            reflux_ratio_factor_to_rmin=1.5,
            feed_thermal_condition_q=1.0,
            vapor_density_kg_m3=1.6,
            souders_brown_factor_m_s=1.0,
        )
        result = calculate(data)
        self.assertAlmostEqual(result.distillate_flow_kmol_h + result.bottoms_flow_kmol_h, data.feed_flow_kmol_h, places=9)
        self.assertGreater(result.n_theoretical, 0)


if __name__ == "__main__":
    unittest.main()
