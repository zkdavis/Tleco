import unittest
import numpy as np
import paramo as pp

class TestPyParamo(unittest.TestCase):
    def test_v_rela_s(self):
        """
        Test relativistic speed function with scalar argument
        """
        beta = 0.1
        gamma = 1.005
        self.assertAlmostEqual(pp.gofb(beta), gamma, 3)

    def test_v_rela_v(self):
        """
        Test relativistic speed function with array argument
        """
        beta = np.array([0.1, 0.25, 0.5])
        gamma = np.array([1.005, 1.033, 1.155])
        self.assertIsNone(np.testing.assert_allclose(pp.gofb(beta), gamma, atol=1))

if __name__ == "__main__":
    unittest.main()
