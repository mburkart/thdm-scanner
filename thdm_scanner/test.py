#!/usr/bin/env python

import unittest
import math

from thdm_scanner import lha_utils
from thdm_scanner.utility import utils


class DToSConversionUtility(unittest.TestCase):

    def test_double_conversion(self):
        self.assertEqual(utils.d_to_fortran_s(0.1), "0.1d0")
        self.assertEqual(utils.d_to_fortran_s(125.), "125.d0")

    def test_string_conversion(self):
        self.assertEqual(utils.fortran_s_to_d("0.1d0"), 0.1)
        self.assertEqual(utils.fortran_s_to_d("125.d0"), 125.)
        self.assertEqual(utils.fortran_s_to_d("125.d-2"), 1.25)
        self.assertEqual(utils.fortran_s_to_d("1.25d2"), 125)

    def test_trim_zeros(self):
        self.assertEqual(utils.trim_zeros_from_dstr("125.0"), "125.")
        self.assertEqual(utils.trim_zeros_from_dstr("125.0000"), "125.")
        self.assertEqual(utils.trim_zeros_from_dstr("125.05"), "125.05")
        self.assertEqual(utils.trim_zeros_from_dstr("125.050"), "125.05")


class SusHiInputTestCase(unittest.TestCase):

    def setUp(self):
        self.sushi_input = lha_utils.SusHiInput("test_sushi_input.out")

    def test_higgs_boson_getter(self):
        self.assertEqual(self.sushi_input.higgs_boson, 11,
                         "incorrect Higgs boson stored")

    def test_higgs_boson_setter(self):
        self.sushi_input.higgs_boson = 12
        self.assertEqual(self.sushi_input.higgs_boson, 12,
                         "wrong Higgs boson after setting with int")
        self.sushi_input.higgs_boson = "21"
        self.assertEqual(self.sushi_input.higgs_boson, 21,
                         "wrong Higgs boson after setting with string")

    def test_higgs_boson_setter_oor(self):
        with self.assertRaises(ValueError):
            self.sushi_input.higgs_boson = 15

    def test_thdm_type_getter(self):
        self.assertEqual(self.sushi_input.thdm_type, 2,
                         "incorrect THDM type stored")

    def test_thdm_type_setter(self):
        self.sushi_input.thdm_type = 1
        self.assertEqual(self.sushi_input.thdm_type, 1,
                         "wrong THDM type after setting with int")
        self.sushi_input.thdm_type = "2"
        self.assertEqual(self.sushi_input.thdm_type, 2,
                         "wrong THDM type after setting with string")

    def test_thdm_type_setter_oor(self):
        with self.assertRaises(ValueError):
            self.sushi_input.thdm_type = 15

    def test_tanb_getter(self):
        self.assertEqual(self.sushi_input.tanb, 10.,
                         "incorrect tanb stored")

    def test_tanb_setter(self):
        self.sushi_input.tanb = 2.
        self.assertEqual(self.sushi_input.tanb, 2.,
                         "wrong tanb after setting")

    def test_mh_getter(self):
        self.assertEqual(self.sushi_input.mh, 125.,
                         "incorrect mh stored")

    def test_mh_setter(self):
        self.sushi_input.mh = 125.38
        self.assertEqual(self.sushi_input.mh, 125.38,
                         "wrong mh after setting")

    def test_mH_getter(self):
        self.assertEqual(self.sushi_input.mH, 200.,
                         "incorrect mH stored")

    def test_mH_setter(self):
        self.sushi_input.mH = 250.
        self.assertEqual(self.sushi_input.mH, 250.,
                         "wrong mH after setting")

    def test_sin_betal_getter(self):
        self.assertEqual(self.sushi_input.sin_betal, 0.5,
                         "incorrect sin_betal stored")

    def test_sin_betal_setter(self):
        self.sushi_input.sin_betal = 0.7
        self.assertAlmostEqual(self.sushi_input.sin_betal, 0.7,
                               msg="wrong sin_betal after setting")

    def test_cos_betal_getter(self):
        self.assertAlmostEqual(self.sushi_input.cos_betal,
                               0.8660254037844386,
                               "incorrect cos_betal stored")

    def test_cos_betal_setter(self):
        self.sushi_input.cos_betal = 0.7
        self.assertAlmostEqual(self.sushi_input.cos_betal, 0.7,
                               msg="wrong cos_betal after setting")

    def test_cos_betal_getter_from_sin_betal(self):
        self.sushi_input.sin_betal = math.sin(math.pi/4)
        self.assertAlmostEqual(self.sushi_input.cos_betal, math.cos(math.pi/4),
                               msg="wrong cos_betal after setting sin_betal")

    def test_Z4_getter(self):
        self.assertEqual(self.sushi_input.Z4, 0.1,
                         "incorrect Z4 stored")

    def test_Z4_setter(self):
        self.sushi_input.Z4 = 0.3
        self.assertEqual(self.sushi_input.Z4, 0.3,
                         "wrong Z4 after setting")

    def test_Z5_getter(self):
        self.assertEqual(self.sushi_input.Z5, 0.1,
                         "incorrect Z5 stored")

    def test_Z5_setter(self):
        self.sushi_input.Z5 = 0.2
        self.assertEqual(self.sushi_input.Z5, 0.2,
                         "wrong Z5 after setting")

    def test_Z7_getter(self):
        self.assertEqual(self.sushi_input.Z7, 0.1,
                         "incorrect Z7 stored")

    def test_Z7_setter(self):
        self.sushi_input.Z7 = 0.
        self.assertEqual(self.sushi_input.Z7, 0.,
                         "wrong Z7 after setting")


class THDMCOutputTestCase(unittest.TestCase):

    def setUp(self):
        self.thdmc_out = lha_utils.THDMCOutput("/work/mburkart/Run2Legacy/2HDM_Theory/THDM_Parameter_Scans/Test_CalcHybrid_Test.out")  # noqa: E501

    def test_valid_model(self):
        self.assertTrue(self.thdmc_out.valid_model)

    def test_mh(self):
        self.assertEqual(self.thdmc_out.mh, 125.)

    def test_mH(self):
        self.assertEqual(self.thdmc_out.mH, 200.)

    def test_mA(self):
        self.assertEqual(self.thdmc_out.mA, 1.66864595e+02)

    def test_mHp(self):
        self.assertEqual(self.thdmc_out.mHp, 1.66864595e+02)

    def test_br_htautau(self):
        self.assertEqual(self.thdmc_out.br_htautau, 8.66264746e-02)

    def test_br_Htautau(self):
        self.assertEqual(self.thdmc_out.br_Htautau, 4.95716970e-02)

    def test_br_Atautau(self):
        self.assertEqual(self.thdmc_out.br_Atautau, 9.28157342e-02)


class SusHiOutputTestCase(unittest.TestCase):

    def setUp(self):
        self.sushi_out = lha_utils.SusHiOutput("/work/mburkart/Run2Legacy/2HDM_Theory/SusHi-1.7.0/bin/2HDMC_H2basis.out")  # noqa: E501

    def test_xs_ggPhi(self):
        self.assertEqual(self.sushi_out.xs_ggPhi, 4.80455095E+01)

    def test_xs_bbPhi(self):
        self.assertEqual(self.sushi_out.xs_bbPhi, 3.50229376E+01)

    def test_xs_ggPhi_scale(self):
        self.assertEqual(self.sushi_out.xs_ggPhi_scale_up, 4.78036642E+00)
        self.assertEqual(self.sushi_out.xs_ggPhi_scale_down, -5.07690730E+00)

    def test_mPhi(self):
        self.assertEqual(self.sushi_out.mPhi, 1.25000000E+02)


if __name__ == "__main__":
    unittest.main()
