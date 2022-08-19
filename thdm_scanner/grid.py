#!/usr/bin/env python

import logging
import math

import numpy as np

import ROOT
ROOT.gROOT.SetBatch()


logger = logging.getLogger(__name__)


class HiggsProperties(object):

    def __init__(self, name):
        self._name = name
        self._mass = 0.
        self._gg_xs = 0.
        self._bb_xs = 0.
        self._br_tautau = 0.
        self._gg_xs_scale_unc = (0., 0.)
        self._yukawa_t = 0
        self._yukawa_b = 0
        self._gg_xs_pdfas_unc = (0., 0.)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, mass):
        self._mass = mass

    @property
    def gg_xs(self):
        return self._gg_xs

    @gg_xs.setter
    def gg_xs(self, xsec):
        self._gg_xs = xsec

    @property
    def bb_xs(self):
        return self._bb_xs

    @bb_xs.setter
    def bb_xs(self, xsec):
        self._bb_xs = xsec

    @property
    def br_tautau(self):
        return self._br_tautau

    @br_tautau.setter
    def br_tautau(self, br):
        self._br_tautau = br

    @property
    def gg_xs_scale_unc(self):
        return self._gg_xs_scale_unc

    @gg_xs_scale_unc.setter
    def gg_xs_scale_unc(self, unc):
        if len(unc) != 2:
            raise ValueError("Uncertainty on gg{} xsec must be given as tuple.".format(self.name))  # noqa: E501
        self._gg_xs_scale_unc = unc

    @property
    def gt(self):
        return self._yukawa_t

    @gt.setter
    def gt(self, yukawa):
        self._yukawa_t = yukawa

    @property
    def gb(self):
        return self._yukawa_b

    @gb.setter
    def gb(self, yukawa):
        self._yukawa_b = yukawa

    @property
    def gg_xs_pdfas_unc(self):
        return self._gg_xs_pdfas_unc

    @gg_xs_pdfas_unc.setter
    def gg_xs_pdfas_unc(self, unc):
        if len(unc) != 2:
            raise ValueError("Uncertainty on gg{} xsec must be given as tuple.".format(self.name))  # noqa: E501
        self._gg_xs_pdfas_unc = unc


class THDMPoint(object):

    def __init__(self, par_names, par_values,
                 h=None, H=None, A=None):
        """Create class scoring Benchmark Parameter point.

        @param: par_names
            Name of the free parameters in the scan
        @param: par_values
            Parameter values of the free parameters for this specific point
        """
        if len(par_names) != 2 or len(par_values) != 2:
            raise ValueError("More than two names or parameters "
                             "specified")
        self.parameter_names = par_names
        self.parameter_point = par_values
        self._is_valid_model = False
        self._h = h
        self._H = H
        self._A = A
        self._cos_betal = 0.

    @property
    def is_valid_model(self):
        return self._is_valid_model

    @is_valid_model.setter
    def is_valid_model(self, model_valid):
        self._is_valid_model = model_valid

    @property
    def h(self):
        return self._h

    @h.setter
    def h(self, h_prop):
        self._h = h_prop

    @property
    def H(self):
        return self._H

    @H.setter
    def H(self, H_prop):
        self._H = H_prop

    @property
    def A(self):
        return self._A

    @A.setter
    def A(self, A_prop):
        self._A = A_prop

    @property
    def cos_betal(self):
        return self._cos_betal

    @cos_betal.setter
    def cos_betal(self, cos_betal):
        # Set cos(beta-alpha) always directly from inputs
        # We can be save that we use convention B then and
        # translate to sin(beta-alpha) here.
        # sin(beta-alpha) should not be set but only retrieved.
        self._cos_betal = cos_betal

    @property
    def sin_betal(self):
        return math.sin(math.acos(self._cos_betal))

    @sin_betal.setter
    def sin_betal(self, sin_betal):
        raise NotImplementedError("Passing a value for sin(beta-alpha) "
                                  "is only allow via the cos_betal member")


class THDMModel(object):

    def __init__(self, model_name, scan_pars, model_pars,
                 include_pdfas_unc=False):
        self._name = model_name
        if len(scan_pars.keys()) != 2:
            raise ValueError("Too few or too less scan parameters given")
        self._scan_parameters = scan_pars
        self._model_parameters = model_pars
        self._model_points = []
        self._pdfas_unc = include_pdfas_unc

    @property
    def name(self):
        return self._name

    @property
    def fixed_model_params(self):
        return self._model_parameters

    def parameter_points(self):
        par_1, par_2 = self._scan_parameters.keys()
        for x in np.round(np.arange(*self._scan_parameters[par_1]),
                          decimals=2):
            for y in np.round(np.arange(*self._scan_parameters[par_2]),
                              decimals=2):
                yield ((par_1, x), (par_2, y))

    def add_point(self, point):
        self._model_points.append(point)

    def write_to_root(self):
        # Quantities to be written:
        # Masses, cross sections and branching ratios for every Higgs boson
        # model validity
        logger.info("Writing created grid outputs to root file {}".format(
                        self.name + ".root"))
        # Open root file
        output = ROOT.TFile(self.name + ".root", "recreate")
        # Determine binning
        par_1, par_2 = self._scan_parameters.keys()
        xlow = (self._scan_parameters[par_1][0]
                - self._scan_parameters[par_1][2] / 2.)
        xup = (self._scan_parameters[par_1][1]
               - self._scan_parameters[par_1][2] / 2.)
        xbins = int(round((self._scan_parameters[par_1][1] - self._scan_parameters[par_1][0])  # noqa: E501
                          / self._scan_parameters[par_1][2]))
        ylow = (self._scan_parameters[par_2][0]
                - self._scan_parameters[par_2][2] / 2.)
        yup = (self._scan_parameters[par_2][1]
               - self._scan_parameters[par_2][2] / 2.)
        ybins = int(round((self._scan_parameters[par_2][1] - self._scan_parameters[par_2][0])  # noqa: E501
                          / self._scan_parameters[par_2][2]))
        # Create histograms for each quantity to be written
        hists = {}
        for boson in ["h", "H", "A"]:
            for quant in map(lambda x: x.format(boson),
                             ["m_{}",
                              "xs_gg{}",
                              "xs_bb{}",
                              "br_{}tautau",
                              "xs_gg{}_scale_down",
                              "xs_gg{}_scale_up",
                              "gt_{}",
                              "gb_{}"]):
                hists["{}".format(quant)] = ROOT.TH2D("{}".format(quant),
                                                      "{}".format(quant),
                                                      xbins,
                                                      xlow,
                                                      xup,
                                                      ybins,
                                                      ylow,
                                                      yup)
            if self._pdfas_unc:
                for quant in map(lambda x: x.format(boson),
                                 ["xs_gg{}_pdfas_down",
                                  "xs_gg{}_pdfas_up",
                                  ]):
                    hists["{}".format(quant)] = ROOT.TH2D("{}".format(quant),
                                                          "{}".format(quant),
                                                          xbins,
                                                          xlow,
                                                          xup,
                                                          ybins,
                                                          ylow,
                                                          yup)
        hists["model_validity"] = ROOT.TH2D("model_validity",
                                            "model_validity",
                                            xbins,
                                            xlow,
                                            xup,
                                            ybins,
                                            ylow,
                                            yup)
        hists["cos_betal"] = ROOT.TH2D("cos(beta-alpha)",
                                       "cos(beta-alpha)",
                                       xbins,
                                       xlow,
                                       xup,
                                       ybins,
                                       ylow,
                                       yup)
        hists["sin_betal"] = ROOT.TH2D("sin(beta-alpha)",
                                       "sin(beta-alpha)",
                                       xbins,
                                       xlow,
                                       xup,
                                       ybins,
                                       ylow,
                                       yup)
        # Fill the histograms per model point
        for point in self._model_points:
            x_val, y_val = (point.parameter_point
                            if point.parameter_names[0] == par_1
                            else reversed(point.parameter_point))
            hists["model_validity"].Fill(x_val, y_val, point.is_valid_model)
            hists["cos_betal"].Fill(x_val, y_val, point.cos_betal)
            hists["sin_betal"].Fill(x_val, y_val, point.sin_betal)
            for boson in ["h", "H", "A"]:
                hists["m_{}".format(boson)].Fill(x_val, y_val,
                                                 getattr(point, boson).mass)
                hists["xs_gg{}".format(boson)].Fill(
                        x_val, y_val,
                        getattr(point, boson).gg_xs)
                hists["xs_bb{}".format(boson)].Fill(
                        x_val, y_val,
                        getattr(point, boson).bb_xs)
                hists["br_{}tautau".format(boson)].Fill(
                        x_val, y_val,
                        getattr(point, boson).br_tautau)
                hists["xs_gg{}_scale_down".format(boson)].Fill(
                        x_val, y_val,
                        getattr(point, boson).gg_xs_scale_unc[0])
                hists["xs_gg{}_scale_up".format(boson)].Fill(
                        x_val, y_val,
                        getattr(point, boson).gg_xs_scale_unc[1])
                hists["gt_{}".format(boson)].Fill(
                        x_val, y_val,
                        getattr(point, boson).gt)
                hists["gb_{}".format(boson)].Fill(
                        x_val, y_val,
                        getattr(point, boson).gb)
                if self._pdfas_unc:
                    print("Filling histograms for uncertainties")
                    hists["xs_gg{}_pdfas_down".format(boson)].Fill(
                            x_val, y_val,
                            getattr(point, boson).gg_xs_pdfas_unc[0])
                    hists["xs_gg{}_pdfas_up".format(boson)].Fill(
                            x_val, y_val,
                            getattr(point, boson).gg_xs_pdfas_unc[1])

        # Write and close the root file
        output.Write()
        output.Close()


class THDMInput(object):

    def __init__(self,
                 mh=125.,
                 mH=200.,
                 cos_betal=0.0,
                 Z4=0.1,
                 Z5=0.1,
                 Z7=0.,
                 tanb=5.,
                 thdm_type=2):
        self._mh = mh
        self._mH = mH
        if isinstance(cos_betal, str):
            self._cos_betal = eval(cos_betal)
        else:
            self._cos_betal = cos_betal
        self._Z4 = Z4
        self._Z5 = Z5
        self._Z7 = Z7
        self._tanb = tanb
        self._type = thdm_type

    @property
    def mh(self):
        return self._mh

    @mh.setter
    def mh(self, mh):
        self._mh = mh

    @property
    def mH(self):
        return self._mH

    @mH.setter
    def mH(self, mH):
        self._mH = mH

    @property
    def cos_betal(self):
        return self._cos_betal

    @cos_betal.setter
    def cos_betal(self, cos_betal):
        self._cos_betal = cos_betal

    @property
    def sin_betal(self):
        # Assume that we always want to convert to the
        # SusHi and 2HDMC convention for values of beta-alpha
        # Add correct rescaling of angles, transformation of
        # conventions taken from https://sushi.hepforge.org/manual.html
        # 0 <= (beta-alpha)_B <= pi/2: no conversion of angle necessary
        # (beta-alpha)_B > pi/2: (beta-alpha)_B = (beta-alpha)_A - pi
        betal = math.acos(self._cos_betal)
        if betal > math.pi/2:
            betal -= math.pi
        sin_betal = math.sin(betal)
        return sin_betal

    @property
    def Z4(self):
        return self._Z4

    @Z4.setter
    def Z4(self, Z4):
        self._Z4 = Z4

    @property
    def Z5(self):
        return self._Z5

    @Z5.setter
    def Z5(self, Z5):
        self._Z5 = Z5

    @property
    def Z7(self):
        return self._Z7

    @Z7.setter
    def Z7(self, Z7):
        self._Z7 = Z7

    @property
    def tanb(self):
        return self._tanb

    @tanb.setter
    def tanb(self, tanb):
        self._tanb = tanb

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, thdm_type):
        self._type = thdm_type


class THDMPhysicsInput(object):

    def __init__(self,
                 mh=125.,
                 mH=200.,
                 mA=300.,
                 mHp=300.,
                 cos_betal=0.0,
                 lambda6=0.0,
                 lambda7=0.0,
                 tanb=5.,
                 m12_square=None,
                 thdm_type=2):
        # TODO: For now calculate m12_square hard coded
        # here
        self._mh = mh
        self._mH = mH
        if isinstance(mA, str):
            self._mA = eval(mA)
        else:
            self._mA = mA
        if isinstance(mHp, str):
            self._mHp = eval(mHp)
        else:
            self._mHp = mHp
        if isinstance(cos_betal, str):
            self._cos_betal = eval(cos_betal)
        else:
            self._cos_betal = cos_betal
        self._lambda6 = lambda6
        self._lambda7 = lambda7
        self._tanb = tanb
        self._type = thdm_type
        # Calculate m12_square from the given values
        if m12_square is None:
            Z5 = self._mH**2 * (1 - self._cos_betal**2) \
                 + mh**2 * self._cos_betal**2 \
                 - self._mA**2
            Z6 = (mh**2 - self._mH**2) \
                  * self._cos_betal*math.sin(math.acos(self._cos_betal))
            lambda5 = Z5 + 0.5 * Z6 * math.tan(2 * math.atan(tanb))
            logger.debug("Parameters in physical basis:")
            logger.debug("Z5: ", Z5)
            logger.debug("Z6: ", Z6)
            logger.debug("lamda5: ", lambda5)
            logger.debug("beta: ", math.atan(tanb))
            logger.debug("sin(beta-alpha): ", math.sin(math.acos(self._cos_betal)))
            self._m12_square = max(1-1./(tanb**2), 0) \
                                * 0.5 * math.sin(2 * math.atan(tanb))*(self._mA**2 + lambda5)
        elif isinstance(m12_square, str):
            self._m12_square = eval(m12_square)
        else:
            self._m12_square = m12_square

    @property
    def mh(self):
        return self._mh

    @mh.setter
    def mh(self, mh):
        self._mh = mh

    @property
    def mH(self):
        return self._mH

    @mH.setter
    def mH(self, mH):
        self._mH = mH

    @property
    def mA(self):
        return self._mA

    @mA.setter
    def mA(self, mA):
        self._mA = mA

    @property
    def mHp(self):
        return self._mHp

    @mHp.setter
    def mHp(self, mHp):
        self._mHp = mHp

    @property
    def cos_betal(self):
        return self._cos_betal

    @cos_betal.setter
    def cos_betal(self, cos_betal):
        self._cos_betal = cos_betal

    @property
    def sin_betal(self):
        # Assume that we always want to convert to the
        # SusHi and 2HDMC convention for values of beta-alpha
        # Add correct rescaling of angles, transformation of
        # conventions taken from https://sushi.hepforge.org/manual.html
        # 0 <= (beta-alpha)_B <= pi/2: no conversion of angle necessary
        # (beta-alpha)_B > pi/2: (beta-alpha)_B = (beta-alpha)_A - pi
        betal = math.acos(self._cos_betal)
        if betal > math.pi/2:
            betal -= math.pi
        sin_betal = math.sin(betal)
        return sin_betal

    @sin_betal.setter
    def sin_betal(self, sin_betal):
        # Convert sin_betal inputs such that we always store cos(bet-al)
        # in convention B
        betal = math.asin(sin_betal)
        if betal < 0:
            betal = betal + math.pi
        self._cos_betal = math.cos(betal)

    @property
    def lambda6(self):
        return self._lambda6

    @lambda6.setter
    def lambda6(self, lambda6):
        self._lambda6 = lambda6

    @property
    def lambda7(self):
        return self._lambda7

    @lambda7.setter
    def lambda7(self, lambda7):
        self._lambda7 = lambda7

    @property
    def tanb(self):
        return self._tanb

    @tanb.setter
    def tanb(self, tanb):
        self._tanb = tanb

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, thdm_type):
        self._type = thdm_type

    @property
    def m12_square(self):
        return self._m12_square

    @m12_square.setter
    def m12_square(self, m12_square):
        self._m12_square = m12_square
