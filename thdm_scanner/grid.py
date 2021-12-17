#!/usr/bin/env python

import logging

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


class THDMModel(object):

    def __init__(self, model_name, scan_pars, model_pars):
        self._name = model_name
        if len(scan_pars.keys()) != 2:
            raise ValueError("Too few or too less scan parameters given")
        self._scan_parameters = scan_pars
        self._model_parameters = model_pars
        self._model_points = []

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
        hists["model_validity"] = ROOT.TH2D("model_validity",
                                            "model_validity",
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
