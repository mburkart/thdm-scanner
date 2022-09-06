#!/usr/bin/env python

from abc import ABCMeta, abstractmethod
import os
import subprocess
import multiprocessing
import logging
import math

import numpy as np

import thdm_scanner
import thdm_scanner.lha_utils as lha_utils


logger = logging.getLogger(__name__)


class THDMHarvesterABC(object):
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def set_inputs(self):
        pass

    @abstractmethod
    def harvest_output(self, model_point):
        pass


class THDMCHarvester(THDMHarvesterABC):
    """Collect results from 2HDMC calculations"""

    def __init__(self, outpath="output",
                 scan_parameters=("mH", "tanb")):
        self._outputfile = os.path.join(outpath, "2HDMC_output.out")
        self._scan_params = scan_parameters

    def set_inputs(self, inputs):
        self._input = inputs

    def harvest_output(self, model_point):
        outname = self._outputfile.replace(".out", ".{}.{}.{}.{}.out".format(
            self._scan_params[0],
            getattr(self._input, self._scan_params[0]),
            self._scan_params[1],
            getattr(self._input, self._scan_params[1])
            )
        )
        outfile = lha_utils.THDMCOutput(outname)
        model_point.is_valid_model = outfile.valid_model
        # Get properties of little h boson from output
        h = thdm_scanner.grid.HiggsProperties("h")
        h.mass = outfile.mh
        h.br_tautau = outfile.br_htautau
        # Get properties of heavy H boson from output
        H = thdm_scanner.grid.HiggsProperties("H")
        H.mass = outfile.mH
        H.br_tautau = outfile.br_Htautau
        # Get properties of pseudoscalar A boson from output
        A = thdm_scanner.grid.HiggsProperties("A")
        A.mass = outfile.mA
        A.br_tautau = outfile.br_Atautau
        # Add Yukawa couplings for each Higgs boson based on
        # input parameters
        beta = math.atan(self._input.tanb)  # atan return range: -pi/2(0), pi/2
        alpha = beta - math.acos(self._input.cos_betal)  # acos ret.r: 0, pi
        if self._input.type == 1:
            h.gt = math.cos(alpha) / math.sin(beta)
            h.gb = math.cos(alpha) / math.sin(beta)
            H.gt = math.sin(alpha) / math.sin(beta)
            H.gb = math.sin(alpha) / math.sin(beta)
            A.gt = -1. / math.tan(beta)
            A.gb = 1. / math.tan(beta)
        else:
            h.gt = math.cos(alpha) / math.sin(beta)
            h.gb = -math.sin(alpha) / math.cos(beta)
            H.gt = math.sin(alpha) / math.sin(beta)
            H.gb = math.cos(alpha) / math.cos(beta)
            A.gt = -1. / math.tan(beta)
            A.gb = -math.tan(beta)
        # Add the correct Higgs properties to the model point
        model_point.h = h
        model_point.H = H
        model_point.A = A
        return model_point


class SusHiHarvester(THDMHarvesterABC):

    def __init__(self, outpath="output",
                 scan_parameters=("mH", "tanb"),
                 run_pdf_uncerts=False):
        self._outputfile = os.path.join(outpath, "SusHi.out")
        self._scan_params = scan_parameters
        self._run_uncerts = run_pdf_uncerts

    def set_inputs(self, inputs):
        self._input = inputs

    def harvest_output(self, model_point):
        higgs_dict = {
                11: "h",
                12: "H",
                21: "A"
        }
        for higgs in [11, 12, 21]:
            outfile = lha_utils.SusHiOutput(
                    self._outputfile.replace(".out",
                                             ".{}.{}.{}.{}.H{}.out".format(
                                                 self._scan_params[0],
                                                 getattr(self._input,
                                                         self._scan_params[0]),
                                                 self._scan_params[1],
                                                 getattr(self._input,
                                                         self._scan_params[1]),
                                                 higgs)))
            getattr(model_point, higgs_dict[higgs]).gg_xs = outfile.xs_ggPhi
            getattr(model_point, higgs_dict[higgs]).bb_xs = outfile.xs_bbPhi
            getattr(model_point, higgs_dict[higgs]).bb_xs = outfile.xs_bbPhi
            getattr(model_point, higgs_dict[higgs]).gg_xs_scale_unc = (outfile.xs_ggPhi_scale_down, outfile.xs_ggPhi_scale_up)  # noqa: E501
            # Check mass of considered Higgs boson between SusHi and 2HDMC
            # as cross check.
            if outfile.mPhi != getattr(model_point, higgs_dict[higgs]).mass:
                raise RuntimeError("Mass calculated from SusHi and 2HDMC "
                                   "does not agree. "
                                   "Values are {} and {}".format(
                                       outfile.mPhi,
                                       getattr(model_point,
                                               higgs_dict[higgs]).mass))
            # Read pdf and alpha_s uncertainties if requested
            if self._run_uncerts:
                # Collect all calculated cross sections for pdf and
                # alpha_s variations.
                ggPhi_xsections = []
                bbPhi_xsections = []
                for pdf_member in range(1, 103):
                    outname = self._outputfile.replace(".out",
                                                       ".{}.{}.{}.{}.H{}.pdf{}.out".format(
                                                           self._scan_params[0],
                                                           getattr(self._input,
                                                                   self._scan_params[0]),
                                                           self._scan_params[1],
                                                           getattr(self._input,
                                                                   self._scan_params[1]),
                                                           higgs,
                                                           pdf_member))
                    outfile = lha_utils.SusHiOutput(outname)
                    ggPhi_xsections.append(outfile.xs_ggPhi)
                    bbPhi_xsections.append(outfile.xs_bbPhi)
                # Calculate the uncertainties from the collected values
                # For the pdf uncertainties there are two different
                # possibilieties to calculate the uncertainties (arXiv:)
                # Possibility 1:
                pdf_unc = np.std(ggPhi_xsections[:-2], ddof=1)
                pdf_unc_bbPhi = np.std(bbPhi_xsections[:-2], ddof=1)
                # Possibility 2 (asymmetric non-Gaussian case):
                # sorted_pdf_uncs = sorted(ggPhi_xsections[:-2])
                # pdf_unc = (sorted_pdf_uncs[83]-sorted_pdf_uncs[15]) / 2.
                # sorted_pdf_uncs = sorted(bbPhi_xsections[:-2])
                # pdf_unc_bbPhi = (sorted_pdf_uncs[83]-sorted_pdf_uncs[15]) / 2.

                # Calculate alpha_s uncertainty from remaining variations
                alphas_unc = (ggPhi_xsections[-1] - ggPhi_xsections[-2]) / 2.
                pdf_as_unc = math.sqrt(pdf_unc**2 + alphas_unc**2)
                getattr(model_point, higgs_dict[higgs]).gg_xs_pdfas_unc = (-pdf_as_unc, pdf_as_unc)
                alphas_unc_bbPhi = (bbPhi_xsections[-1] - bbPhi_xsections[-2]) / 2.
                pdf_as_unc_bbPhi = math.sqrt(pdf_unc_bbPhi**2 + alphas_unc_bbPhi**2)
                getattr(model_point, higgs_dict[higgs]).bb_xs_pdfas_unc = (-pdf_as_unc_bbPhi, pdf_as_unc_bbPhi)
        return
