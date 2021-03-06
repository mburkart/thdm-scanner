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


class THDMRunnerABC(object):
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def set_inputs(self):
        pass

    @abstractmethod
    def run(self):
        pass


class THDMCRunner(THDMRunnerABC):
    """Run THDMC with inputs in Hybrid basis"""

    def __init__(self, outpath="output",
                 scan_parameters=("mH", "tanb")):
        self._outputfile = os.path.join(outpath, "2HDMC_output.out")
        self._scan_params = scan_parameters

    def set_inputs(self, inputs):
        self._input = inputs

    def run(self, hybrid_basis=True):
        self._outputfile = self._outputfile.replace(
                ".out",
                ".{}.{}.{}.{}.out".format(self._scan_params[0],
                                          getattr(self._input,
                                                  self._scan_params[0]),
                                          self._scan_params[1],
                                          getattr(self._input,
                                                  self._scan_params[1])
                                          ))
        # Run 2HDMC code with correct inputs
        if hybrid_basis:
            logger.info("Running 2HDMC with command {}".format(
                            " ".join([os.path.join(os.environ["THEORY_CODE_PATH"],
                                                   "2HDMC-1.8.0",
                                                   "CalcHybrid"),
                                      str(self._input.mh),
                                      str(self._input.mH),
                                      str(self._input.cos_betal),
                                      str(self._input.Z4),
                                      str(self._input.Z5),
                                      str(self._input.Z7),
                                      str(self._input.tanb),
                                      str(self._input.type),
                                      self._outputfile])))
            subprocess.run([os.path.join(os.environ["THEORY_CODE_PATH"],
                                         "2HDMC-1.8.0",
                                         "CalcHybrid"),
                            str(self._input.mh),
                            str(self._input.mH),
                            str(self._input.cos_betal),
                            str(self._input.Z4),
                            str(self._input.Z5),
                            str(self._input.Z7),
                            str(self._input.tanb),
                            str(self._input.type),
                            self._outputfile],
                           stdout=subprocess.DEVNULL)
        else:
            logger.info("Running 2HDMC with command {}".format(
                            " ".join([os.path.join(os.environ["THEORY_CODE_PATH"],
                                                   "2HDMC-1.8.0",
                                                   "CalcPhys"),
                                      str(self._input.mh),
                                      str(self._input.mH),
                                      str(self._input.mA),
                                      str(self._input.mHp),
                                      str(self._input.sin_betal),  # Conversion between conventions done in input class
                                      str(self._input.lambda6),
                                      str(self._input.lambda7),
                                      str(self._input.m12_square),
                                      str(self._input.tanb),
                                      str(self._input.type),
                                      self._outputfile])))
            subprocess.run([os.path.join(os.environ["THEORY_CODE_PATH"],
                                         "2HDMC-1.8.0",
                                         "CalcPhys"),
                            str(self._input.mh),
                            str(self._input.mH),
                            str(self._input.mA),
                            str(self._input.mHp),
                            str(self._input.sin_betal),  # Conversion between conventions done in input class
                            str(self._input.lambda6),
                            str(self._input.lambda7),
                            str(self._input.m12_square),
                            str(self._input.tanb),
                            str(self._input.type),
                            self._outputfile],
                           stdout=subprocess.DEVNULL)
        logger.debug("Check if output file has been successfully written.")
        if not os.path.exists(self._outputfile):
            raise RuntimeError("Output file of 2HDMC run not created.")
        else:
            logger.debug("Output file {} has been successfully written."
                         .format(self._outputfile))
        return


class SusHiRunner(THDMRunnerABC):

    def __init__(self, outpath="output",
                 scan_parameters=("mH", "tanb"),
                 run_pdf_uncerts=False):
        self._inputfile = os.path.join(outpath, "SusHi_input.in")
        self._outputfile = os.path.join(outpath, "SusHi.out")
        self._scan_params = scan_parameters
        self._run_uncerts = run_pdf_uncerts

    def set_inputs(self, inputs):
        self._input = inputs

    def run(self, multiproc=True,
            hybrid_basis=True):
        # For each Higgs boson
        if multiproc:
            with multiprocessing.Pool(3) as pool:
                pool.map(self._run_single_higgs, [11, 12, 21],
                         hybrid_basis)
        else:
            for higgs in [11, 12, 21]:
                self._run_single_higgs(higgs, hybrid_basis)
        return

    def _run_single_higgs(self, higgs,
                          hybrid_basis=True):
        # Prepare input file for Higgs boson
        infile_name = self._inputfile.replace(".in",
                                              ".{}.{}.{}.{}.H{}.in".format(
                                                  self._scan_params[0],
                                                  getattr(self._input,
                                                          self._scan_params[0]),  # noqa: E501
                                                  self._scan_params[1],
                                                  getattr(self._input,
                                                          self._scan_params[1]),  # noqa: E501
                                                  higgs))
        infile = lha_utils.SusHiInput(infile_name)
        infile.higgs_boson = higgs
        infile.mh = self._input.mh
        infile.mH = self._input.mH
        if hybrid_basis:
            infile.sin_betal = self._input.sin_betal
            infile.Z4 = self._input.Z4
            infile.Z5 = self._input.Z5
            infile.Z7 = self._input.Z7
        else:
            # Conversion to correct angle convention done in input class
            infile.sin_betal = self._input.sin_betal
            infile.mA = self._input.mA
            infile.mHp = self._input.mHp
            infile.lambda6 = self._input.lambda6
            infile.lambda7 = self._input.lambda7
            infile.m12 = math.sqrt(self._input.m12_square)
            infile.thdm_basis = 2
        infile.tanb = self._input.tanb
        infile.thdm_type = self._input.type
        # Write prepared input to input file
        infile.write_file()
        # Run Sushi with prepared input
        subprocess.run([os.path.join(os.environ["THEORY_CODE_PATH"],
                                      "SusHi-1.7.0",
                                      "bin",
                                      "sushi"),
                         infile_name,
                         infile_name.replace("in", "out").replace("_output", "")])  # noqa: E501
        # Perform a separate run for each PDF and alpha_s
        # replica to calculate uncertainty
        if self._run_uncerts:
            for pdf_member in range(1, 103):
                inname_update = infile_name.replace(
                    ".in",
                    ".pdf{}.in".format(pdf_member))
                # Write prepared input to input file
                infile.pdf_set = pdf_member
                # TODO: alpha_s needs to be set correctly for the two alpha_s variations
                if pdf_member == 101:
                    infile.alpha_s = 0.1165
                elif pdf_member == 102:
                    infile.alpha_s = 0.1195
                infile.write_file(inname_update)
                # Run Sushi with prepared input
                subprocess.run([os.path.join(os.environ["THEORY_CODE_PATH"],
                                              "SusHi-1.7.0",
                                              "bin",
                                              "sushi"),
                                 inname_update,
                                 inname_update.replace("in", "out").replace("_output", "")])  # noqa: E501
        return

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
                                   "do not agree. "
                                   "Values are {} and {}".format(
                                       outfile.mPhi,
                                       getattr(model_point,
                                               higgs_dict[higgs]).mass))
            # Read pdf and alpha_s uncertainties if requested
            if self._run_uncerts:
                # Collect all calculated cross sections for pdf and
                # alpha_s variations.
                ggPhi_xsections = []
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
                    # Remove all files related to the calculation
                    # of pdf and alpha_s uncertainties
                    os.remove(outname)
                    os.remove(outname.replace(".out", ".in").replace("SusHi", "SusHi_input"))
                    os.remove(outname.replace(".out", ".out_murdep"))
                # Calculate the uncertainties from the collected values
                # For the pdf uncertainties there are two different
                # possibilieties to calculate the uncertainties (arXiv:)
                # Possibility 1:
                pdf_unc = np.std(ggPhi_xsections[:-2], ddof=1)
                # Possibility 2 (asymmetric non-Gaussian case):
                # sorted_pdf_uncs = sorted(ggPhi_xsections[:-2])
                # pdf_unc = (sorted_pdf_uncs[83]-sorted_pdf_uncs[15]) / 2.

                # Calculate alpha_s uncertainty from remaining variations
                alphas_unc = (ggPhi_xsections[-1] - ggPhi_xsections[-2]) / 2.
                pdf_as_unc = math.sqrt(pdf_unc**2 + alphas_unc**2)
                getattr(model_point, higgs_dict[higgs]).gg_xs_pdfas_unc = (-pdf_as_unc, pdf_as_unc)
        return
