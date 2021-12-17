#!/usr/bin/env python

from abc import ABCMeta, abstractmethod
import os
import subprocess
import multiprocessing
import logging

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

    @abstractmethod
    def harvest_output(self, model_point):
        pass


class THDMCRunner(THDMRunnerABC):
    """Run THDMC with inputs in Hybrid basis"""

    def __init__(self, outpath="output",
                 scan_parameters=("mH", "tanb")):
        self._outputfile = os.path.join(outpath, "2HDMC_output.out")
        self._scan_params = scan_parameters

    def set_inputs(self, inputs):
        self._input = inputs

    def run(self):
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
        logger.debug("Check if output file has been successfully written.")
        if not os.path.exists(self._outputfile):
            raise RuntimeError("Output file of 2HDMC run not created.")
        else:
            logger.debug("Output file {} has been successfully written."
                         .format(self._outputfile))
        return

    def harvest_output(self, model_point):
        outfile = lha_utils.THDMCOutput(self._outputfile)
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
        # Add the correct Higgs properties to the model point
        model_point.h = h
        model_point.H = H
        model_point.A = A
        return model_point


class SusHiRunner(THDMRunnerABC):

    def __init__(self, outpath="output",
                 scan_parameters=("mH", "tanb")):
        self._inputfile = os.path.join(outpath, "SusHi_input.in")
        self._outputfile = os.path.join(outpath, "SusHi_output.out")
        self._scan_params = scan_parameters

    def set_inputs(self, inputs):
        self._input = inputs

    def run(self, multiproc=True):
        # For each Higgs boson
        if multiproc:
            with multiprocessing.Pool(3) as pool:
                pool.map(self._run_single_higgs, [11, 12, 21])
        else:
            for higgs in [11, 12, 21]:
                self._run_single_higgs(higgs)
        return

    def _run_single_higgs(self, higgs):
        # Prepare input file for Higgs boson
        infile_name = self._inputfile.replace(".in",
                                              "{}.{}.{}.{}.H{}.in".format(
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
        infile.cos_betal = self._input.cos_betal
        infile.Z4 = self._input.Z4
        infile.Z5 = self._input.Z5
        infile.Z7 = self._input.Z7
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
                         infile_name.replace("in", "out")])  # noqa: E501
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
                                             "{}.{}.{}.{}.H{}.out".format(
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
        return
