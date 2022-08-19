# coding: utf-8

"""
Law example tasks to demonstrate HTCondor workflows at CERN.

The actual payload of the tasks is rather trivial.
"""

import yaml
import os
import shutil
import logging
logging.basicConfig(level=logging.INFO)
from itertools import repeat
from multiprocessing import Pool

import law
import luigi
import thdm_scanner


# import our "framework" tasks
from analysis.framework import Task, HTCondorWorkflow


def load_tarballs(infile, outpath):
    infile.load(outpath)


def create_inputs(mod_pars, model,
                  hybrid_basis=True):
    # Create a model point in the 2D model plane
    (par_1, val_1), (par_2, val_2) = mod_pars
    # Build the dictionary of the inputs from the
    # fixed parameters of the model and the scanned values
    # to initialize the input object.
    dicts = [model.fixed_model_params,
             {par_1: float(val_1), par_2: float(val_2)}]
    input_dict = {k: v for di in dicts for k, v in di.items()}
    if hybrid_basis:
        inputs = thdm_scanner.THDMInput(**input_dict)
    else:
        inputs = thdm_scanner.THDMPhysicsInput(**input_dict)
    return inputs


class PerformScan(Task, HTCondorWorkflow, law.LocalWorkflow):
    """Some docstring

    Some more information.
    """

    def htcondor_job_config(self, config, job_num, branches):
        config = super().htcondor_job_config(config, job_num, branches)
        config.custom_content.append(("JobBatchName", "-".join(self.store_parts())))
        return config


    def create_branch_map(self):
        with open(os.getenv("ANALYSIS_PATH")+"/models/{}.yaml".format(self.scenario),
                  "r") as fi:
            config = yaml.safe_load(fi)
        model = thdm_scanner.THDMModel(
                    config["name"],
                    config["scan_parameter"],
                    config["benchmark_parameter"],
                    include_pdfas_unc=self.run_pdfas_uncerts
                    )
        # map branch indexes to ascii numbers from 97 to 122 ("a" to "z")
        return {i: create_inputs(mod_pars, model, hybrid_basis=not "physical_basis" in config.keys())
                for i, mod_pars in enumerate(model.parameter_points())}

    def output(self):
        # it's best practice to encode the branch number into the output target
        return self.local_target("output_{}.tar.gz".format(self.branch))

    def run_point(self, inputs, outpath="output",
                  hybrid_basis=True, run_pdf_uncerts=False):
        # Run 2HDMC calculations
        with open(os.getenv("ANALYSIS_PATH")+"/models/{}.yaml".format(self.scenario),
                  "r") as fi:
            par_1, par_2 = yaml.safe_load(fi)["scan_parameter"].keys()
        thdm_runner = thdm_scanner.THDMCRunner(outpath=outpath,
                                               scan_parameters=(par_1, par_2))
        thdm_runner.set_inputs(inputs)
        thdm_runner.run(hybrid_basis)

        # Run SusHi calculation separately for each Higgs boson
        sushi_runner = thdm_scanner.SusHiRunner(
                outpath=outpath,
                scan_parameters=(par_1, par_2),
                run_pdf_uncerts=self.run_pdfas_uncerts
                )
        sushi_runner.set_inputs(inputs)
        sushi_runner.run(multiproc=False, hybrid_basis=hybrid_basis)
        return

    def run(self):
        # the branch data holds the inputs for the specific grid point
        inp = self.branch_data

        if self.workflow == "local":
            os.chdir(os.getenv("ANALYSIS_PATH"))
            outpath = "output"
        else:
            outpath = os.getcwd()
        # actual payload: run 2HDMC and SusHi to get the values to collect
        self.run_point(inp, outpath=outpath,
                       hybrid_basis=not isinstance(inp, thdm_scanner.THDMPhysicsInput),
                       run_pdf_uncerts=self.run_pdfas_uncerts)

        # use target formatters (implementing dump and load,
        # based on the file extension)
        # to write the output target
        output = self.output()
        with open(os.getenv("ANALYSIS_PATH")+"/models/{}.yaml".format(self.scenario),
                  "r") as fi:
            par_1, par_2 = yaml.safe_load(fi)["scan_parameter"].keys()
        point_str = "{}.{}.{}.{}".format(par_1,
                                         getattr(inp, par_1),
                                         par_2,
                                         getattr(inp, par_2))
        file_list = ["2HDMC_output.{}.out".format(point_str),
                     "SusHi.{}.H11.out".format(point_str),
                     "SusHi.{}.H12.out".format(point_str),
                     "SusHi.{}.H21.out".format(point_str),
                    ]
        if self.run_pdfas_uncerts:
            for h in ["11", "12", "21"]:
                file_list.extend(["SusHi.{}.H{}.pdf{}.out".format(point_str, h, i) for i in range(1, 103)])
        # Check if all expected files have been written
        for fi in file_list:
            if not os.path.exists(os.path.join(outpath, fi)):
                raise Exception("Not all outputs have been written...")
        output.dump(map(os.path.join, repeat(outpath), file_list))


class CollectScan(Task):
    """
    Some docstring.
    """

    skip_unpacking = luigi.BoolParameter(description="Skip unpacking of tarballs for crashed collection")

    def requires(self):
        # req() is defined on all tasks and handles the passing of all
        # parameter values that are common between the required task
        # and the instance (self)
        # note that the workflow is required (branch -1, the default),
        # not the particular branch
        return PerformScan.req(self)

    def output(self):
        # output a plain text file
        return self.local_target(self.scenario + ".root")

    def run(self):
        with open(os.getenv("ANALYSIS_PATH")+"/models/{}.yaml".format(self.scenario),
                  "r") as fi:
            config = yaml.safe_load(fi)
        model = thdm_scanner.THDMModel(
                config["name"],
                config["scan_parameter"],
                config["benchmark_parameter"],
                include_pdfas_unc=self.run_pdfas_uncerts
                )
        print("Status: Unzipping tarball results...")
        # Running with pool not possible as class/instance methods are used
        if self.skip_unpacking:
            pass
        else:
            with Pool(2) as pool:
                pool.starmap(load_tarballs,
                             zip(self.input()["collection"].targets.values(),
                                 repeat("tarballs_results/{}".format(self.scenario)))
                             )
        # for inp in self.input()["collection"].targets.values():
        #     inp.load("tarballs_results/{}".format(self.scenario))
        print("Status: Parsing unzipped results...")
        for mod_pars in model.parameter_points():
            (par_1, val_1), (par_2, val_2) = mod_pars
            inputs = create_inputs(mod_pars, model,
                                   hybrid_basis="physical_basis" not in config)
            # Collect results from written output files
            model_point = thdm_scanner.THDMPoint((par_1, par_2),
                                                 (val_1, val_2))
            model_point.cos_betal = inputs.cos_betal
            # First collect results from 2HDMC calculations
            thdm_harvester = thdm_scanner.THDMCHarvester(
                    # outpath=os.getenv("ANALYSIS_DATA_PATH"),
                    outpath=os.path.join(os.getcwd(), "tarballs_results/{}".format(self.scenario)),
                    scan_parameters=(par_1, par_2)
                    )
            thdm_harvester.set_inputs(inputs)
            thdm_harvester.harvest_output(model_point)
            # Collect results from SusHi calculations
            sushi_harvester = thdm_scanner.SusHiHarvester(
                    # outpath=os.getenv("ANALYSIS_DATA_PATH"),
                    outpath=os.path.join(os.getcwd(), "tarballs_results/{}".format(self.scenario)),
                    scan_parameters=(par_1, par_2),
                    run_pdf_uncerts=self.run_pdfas_uncerts)
            sushi_harvester.set_inputs(inputs)
            sushi_harvester.harvest_output(model_point)

            # Add the harvested point to the model
            model.add_point(model_point)
        model.write_to_root()
        self.publish_message("\nbuilt THDM scan for model: {}\n".format(self.scenario))
        output = self.output()
        output.move_from_local(self.scenario + ".root")
        # Clean up directory with unpacked tar files
        print("Status: Cleaning up unzipped results...")
        shutil.rmtree("tarballs_results/{}".format(self.scenario))
