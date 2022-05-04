# coding: utf-8

"""
Law example tasks to demonstrate HTCondor workflows at CERN.

The actual payload of the tasks is rather trivial.
"""

import yaml

import six
import law
import thdm_scanner


# import our "framework" tasks
from analysis.framework import Task, HTCondorWorkflow


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
    print(input_dict)
    if hybrid_basis:
        inputs = thdm_scanner.THDMInput(**input_dict)
    else:
        inputs = thdm_scanner.THDMPhysicsInput(**input_dict)
    return inputs


class PerformScan(Task, HTCondorWorkflow, law.LocalWorkflow):
    """Some docstring

    Some more information.
    """

    def create_branch_map(self):
        with open("???/{}.yaml".format(self.scenario), "r") as fi:
            config = yaml.safe_load(fi)
        model = thdm_scanner.THDMModel(config["name"],
                                       config["scan_parameter"],
                                       config["benchmark_parameter"],
                                       include_pdfas_unc=self.run_pdfas_uncerts)
        # map branch indexes to ascii numbers from 97 to 122 ("a" to "z")
        return {i: create_inputs(mod_pars, model) for i, mod_pars in enumerate(model.parameter_points())}

    def output(self):
        # it's best practice to encode the branch number into the output target
        return self.local_target("output_{}.txt".format(self.branch))

    def run_point(inputs, outpath="output",
                  hybrid_basis=True, run_pdf_uncerts=False):
        # Run 2HDMC calculations
        thdm_runner = thdm_scanner.THDMCRunner(outpath=outpath,
                                               scan_parameters=(par_1, par_2))
        thdm_runner.set_inputs(inputs)
        thdm_runner.run(hybrid_basis)

        # Run SusHi calculation separately for each Higgs boson
        sushi_runner = thdm_scanner.SusHiRunner(outpath=outpath,
                                                scan_parameters=(par_1, par_2),
                                                run_pdf_uncerts=run_pdf_uncerts)
        sushi_runner.set_inputs(inputs)
        sushi_runner.run(multiproc=False, hybrid_basis=hybrid_basis)
        return

    def run(self):
        # the branch data holds the integer number to convert
        inp = self.branch_data

        # actual payload: convert to char
        run_point(inp, not isinstance(inp, THDMPhysicsInput),
                  run_pdf_uncerts=self.run_pdfas_uncerts)

        # use target formatters (implementing dump and load, based on the file extension)
        # to write the output target
        output = self.output()
        output.dump({"num": num, "char": char})


class CollectScan(Task):
    """
    Some docstring.
    """

    def requires(self):
        # req() is defined on all tasks and handles the passing of all parameter values that are
        # common between the required task and the instance (self)
        # note that the workflow is required (branch -1, the default), not the particular branch
        # tasks (branches [0, 26))
        return PerformScan.req(self)

    def output(self):
        # output a plain text file
        return self.local_target("???.root")

    def run(self):
        with open("???/{}.yaml".format(self.scenario), "r") as fi:
            config = yaml.safe_load(fi)
        model = thdm_scanner.THDMModel(config["name"],
                                       config["scan_parameter"],
                                       config["benchmark_parameter"],
                                       include_pdfas_unc=self.run_pdfas_uncerts)
        for mod_pars in model.parameter_points():
            (par_1, val_1), (par_2, val_2) = mod_pars
            inputs = create_inputs(mod_pars, model, hybrid_basis="physical_basis" not in config)
            # Collect results from written output files
            model_point = thdm_scanner.THDMPoint((par_1, par_2),
                                                 (val_1, val_2))
            model_point.cos_betal = inputs.cos_betal
            # First collect results from 2HDMC calculations
            thdm_harvester = THDMCHarvester(outpath=outpath,
                                            scan_parameters=(par_1, par_2))
            thdm_harvester.set_inputs(inputs)
            thdm_harvester.harvest_output(model_point)
            # Collect results from SusHi calculations
            sushi_harvester = SusHiHarvester(outpath=outpath,
                                             scan_parameters=(par_1, par_2),
                                             run_pdf_uncerts=self.run_pdfas_uncerts)
            sushi_harvester.set_inputs(inputs)
            sushi_harvester.harvest_output(model_point)

            # Add the harvested point to the model
            model.add_point(model_point)
        model.write_to_root()
        self.publish_message("\nbuilt THDM scan for model: {}\n".format(self.scenario))
