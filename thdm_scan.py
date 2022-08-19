#!/usr/bin/env python

import yaml
import argparse
import logging
import os
from itertools import repeat
import multiprocessing

import thdm_scanner


logger = logging.getLogger("")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "-c", "--config",
            type=str,
            required=True,
            help="Input config file"
    )
    parser.add_argument(
            "-o", "--output-path",
            type=str,
            default=None,
            help="Optional path to output directory. Will be created if not present."
    )
    parser.add_argument(
            "-v", "--verbose",
            action="store_true",
            help="Enable debug logging"
    )
    parser.add_argument(
            "-p", "--processes",
            dest="procs",
            type=int,
            default=1,
            help="Number of processes used to "
            "compute the grid of points."
    )
    parser.add_argument(
            "-u", "--run-pdfas-uncerts",
            action="store_true",
            help="Run variations for pdf and alpha_s uncertainties"
    )
    return parser.parse_args()


def setup_logging(level=logging.INFO):
    logger.setLevel(level)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)


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


def collect_result(inputs, outpath="output",
                   run_pdf_uncerts=False):
    (par_1, val_1), (par_2, val_2) = mod_pars
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
                                     run_pdf_uncerts=run_pdf_uncerts)
    sushi_harvester.set_inputs(inputs)
    sushi_harvester.harvest_output(model_point)
    return model_point


def main(args):
    with open(args.config, "r") as fi:
        config = yaml.load(fi, Loader=yaml.SafeLoader)

    print(config)
    if args.output_path is None:
        # By default write into the current directory
        output_path = os.path.basename(args.config).replace(".yaml", "")
    else:
        output_path = os.path.join(args.output_path,
                                   os.path.basename(args.config).replace(".yaml", "")
        )
    # Check if output directory exists, create it in case not
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    # For each parameter point
    model = thdm_scanner.THDMModel(config["name"],
                                   config["scan_parameter"],
                                   config["benchmark_parameter"],
                                   include_pdfas_unc=args.run_pdfas_uncerts)

    if "physical_basis" not in config:
        hybrid_basis = True
    else:
        hybrid_basis = False
    if args.procs > 1:
        with multiprocessing.Pool(args.procs) as pool:
            pool.starmap(run_point,
                         zip(map(create_inputs, model.parameter_points(), repeat(model)),
                         repeat(output_path),
                         repeat(hybrid_basis),
                         repeat(args.run_pdfas_uncerts)))
    else:
        for mod_pars in model.parameter_points():
            inputs = create_inputs(mod_pars, model, hybrid_basis=hybrid_basis)
            run_point(inputs,
                      output_path, hybrid_basis,
                      args.run_pdfas_uncerts))

    for mod_pars in model.parameter_points():
        # Add model point to grid.
        model.add_point(collect_result(mod_pars, model,
                                       output_path, hybrid_basis,
                                       args.run_pdfas_uncerts))
    # Write grid of points to root file
    model.write_to_root()
    return


if __name__ == "__main__":
    args = parse_args()
    setup_logging(level=logging.DEBUG if args.verbose else logging.INFO)
    main(args)
