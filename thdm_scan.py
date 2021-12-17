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
    return parser.parse_args()


def setup_logging(level=logging.INFO):
    logger.setLevel(level)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def run_point(mod_pars, model=None, outpath="output"):
    (par_1, val_1), (par_2, val_2) = mod_pars
    model_point = thdm_scanner.THDMPoint((par_1, par_2),
                                         (val_1, val_2))
    # print(model.fixed_model_params)
    dicts = [model.fixed_model_params,
             {par_1: float(val_1), par_2: float(val_2)}]
    input_dict = {k: v for di in dicts for k, v in di.items()}
    print(input_dict)
    inputs = thdm_scanner.THDMInput(**input_dict)
    # Run 2HDMC calculations
    thdm_runner = thdm_scanner.THDMCRunner(outpath=outpath,
                                           scan_parameters=(par_1, par_2))
    thdm_runner.set_inputs(inputs)
    thdm_runner.run()
    thdm_runner.harvest_output(model_point)

    # Run SusHi calculation separately for each Higgs boson
    sushi_runner = thdm_scanner.SusHiRunner(outpath=outpath,
                                            scan_parameters=(par_1, par_2))
    sushi_runner.set_inputs(inputs)
    sushi_runner.run(multiproc=False)
    sushi_runner.harvest_output(model_point)
    return model_point


def main(args):
    with open(args.config, "r") as fi:
        config = yaml.load(fi, Loader=yaml.SafeLoader)

    print(config)
    # Check if output directory exists, create it in case not
    output_path = os.path.basename(args.config).replace(".yaml", "")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    # For each parameter point
    model = thdm_scanner.THDMModel(config["name"],
                                   config["scan_parameter"],
                                   config["benchmark_parameter"])

    if args.procs > 1:
        with multiprocessing.Pool(args.procs) as pool:
            model_points = pool.starmap(run_point,
                                        zip(model.parameter_points(),
                                            repeat(model),
                                            repeat(output_path)))
    else:
        model_points = []
        for mod_pars in model.parameter_points():
            model_points.append(run_point(mod_pars, model))

    for model_point in model_points:
        # Add model point to grid.
        model.add_point(model_point)
    # Write grid of points to root file
    model.write_to_root()
    return


if __name__ == "__main__":
    args = parse_args()
    setup_logging(level=logging.DEBUG if args.verbose else logging.INFO)
    main(args)
