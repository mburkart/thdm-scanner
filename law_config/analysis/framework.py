# coding: utf-8

"""
Law example tasks to demonstrate HTCondor workflows at CERN.

In this file, some really basic tasks are defined that can be inherited by
other tasks to receive the same features. This is usually called "framework"
and only needs to be defined once per user / group / etc.
"""


import os
import math

import luigi
import law


# the htcondor workflow implementation is part of a law contrib package
# so we need to explicitly load it
law.contrib.load("htcondor")


class Task(law.Task):
    """
    Base task that we use to force a version parameter on all inheriting tasks, and that provides
    some convenience methods to create local file and directory targets at the default data path.
    """

    version = luigi.Parameter()
    scenario = luigi.ChoiceParameter(str,
                                     choices=[
                                         "BP1_Type1_mod",
                                         "BP1_Type2",
                                         "HWWLike_Type1",
                                         "HWWLike_Type2",
                                         "PASComp_Type1",
                                         "PASComp_Type2",
                                         ],
                                     description="Name of the scenario that is run."
            "This name is also used as the identifier for the configuration file")
    run_pdfas_uncerts = luigi.BoolParameter(description="Run extra variations for pdf and alpha_s uncertainties.")

    def store_parts(self):
        return (self.__class__.__name__, self.scenario, self.version)

    def local_path(self, *path):
        # ANALYSIS_DATA_PATH is defined in setup.sh
        parts = (os.getenv("ANALYSIS_DATA_PATH"),) + self.store_parts() + path
        return os.path.join(*parts)

    def local_target(self, *path):
        return law.LocalFileTarget(self.local_path(*path))

    # def remote_path(self, *path):
    #     # ANALYSIS_DATA_PATH is defined in setup.sh
    #     parts = (os.getenv("ANALYSIS_DATA_PATH"),) + self.store_parts() + path
    #     return os.path.join(*parts)

    # def remote_target(self, *path):
    #     return law.RemoteFileTarget(self.local_path(*path))


class HTCondorWorkflow(law.htcondor.HTCondorWorkflow):
    """
    Batch systems are typically very heterogeneous by design, and so is HTCondor. Law does not aim
    to "magically" adapt to all possible HTCondor setups which would certainly end in a mess.
    Therefore we have to configure the base HTCondor workflow in law.contrib.htcondor to work with
    the CERN HTCondor environment. In most cases, like in this example, only a minimal amount of
    configuration is required.
    """

    htcondor_accounting_group = luigi.Parameter()
    htcondor_requirements = luigi.Parameter()
    htcondor_remote_job = luigi.Parameter()
    htcondor_walltime = luigi.Parameter()
    htcondor_request_cpus = luigi.Parameter()
    htcondor_request_memory = luigi.Parameter()
    htcondor_universe = luigi.Parameter()
    htcondor_docker_image = luigi.Parameter()
    htcondor_request_disk = luigi.Parameter()

    def htcondor_output_directory(self):
        # the directory where submission meta data should be stored
        return law.LocalDirectoryTarget(self.local_path())

    def htcondor_bootstrap_file(self):
        # each job can define a bootstrap file that is executed prior to the actual job
        # in order to setup software and environment variables
        return law.util.rel_path(__file__, "bootstrap.sh")

    def htcondor_job_config(self, config, job_num, branches):
        # render_variables are rendered into all files sent with a job
        config.render_variables["analysis_path"] = os.getenv("ANALYSIS_PATH")

        # force to run on CC7, http://batchdocs.web.cern.ch/batchdocs/local/submit.html#os-choice
        config.custom_content.append(("Requirements", self.htcondor_requirements))
        config.custom_content.append(("+RemoteJob", self.htcondor_remote_job))
        config.custom_content.append(("accounting_group", self.htcondor_accounting_group))
        config.custom_content.append(("universe", self.htcondor_universe))
        config.custom_content.append(("docker_image", self.htcondor_docker_image))
        config.custom_content.append(("+RequestWalltime", self.htcondor_walltime))
        config.custom_content.append(("request_cpus", self.htcondor_request_cpus))
        config.custom_content.append(("RequestMemory", self.htcondor_request_memory))
        config.custom_content.append(("RequestDisk", self.htcondor_request_disk))

        # copy the entire environment
        config.custom_content.append(("getenv", "true"))

        # the CERN htcondor setup requires a "log" config, but we can safely set it to /dev/null
        # if you are interested in the logs of the batch system itself, set a meaningful value here
        config.custom_content.append(("log", "/dev/null"))

        return config
