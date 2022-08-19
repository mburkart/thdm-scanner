#!/usr/bin/env bash

action() {
    # determine the directory of this file
    local this_file="$( [ ! -z "$ZSH_VERSION" ] && echo "${(%):-%x}" || echo "${BASH_SOURCE[0]}" )"
    local this_dir="$( cd "$( dirname "$this_file" )" && pwd )"

    export PYTHONPATH="$this_dir:$PYTHONPATH"
    export LAW_HOME="$this_dir/.law"
    export LAW_CONFIG_FILE="$this_dir/law.cfg"
    export LUIGI_CONFIG_PATH="$this_dir/luigi.cfg"

    export ANALYSIS_PATH="$this_dir"
    export ANALYSIS_DATA_PATH="$ANALYSIS_PATH/data"

    source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc8-opt/setup.sh
    export LHAPATH=/cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc8-opt/lib
    export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/:/cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc8-opt/share/LHAPDF/
    export PYTHONPATH="${this_dir}/../:$PYTHONPATH"
    export THEORY_CODE_PATH=/work/mburkart/Run2Legacy/2HDM_Theory/
    source "$( law completion )" ""
}
action
