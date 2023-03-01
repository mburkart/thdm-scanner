# thdm-scanner

Code and workflow to perform scans of 2HDM parameter spaces.
The theoretical calculations are performed using the publicly available theory codes [2HDMC](https://2hdmc.hepforge.org/) and [SusHi](https://sushi.hepforge.org/) in the versions 1.8.0 and 1.7.0. The workflow creates configurations for these codes, runs them and collects the inputs into a single root files containing 2D histograms as function of two model parameters. The workflow is implemented within the law framework.

## Checkout instructions
We start by cloning this repository
```bash
git clone https://github.com/mburkart/thdm-scanner.git
cd thdm-scanner/
```
Afterwards we need to get and install the dependencies of the package to be able to run it. To provide some of the dependencies we rely on LCG stacks, thus we start by setting up the LCG stack.
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc8-opt/setup.sh
```
Afterwards we set up a virtual environment and install law in this environment via pip:
```bash
python3 -m venv law_env
source law_env/bin/activate
python3 -m pip install law
```
Finally we need to set up the 2HDMC and SusHi tools to be able to run the workflow. Both tools are available from the websites of the corresponding projects.
```bash
wget https://2hdmc.hepforge.org/downloads/?f=2HDMC-1.8.0.tar.gz --no-check-certificate
```
