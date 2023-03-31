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
Finally we need to set up the 2HDMC and SusHi tools to be able to run the workflow. Both tools are available from the websites of the corresponding projects. So, we download the tarballs containing the source code and unpack them to the local directory.
```bash
wget https://2hdmc.hepforge.org/downloads/?f=2HDMC-1.8.0.tar.gz --no-check-certificate -O 2HDMC-1.8.0.tar.gz && tar xzf 2HDMC-1.8.0.tar.gz
wget https://sushi.hepforge.org/downloads/?f=SusHi-1.7.0.tar.gz --no-check-certificate -O SusHi-1.7.0.tar.gz && tar xzf SusHi-1.7.0.tar.gz
```
A few additional changes are necessary to be able to compile the 2HDMC source code. The Makefile does not find the libraries from the LCG stack directory, so we have to add them manually and also add an executable to the list of executables to compile.
```bash
sed -i '/^SOURCES=/ s/$/ runTHDM.cpp/' 2HDMC-1.8.0/Makefile
sed -i '21 c\LDFLAGS+=-L$(LIBDIR) -L/cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc8-opt/lib -l2HDMC -lgsl -lgslcblas -lm' 2HDMC-1.8.0/Makefile
```
Additionally we set the values for SM input parameters for which the recommendation in YR4 differs from the default values in the exectuables we want to use:
```bash
sed -i '/SM sm;/ c\  SM sm = model.get_SM();\n  sm.set_alpha_s(0.118);\n  sm.set_MZ(91.1876);\n  sm.set_qmass_msbar(5, 4.18);\n  sm.set_qmass_pole(6, 172.5);\n  \n  model.set_SM(sm);' 2HDMC-1.8.0/src/CalcHybrid.cpp
sed -i '/SM sm;/ c\  SM sm = model.get_SM();\n  sm.set_alpha_s(0.118);\n  sm.set_MZ(91.1876);\n  sm.set_qmass_msbar(5, 4.18);\n  sm.set_qmass_pole(6, 172.5);\n  \n  model.set_SM(sm);' 2HDMC-1.8.0/src/CalcPhys.cpp
```
Afterwards we can compile the executables for 2HDMC.
```bash
pushd 2HDMC-1.8.0
make
popd
```
Then we set up SusHi by specifying the path to the 2HDMC directory and the correct 2HDMC version number. Additionally, we have to provide the location of the LHA library and can compile SusHi afterwards.
```bash
sed -i "23 c\2HDMCPATH = ${PWD}/2HDMC-1.8.0" SusHi-1.7.0/Makefile
sed -i "24 c\2HDMCVERSION = 1.8.0" SusHi-1.7.0/Makefile
sed -i "35 c\LHAPATH = /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc8-opt/lib" SusHi-1.7.0/Makefile
# Start compiling after the edits
pushd SusHi-1.7.0
./configure
make predef=2HDMC
popd
```

## Running the scans
Before running the workflow with law add the path to the installation directory of 2HDMC and SusHi line 20 in the setup script for law in the `law_config/setup.sh` script and source the script afterwards. The workflow can then be run with:

```bash
law run CollectScan --version v1 --scenario $scenario --run-pdfas-uncerts
```
where `$scenario` is one of the scenarios provided in `law_config/models` without the file extension.
