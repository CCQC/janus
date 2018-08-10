[![Build Status](https://travis-ci.com/bzhang25/janus.svg?token=6xT7vBmfnsKxZPabRnWW&branch=master)](https://travis-ci.com/bzhang25/janus)
[![codecov](https://codecov.io/gh/bzhang25/janus/branch/master/graph/badge.svg?token=oncB2345LQ)](https://codecov.io/gh/bzhang25/janus)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# Janus
A Python library for adaptive QM/MM methods 

## Set up steps: 
* Download miniconda [installer](https://conda.io/docs/user-guide/install/macos.html)
* Go to folder where this is downloaded and run: `bash Miniconda3-latest-OSX-x86_64.sh -b -p $HOME/miniconda`
* To set conda path to miniconda: `echo PATH="\$HOME/miniconda/bin:\$PATH" >> ~/.bash_profile`
* Create environment in miniconda: `conda create -n janus python=3.6 psi4 psi4-rt numpy openmm -c psi4 -c omnia`
* To activate the environment: `source activate janus`
* To install locally: `pip install -e .`
* For reading external datafiles with pytest within the janus environment: `pip install pytest-datafiles` 
* Adding mendeleev package: `pip install mendeleev` 
* To add MDtraj: conda install -c conda-forge mdtraj 
