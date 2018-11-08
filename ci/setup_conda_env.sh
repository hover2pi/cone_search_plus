#!/bin/bash

echo "Creating a Python $PYTHON_VERSION environment"
conda create -n csp python=$PYTHON_VERSION || exit 1
source activate csp

echo "Installing packages..."
conda install flake8 beautifulsoup4 lxml numpy astropy matplotlib
pip install pytest pytest-cov coveralls