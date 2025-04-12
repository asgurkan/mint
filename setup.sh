#!/bin/bash

set -o errexit   # Exit immediately if a command exits with a non-zero status
set -o nounset   # Treat unset variables as an error when substituting
set -o pipefail  # Pipeline returns the exit status of the last command that fails

conda env create -f mint_setup.yaml  
conda activate mint 
