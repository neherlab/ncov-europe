#!/bin/env bash


PATH=$PATH:$HOME/../GROUP/bin/aws

# activate augur conda env
source $HOME/../neher/miniconda3/etc/profile.d/conda.sh
conda activate nextstrain

nohup snakemake $@ > snakemake_log &


