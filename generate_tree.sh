#!/usr/bin/env bash
input_dir=${1:-"$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"}
idx=${2:-'0'}
source /hps/software/opt/linux-rocky8-cascadelake/anaconda3/etc/profile.d/conda.sh
conda activate compat
compat ${input_dir}/clusters/cluster.$idx.mfa ${input_dir}/results/cluster.$idx.netwick