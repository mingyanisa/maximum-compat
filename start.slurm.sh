#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
for i in {1..20}; do
    sbatch --export ALL -N 2 -p standard --mem 4096 -t 01:00:00 "${DIR}/generate_tree.sh" $DIR $i
done