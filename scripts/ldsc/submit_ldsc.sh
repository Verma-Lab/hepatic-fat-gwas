#!/usr/bin/env bash
#BSUB -q epistasis_normal
#BSUB -J ldsc_heritability
#BSUB -o ldsc_heritability.%J.out
#BSUB -e ldsc_heritability.%J.err
#BSUB -n 1
#BSUB -W 6:00
#BSUB -M 16000
#BSUB -R "rusage[mem=16000]"

set -euo pipefail
cd /home/agaro/verma_shared/projects/Liver_IDPs/scripts/ldsc
./run_ldsc_heritability.sh