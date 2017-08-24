#!/bin/bash

#
# Usage: ./submit_quant.sh
#
set -euo pipefail 
set -x

# Google project 
GS_PROJECT=${GS_PROJECT:-XXXXXXX}

# Docker image to use
DOCKER_IMAGE=nareshr/kallisto:v0.43

# Reference database
REF=Homo_sapiens.GRCh37.cdna.all.fa.gz
REFURL=ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/cdna/${REF}

# Storage locations
GS_BUCKET=${GS_BUCKET:-cgc-kallisto-example}
GS_FAQ=${GS_BUCKET}/All_CPTAC_customDB.fastq
GS_IDX=${GS_BUCKET}/${REF%%fa.gz}kal.idx
GS_OUT=${GS_BUCKET}/output
GS_LOG=${GS_BUCKET}/log

# Run kallisto quantification 
dsub \
   --name kallisto_quant \
   --project ${GS_PROJECT} \
   --zones 'us-*' \
   --image "${DOCKER_IMAGE}" \
   --input "IDX=${GS_IDX}" \
   --input "FASTQ=${GS_FAQ}" \
   --output-recursive "KALOUT=${GS_OUT}" \
   --logging ${GS_BUCKET}/log \
   --min-cores 8 \
   --command 'kallisto quant -i ${IDX} -o ${KALOUT} -b 100 --single -l 180 -s 20 -t 8 ${FASTQ}' \
   --wait
