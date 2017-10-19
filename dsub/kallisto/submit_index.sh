#!/bin/bash

#
# Usage: ./submit_index.sh
#
set -euo pipefail 
set -x

# Google project 
#GS_PROJECT=

# Docker image to use
DOCKER_IMAGE=nareshr/kallisto:v0.43

# Reference database
REF=Homo_sapiens.GRCh37.cdna.all.fa.gz
REF_URL=ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/cdna/${REF}

# Storage locations
GS_BUCKET=${GS_BUCKET:-cgc-kallisto-example}
GS_REF=${GS_BUCKET}/${REF}
GS_IDX=${GS_BUCKET}/${REF%%fa.gz}kal.idx
GS_LOG=${GS_BUCKET}/log

# Copy Ensembl cDNA file to Google Storage
curl -L ${REF_URL} | gsutil cp - ${GS_REF}

# Run kallisto index
dsub \
   --name kallisto_index \
   --project ${GS_PROJECT} \
   --zones 'us-*' \
   --image "${DOCKER_IMAGE}" \
   --input "GS_REF=${GS_REF}" \
   --output "GS_IDX=${GS_IDX}" \
   --logging ${GS_LOG} \
   --command 'kallisto index -i ${GS_IDX} ${GS_REF}' \
   --min-ram 16 \
   --wait
