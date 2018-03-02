#!/bin/bash
#
# Usage: 
#
#    dsub_ericscript.sh <BATCH TSV FILE>
#
: ${1?"Missing batch file"}

set -euo pipefail 
set -x

# Google project 
GS_PROJECT=${GS_PROJECT:-XXXXXXX}

# Storage locations
GS_BUCKET=${GS_BUCKET:-isb-cgc-examples}

# Storage prefix
GS_PREFIX=gs://${GS_BUCKET}/ericscript

# Inputs
GS_REF="${GS_REF:-gs://isb-cgc-misc/reference-data/ericscript/db_homosapiens_ensembl73/}"
GS_LOG="${GS_LOG:-${GS_PREFIX}/logs}"

# Docker images
ERICSCRIPT_IMAGE=cgrlab/ericscript

# Run ericscript
dsub \
   --name ericscript_run \
   --project ${GS_PROJECT} \
   --image "${ERICSCRIPT_IMAGE}" \
   --zones 'us-*' \
   --input-recursive "ES_REF=${GS_REF}" \
   --logging "${GS_LOG}/ericscript/logs" \
   --min-cores 8 \
   --min-ram 20 \
   --disk-size 400 \
   --command 'rmdir ${RESULT}; \
              /opt/EricScript/ericscript.pl -db ${ES_REF} -o ${RESULT} ${ES_FQ1_FILE} ${ES_FQ2_FILE} -p 8; \
             ' \
   --tasks $1 \
   --skip 
