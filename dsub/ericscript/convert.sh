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

# Logs
GS_LOG="${GS_PREFIX}/logs"

# Docker images
SAMTOOLS_IMAGE=quay.io/biocontainers/samtools:1.5--2

# Convert bam to fastq
dsub \
   --name ericscript_convert \
   --project "${GS_PROJECT}" \
   --image "${SAMTOOLS_IMAGE}" \
   --zones 'us-*' \
   --logging "${GS_LOG}/convert" \
   --min-cores 4 \
   --command 'samtools fastq -c 6 -@ 4 -1 ${ES_FQ1_FILE} -2 ${ES_FQ2_FILE} ${ES_BAM_FILE}' \
   --skip \
   --tasks $1
