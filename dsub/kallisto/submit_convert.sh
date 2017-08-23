#!/bin/bash

#
# Usage: ./submit_convert.sh
#
set -euo pipefail 
set -x

# Google project 
GS_PROJECT=${GS_PROJECT:-XXXXXXX}

# Docker image to use
DOCKER_IMAGE=nareshr/kallisto:v0.43

# Storage locations
GS_BUCKET=${GS_BUCKET:-gs://cgc-kallisto-example}
GS_LOG=${GS_BUCKET}/log

# BAM 
GS_BAM=${BAM:-gs://isb-cgc-open/CPTAC/Phase_II/TCGA_Colorectal_Cancer_S_022/TCGA_Colorectal_Cancer_proBAM_PSM_genome_mapping_files/All_CPTAC_customDB.bam}
GS_FQ=${GS_BUCKET}/All_CPTAC_customDB.fastq

# Convert using bedtools
dsub \
   --name kallisto_convert \
   --project ${GS_PROJECT} \
   --zones 'us-*' \
   --image "biocontainers/bedtools" \
   --input "BAM=${GS_BAM}" \
   --output "FASTQ=${GS_FQ}" \
   --logging "${GS_LOG}" \
   --command 'bedtools bamtofastq -i ${BAM} -fq ${FASTQ}' \
   --wait
