# Quantification with kallisto

In the following example we use the program [kallisto](https://pachterlab.github.io/kallisto/about.html) by the Pachter Lab to quantify transcript abundances from RNA-Seq data from high-throughput sequencing reads provided in the ISB CGC repository, all being done using the Google Cloud.   The data will reside in Google Cloud Storage buckets and execution of kallisto is performed by submitting a qsub job to Google Pipelines.   The job executes on a Google Compute virtual machine with input files downloaded and output pushed to Google Cloud Storage.

## Setup

* Follow the [dsub getting started](https://github.com/googlegenomics/dsub/blob/master/README.md#getting-started)  instructions.  This can be installed on a local desktop, a Google Compute virtual machine, or in the Google console.
* Next  set environment variables used by commands in this example:

	* ```export GS_BUCKET=YOUR_BUCKET_NAME    # e.g. gs://isb-cgc-kallisto-example```
	* ```export GS_PROJECT=YOUR_PROJECT_ID```

* The last item you will need is a input/output Google Storage bucket for input/output files.  If you don't already have one you can create one using either the Google Console or the command line using the command:

### Build an index
Prior to running kallisto on sequencing reads, one has to first create index of the reference transcriptome file.  This index file is used to perform rapid pseudoalignments on the reads instead of having to align all the reads against the reference.  Once created this index file can be used to analyze as many sample files as desired.  For this example we are going to use the Genome Reference Consortium Human Build 37 from Ensembl.  This file consists of transcript sequences for actual and possible genes, including pseudogenes, NMD and the like, and is available from the FTP site.  

So we'll first need to copy the file into our Google Storage Bucket:
```
curl -L ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz | gsutil cp - ${GS_BUCKET}/Homo_sapiens.GRCh37.cdna.all.fa.gz
```
Next submit a job to index the reference transcriptome file using kallisto.  For this command we'll be using a docker image already created and uploaded in [DockerHub](https://hub.docker.com/r/nareshr/kallisto/).  We'll also request that the job submission wait (using the --wait option) until the job completes before returning back the prompt.  The input fasta file identifed by KALREF will be downloaded from your bucket and the resulting output file uploaded to your bucket in the location identified by KALIDX.

```
dsub \
   --name kallisto_index \
   --project ${GS_PROJECT} \
   --zones 'us-*' \
   --image "nareshr/kallisto:v0.43" \
   --input "KALREF=${GS_BUCKET}/Homo_sapiens.GRCh37.cdna.all.fa.gz" \
   --output "KALIDX=${GS_BUCKET}/Homo_sapiens.GRCh37.cdna.all.kal.idx" \
   --logging ${GS_BUCKET}/log \
   --command 'kallisto index -i ${KALIDX} ${KALREF}' \
   --min-ram 16 \
   --wait
```

You should see output that looks like:

```
Job: kallisto-i--<user-id>--170808-225609-25
Launched job-id: kallisto-i--<user-id>--170808-225609-25
To check the status, run:
  dstat --project <project-id> --jobs kallisto-i--<user-id>--170808-225609-25 --status '*'
To cancel the job, run:
  ddel --project <project-id> --jobs kallisto-i--<user-id>--170808-225609-25
Waiting for job to complete...
Waiting for: kallisto-i--<user-id>--170808-225609-25.
```

### Conversion to fastq

RNA-Seq data in ISB-CGC are stored in bam files, a binary format for storing sequencing data.  For this example we're going to use a publicly available bam files for colorectal cancer from the CPTAC project.  The program kalllisto however only reads fastq files, which is a text-based format for storing both a sequence and its corresponding quality scores.  Therefore before being able to use kallisto we'll first need to convert the bam files into fastq.  For this we'll use the program bedtools which contains the command bedtofastq.  Like before, we'll also use a pre-build docker container from DockerHub and submit the job to Google pipelines.  The command is:

```
dsub \
   --name kallisto_convert \
   --project ${GS_PROJECT} \
   --zones 'us-*' \
   --image "biocontainers/bedtools" \
   --input "BAM=gs://isb-cgc-open/CPTAC/Phase_II/TCGA_Colorectal_Cancer_S_022/TCGA_Colorectal_Cancer_proBAM_PSM_genome_mapping_files/All_CPTAC_customDB.bam" \
   --output "FASTQ=${GS_BUCKET}/All_CPTAC_customDB.fastq" \
   --logging ${GS_BUCKET}/log \
   --command 'bedtools bamtofastq -i ${BAM} -fq ${FASTQ}' \
   --wait
```

### Quantification

Finally to quantify the transcripts, submit a kallisto job as follows:

```
dsub \
   --name kallisto_quant \
   --project ${GS_PROJECT} \
   --zones 'us-*' \
   --image "nareshr/kallisto:v0.43" \
   --input "KALIDX=${GS_BUCKET}/Homo_sapiens.GRCh37.cdna.all.kal.idx" \
   --input "FASTQ=${GS_BUCKET}/All_CPTAC_customDB.fastq" \
   --output-recursive "KALOUT=${GS_BUCKET}/output" \
   --logging ${GS_BUCKET}/log \
   --min-cores 8 \
   --command 'kallisto quant -i ${KALIDX} -o ${KALOUT} -b 100 --single -l 180 -s 20 -t 8 ${FASTQ}' \
   --wait
```

### Check the results

To check the results, first list the output files with the following command:

```
gsutil ls ${GS_BUCKET}/output
```

The output should look like:

```
gs://cgc-kallisto-example/output/abundance.h5
gs://cgc-kallisto-example/output/abundance.tsv
gs://cgc-kallisto-example/output/run_info.json
```

The json file contains details about the kalliso run, the transcript abundances are in the tsv file, and the h5 file can be used to load the data into HDF.  To list the top 10 abundant transcripts called, run:

```
gsutil cat ${GS_BUCKET}/output/abundance.tsv | sort -t$'\t' -k4nr | head -10
```

Which should output:

```
ENST00000511370.1       1519    1340    42340.7 7349.25
ENST00000358087.5       3399    3220    29309.2 2117.09
ENST00000484216.1       403     224     26359.5 27370.3
ENST00000335295.4       628     449     24001.3 12433.1
ENST00000493945.1       1236    1057    23664.5 5207.29
ENST00000396859.1       1256    1077    19926   4303.23
ENST00000224237.5       1868    1689    15417   2123.05
ENST00000375167.1       1373    1194    14634.1 2850.69
ENST00000359193.2       2239    2060    14422.4 1628.39
ENST00000479986.1       522     343     14341.7 9725.14
```

