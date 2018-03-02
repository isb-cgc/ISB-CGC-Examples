# Gene fusion discovery with EricScript


This example shows how one might discover gene fusions in paired end RNA-seq data stored in ISB-CGC cloud storage.  To do this we'll be using the program [EricScript](https://sites.google.com/site/bioericscript/).   The data will reside in Google Cloud Storage buckets and the execution of EricScript is performed by submitting a pipeline to Google Pipelines using dsub.   The job executes on a Google Compute virtual machine with input files downloaded and output pushed to Google Cloud Storage.
<pre>
<b>NOTE:</b> This example will incur some charges as there are multiple computational analyses performed.
</pre>

## Setup

#### Install dsub
Follow the [dsub getting started](https://github.com/googlegenomics/dsub#getting-started)  instructions.  This can be installed on a local desktop, a Google Compute virtual machine, or in the Google Shell.

#### Setup the environment

Next  set environment variables used by commands in this example.  The first is the Google project and the second is the name of the Google bucket to store the data.

```
export GS_PROJECT=<YOUR PROJECT ID>
export GS_BUCKET=<YOUR_BUCKET_NAME>    # e.g. gs://isb-cgc-ericscript-example
```

#### Create the storage bucket
If the storage bucket don't already have one you can create one using either the Google Console or the command line using the command:

```gsutil mb gs://${GS_BUCKET}```

#### Prepare the database

You'll need a reference Ensembl Database in Google storage in order to be able to run EricScript using dsub in the Google cloud.  Usually you'd use the --printdb command of EricScript to list the available databases and the --downdb command to download a particular database and copy this into Google storage.  But in the current version (0.5.5) these don't  work correctly.

The program authors have also placed a [download](https://sites.google.com/site/bioericscript/download) section on their website for retrieving a database.  The files are currently (as of the time of this writing) hosted from a Google Drive and therefore can only easily be downloaded using a browser and with a Google identity.  Therefore we've provided the Human database in a ISB-CGC publicly available  Google storage [bucket](https://console.cloud.google.com/storage/browser/isb-cgc-misc/reference-data/ericscript) to simplify the example.  The easiest way to use this database is to setup a environment variable that'll refer to the reference database as so:

```
export ES_REF=gs://isb-cgc-misc/reference-data/ericscript/db_homosapiens_ensembl73/
```

## Select Some RNA Seq Data

Next we'll need some RNA Seq data. 

#### Using ISB-CGC

Go to the ISB-CGC website [https://isb-cgc.appspot.com](https://isb-cgc.appspot.com) and log on.  On the main page, click on the "hamburger" in Menu in the top right and choose Cohorts > Create Cohort - Filters.   Click on the CCLE DATA tab, and apply the following additional filters:

 - For program select "CCLE"
 - For gender select "Male"
 - For site primary select "Lung"
 - For Sample Type select "Cell Lines"
 - For histological subtype select "Large Cell Carcinoma"

This should limit you to 12 number of samples/cases.  Go ahead and click save to save it as a new cohort named "ericscript-example".  Then when the list of cohorts are displayed click on the cohort name for the cohort you just created.  Click "File Browser" in the list of buttons at the top of the cohort details page.   In the filter area of the file browser, select 'RNA-Seq' under experimental strategy.  This should give you a list of 12 RNA-Seq bam files.  Go ahead and click the CSV button to download the list, which will contain a column labeled "Cloud Storage Location" containing the Google Storage location for the files.  Save the file in your work area as "cohort.tsv".

## Conversion to fastq

RNA-Seq data in ISB-CGC are stored in bam files, a binary format for storing sequencing data.  Like the kalisto [example](https://github.com/isb-cgc/ISB-CGC-Examples/tree/master/dsub/kallisto) we'll first need to convert the bam files into fastq as that's what EricScript uses as input. We'll use using the program bedtools which contains the command bedtofastq via an existing pre-build docker container and submit the job to Google Pipelines using dsub.  However since we have multiple bed files to convert, we'll submit it as a [dsub batch job](https://github.com/googlegenomics/dsub#submitting-a-batch-job) to convert all the files at once instead of having to convert each file one at a time.

#### Creating the bedtools batch file

Using the previously saved file containing your cohort samples, remove the first line and all columns other than the sample column and the cloud storage location column.  Rename the first column from sample to the string "--output ES_FQ1_FILE" and then make an exact copy of this column, naming it "--output ES_FQ2_FILE".   Lastly rename the header Cloud Storage Location to the string "--input ES_BAM_FILE".   

This will be the tab delimited file you'll use to submit the batch job.  For each line in this file dsub will create a new task in its own separate Google VM to execute a given command using the output paths provided in the two output columns and the input file provided in the input column.  In the command the input and output variables will have the values of the corresponding row.  

Tip: the following little command line magic can do all of this for you.
```
awk -v gs="gs://${GS_BUCKET}/ericscript" -F $'\t' \
  'BEGIN { OFS = FS } \
   NR == 1 { print "--input ES_BAM_FILE", \
                   "--output ES_FQ1_FILE", \
                   "--output ES_FQ2_FILE" } \
   NR > 2 { print $8, gs"/"$1"/reads_1.fastq.gz", \
                      gs"/"$1"/reads_2.fastq.gz" } \
  ' cohort.tsv > convert.tsv
```

The end product should look like something like this:
```
--input ES_BAM_FILE     --output ES_FQ1_FILE    --output ES_FQ2_FILE
gs://isb-cgc-open/NCI-GDC/legacy/CCLE/CCLE-LUSC/RNA-Seq/Aligned_reads/01520b32-7a1c-43c5-9471-a7bb2eac9e94/G41716.IA-LM.5.bam   gs:///ericscript/CCLE-IA-LM/reads_1.fastq.gz    gs:///ericscript/CCLE-IA-LM/reads_2.fastq.gz
gs://isb-cgc-open/NCI-GDC/legacy/CCLE/CCLE-LUSC/RNA-Seq/Aligned_reads/f62a7d4b-1d78-4471-b366-71ab9b296720/G41690.NCI-H460.5.bam        gs:///ericscript/CCLE-NCI-H460/reads_1.fastq.gz gs:///ericscript/CCLE-NCI-H460/reads_2.fastq.gz
gs://isb-cgc-open/NCI-GDC/legacy/CCLE/CCLE-LUSC/RNA-Seq/Aligned_reads/10b2db78-b409-4074-935f-bccdc20ac756/G30580.T3M-10.1.bam  gs:///ericscript/CCLE-T3M-10/reads_1.fastq.gz   gs:///ericscript/CCLE-T3M-10/reads_2.fastq.gz
gs://isb-cgc-open/NCI-GDC/legacy/CCLE/CCLE-LUSC/RNA-Seq/Aligned_reads/d77e8809-d688-4bb1-8041-6d774c50c1ca/G28534.NCI-H810.1.bam        gs:///ericscript/CCLE-NCI-H810/reads_1.fastq.gz gs:///ericscript/CCLE-NCI-H810/reads_2.fastq.gz
gs://isb-cgc-open/NCI-GDC/legacy/CCLE/CCLE-LUSC/RNA-Seq/Aligned_reads/cbf2f7c7-ded1-4ff0-afc8-222ec3c404f2/G41704.NCI-H1581.5.bam       gs:///ericscript/CCLE-NCI-H1581/reads_1.fastq.gz        gs:///ericscript/CCLE-NCI-H1581/reads_2.fastq.gz
gs://isb-cgc-open/NCI-GDC/legacy/CCLE/CCLE-LUSC/RNA-Seq/Aligned_reads/427b9dee-4905-4b92-98b0-c7b700e8986c/G28041.LCLC-97TM1.1.bam      gs:///ericscript/CCLE-LCLC-97TM1/reads_1.fastq.gz       gs:///ericscript/CCLE-LCLC-97TM1/reads_2.fastq.gz
gs://isb-cgc-open/NCI-GDC/legacy/CCLE/CCLE-LUSC/RNA-Seq/Aligned_reads/ff04194f-7c26-4cb1-a763-04ef0f377752/G27304.HCC-1438.1.bam        gs:///ericscript/CCLE-HCC-1438/reads_1.fastq.gz gs:///ericscript/CCLE-HCC-1438/reads_2.fastq.gz
gs://isb-cgc-open/NCI-GDC/legacy/CCLE/CCLE-LUSC/RNA-Seq/Aligned_reads/e48ea2ee-1dda-4061-a199-6e22fd2df382/G25212.NCI-H661.1.bam        gs:///ericscript/CCLE-NCI-H661/reads_1.fastq.gz gs:///ericscript/CCLE-NCI-H661/reads_2.fastq.gz
gs://isb-cgc-open/NCI-GDC/legacy/CCLE/CCLE-LUSC/RNA-Seq/Aligned_reads/e6075d40-df96-4457-9350-cc632931ec39/G28061.LCLC-103H.1.bam       gs:///ericscript/CCLE-LCLC-103H/reads_1.fastq.gz        gs:///ericscript/CCLE-LCLC-103H/reads_2.fastq.gz
gs://isb-cgc-open/NCI-GDC/legacy/CCLE/CCLE-LUSC/RNA-Seq/Aligned_reads/6108dcae-c84c-43ca-b33b-065ee422eea0/G41697.COR-L23.5.bam gs:///ericscript/CCLE-COR-L23/reads_1.fastq.gz  gs:///ericscript/CCLE-COR-L23/reads_2.fastq.gz
gs://isb-cgc-open/NCI-GDC/legacy/CCLE/CCLE-LUSC/RNA-Seq/Aligned_reads/112eb28d-52a1-43eb-8805-640b21056879/G28592.NCI-H1155.1.bam       gs:///ericscript/CCLE-NCI-H1155/reads_1.fastq.gz        gs:///ericscript/CCLE-NCI-H1155/reads_2.fastq.gz
gs://isb-cgc-open/NCI-GDC/legacy/CCLE/CCLE-LUSC/RNA-Seq/Aligned_reads/8c0d94ae-a908-47c9-a301-3b5c3d0fafaa/G28021.LU99.1.bam    gs:///ericscript/CCLE-LU99/reads_1.fastq.gz     gs:///ericscript/CCLE-LU99/reads_2.fastq.gz
```
#### Submitting the conversion batch job

Running all the tasks in the batch file only requires one simple command:
```
dsub \
   --name ericscript_convert \
   --project ${GS_PROJECT} \
   --image "quay.io/biocontainers/samtools:1.5--2" \
   --zones 'us-*' \
   --logging "gs://${GS_BUCKET}/ericscript/logs" \
   --min-cores 4 \
   --command 'samtools fastq -c 6 -@ 4 -1 ${ES_FQ1_FILE} -2 ${ES_FQ2_FILE} ${ES_BAM_FILE}' \
   --tasks convert.tsv \
   --skip \
   --wait
```
This command will launch a task for each task in the tsv file, filling in the ES_FQ1_FILE, ES_FQ2_FILE, and ES_BAM_FILE variables in the command for each task.  The command will skip any results already in cloud storage, and will wait till all the tasks are complete before finishing.  Its possible to save a little money by using the '--preemptible' flag, but if a VM gets preempted you'll have to rerun the dsub command over again.


## Running EricScript

#### Creating the EricScript batch file
Similarly to the conversion, we're going to create a tab-delimited task file to invoke EricScript on each of the sets of fasta files created in the conversion step.  Starting from the conversion batch file, change the output fastq files to be inputs, by changing "--output" to "--input".  Remove the bam input column, and duplicate one of the fastq columns.  Change the heading of this new column to "--output-recursive RESULT", and change the fasta file "reads_1.fastq.gz" to be "results".  Some linux command line magic to do all of this is:

```
awk -F $'\t' 'BEGIN { OFS = FS } \
NR == 1 { print "--input ES_FQ1_FILE","--input ES_FQ2_FILE","--output-recursive RESULT" } \
NR > 1 { a=$3; sub(/reads_2.fastq.gz/,"results",a); \
         print $2,$3,a }' convert.tsv > ericscript.tsv
```
When done you should have a file named ericscript.tsv that looks like:

```
--input ES_FQ1_FILE     --input ES_FQ2_FILE     --output-recursive RESULT
gs:///ericscript/CCLE-IA-LM/reads_1.fastq.gz    gs:///ericscript/CCLE-IA-LM/reads_2.fastq.gz    gs:///ericscript/CCLE-IA-LM/results
gs:///ericscript/CCLE-NCI-H460/reads_1.fastq.gz gs:///ericscript/CCLE-NCI-H460/reads_2.fastq.gz gs:///ericscript/CCLE-NCI-H460/results
gs:///ericscript/CCLE-T3M-10/reads_1.fastq.gz   gs:///ericscript/CCLE-T3M-10/reads_2.fastq.gz   gs:///ericscript/CCLE-T3M-10/results
gs:///ericscript/CCLE-NCI-H810/reads_1.fastq.gz gs:///ericscript/CCLE-NCI-H810/reads_2.fastq.gz gs:///ericscript/CCLE-NCI-H810/results
gs:///ericscript/CCLE-NCI-H1581/reads_1.fastq.gz        gs:///ericscript/CCLE-NCI-H1581/reads_2.fastq.gz        gs:///ericscript/CCLE-NCI-H1581/results
gs:///ericscript/CCLE-LCLC-97TM1/reads_1.fastq.gz       gs:///ericscript/CCLE-LCLC-97TM1/reads_2.fastq.gz       gs:///ericscript/CCLE-LCLC-97TM1/results
gs:///ericscript/CCLE-HCC-1438/reads_1.fastq.gz gs:///ericscript/CCLE-HCC-1438/reads_2.fastq.gz gs:///ericscript/CCLE-HCC-1438/results
gs:///ericscript/CCLE-NCI-H661/reads_1.fastq.gz gs:///ericscript/CCLE-NCI-H661/reads_2.fastq.gz gs:///ericscript/CCLE-NCI-H661/results
gs:///ericscript/CCLE-LCLC-103H/reads_1.fastq.gz        gs:///ericscript/CCLE-LCLC-103H/reads_2.fastq.gz        gs:///ericscript/CCLE-LCLC-103H/results
gs:///ericscript/CCLE-COR-L23/reads_1.fastq.gz  gs:///ericscript/CCLE-COR-L23/reads_2.fastq.gz  gs:///ericscript/CCLE-COR-L23/results
gs:///ericscript/CCLE-NCI-H1155/reads_1.fastq.gz        gs:///ericscript/CCLE-NCI-H1155/reads_2.fastq.gz        gs:///ericscript/CCLE-NCI-H1155/results
gs:///ericscript/CCLE-LU99/reads_1.fastq.gz     gs:///ericscript/CCLE-LU99/reads_2.fastq.gz     gs:///ericscript/CCLE-LU99/results
```

#### Submitting the EricScript batch job

Running all EricScript tasks is similar to running the conversion tasks, one simple command:

```
dsub \
   --name ericscript_run \
   --project ${GS_PROJECT} \
   --image "cgrlab/ericscript" \
   --zones 'us-*' \
   --input-recursive "ES_REF=${ES_REF}" \
   --logging "gs://${GS_BUCKET}/ericscript/logs" \
   --min-cores 8 \
   --min-ram 20 \
   --disk-size 400 \
   --command 'rmdir ${RESULT}; \
              /opt/EricScript/ericscript.pl -db ${ES_REF} -o ${RESULT} ${ES_FQ1_FILE} ${ES_FQ2_FILE} -p 8; \
             ' \
   --tasks ericscript.tsv \
   --skip \
   --wait
```

You should see output that looks like:

```
Job: ericscript--<user>--180301-185226-88
Launched job-id: ericscript--<user>--180301-185226-88
12 task(s)
To check the status, run:
  dstat --project >project> --jobs 'ericscript--<user>--180301-185226-88' --status '*'
To cancel the job, run:
  ddel --project <project> --jobs 'ericscript--<user>--180301-185226-88'
ericscript--<user>--180301-185226-88
Waiting for job to complete...
Waiting for: ericscript--jslagel--180301-185226-88
```


### Check the results

To check the results, first list the result folders for all samples with the following command:

```
gsutil ls gs://${GS_BUCKET}/ericscript/*/results
```

The output should look something like:

```
gs://isb-cgc-examples/ericscript/CCLE-IA-LM/results/:
gs://isb-cgc-examples/ericscript/CCLE-IA-LM/results/MyEric.Summary.RData
gs://isb-cgc-examples/ericscript/CCLE-IA-LM/results/MyEric.results.filtered.tsv
gs://isb-cgc-examples/ericscript/CCLE-IA-LM/results/MyEric.results.total.tsv
gs://isb-cgc-examples/ericscript/CCLE-IA-LM/results/aln/
gs://isb-cgc-examples/ericscript/CCLE-IA-LM/results/out/

gs://isb-cgc-examples/ericscript/CCLE-NCI-H460/results/:
gs://isb-cgc-examples/ericscript/CCLE-NCI-H460/results/MyEric.Summary.RData
gs://isb-cgc-examples/ericscript/CCLE-NCI-H460/results/MyEric.results.filtered.tsv
gs://isb-cgc-examples/ericscript/CCLE-NCI-H460/results/MyEric.results.total.tsv
gs://isb-cgc-examples/ericscript/CCLE-NCI-H460/results/aln/
gs://isb-cgc-examples/ericscript/CCLE-NCI-H460/results/out/
...
```

To example the results for a single sample, use something like the following command:
```
gsutil cat gs://${GS_BUCKET}/ericscript/CCLE-IA-LM/results/MyEric.results.total.tsv | head -10
```

Which should output:

```
GeneName1       GeneName2       chr1    Breakpoint1     strand1 chr2    Breakpoint2     strand2 EnsemblGene1    EnsemblGene2    crossingreads   spanningreads    mean.insertsize homology        fusiontype      Blacklist       InfoGene1       InfoGene2       JunctionSequence        GeneExpr1GeneExpr2       GeneExpr_Fused  ES      GJS     US      EricScore
PCMTD2  RCBTB1  20      62919222        +       13      50120753        -       ENSG00000203880 ENSG00000136144 3       38      44.16           inter-chromosomal                protein-L-isoaspartate (D-aspartate) O-methyltransferase domain containing 2 [Source:HGNC Symbol;Acc:15882]     regulator of chromosome condensation (RCC1) and BTB (POZ) domain containing protein 1 [Source:HGNC Symbol;Acc:18243]     aacatttggaattttcgtatacgtggattctagaggggtgacagcgaaacTGGGGCCGTCCAGTTGTAGGAAAACAAGCTCAGGACTCCCACTGATTCTA     39.44   8.89    83.99   0.9052  0.8969  0.0789473684210526      0.355265882261857
GTF2I   STEAP2  7       74171037        +       7       89863859        +       ENSG00000077809 ENSG00000157214 5       4       22.15   More than 30 homologies found: ENSG00000079950 (49%), ENSG00000134253 (48%), ENSG00000175224 (48%), ENSG00000160049 (48%), ENSG00000003756 (44%), ENSG00000262826 (48%), ENSG00000143624 (48%), ENSG00000056277 (42%), ENSG00000130244 (48%), ENSG00000189051 (42%), ENSG00000137814 (36%), ENSG00000186854 (40%), ENSG00000134255 (40%), ENSG00000269198 (38%), ENSG00000102032 (38%), ENSG00000105220 (38%), ENSG00000165323 (26%), ENSG00000187446 (25%), ENSG00000181798 (25%), ENSG00000171777 (25%), ENSG00000163359 (25%), ENSG00000132780 (25%), ENSG00000131375 (25%), ENSG00000126934 (25%), ENSG00000231825 (25%), ENSG00000231370 (25%), ENSG00000230873 (25%), ENSG00000226618 (25%), ENSG00000225748 (25%), ENSG00000225164 (25%)      intra-chromosomal                general transcription factor IIi [Source:HGNC Symbol;Acc:4659]  STEAP family member 2, metalloreductase [Source:HGNC Symbol;Acc:17885]   agttttatggtaattactgtgaattagtgattctgagggacatatcagaaTCTTTTTTTTTTTTTTTTTTTTTAAAGACAGAGTCTTGCTCTGTCACCCA    70.96   5.27    38.58   0.3005   0.2368  0.8     0.433960685058166
LRRFIP2 HTRA4   3       37190385        -       8       38845848        +       ENSG00000093167 ENSG00000169495 5       10      1.5     ENSG00000154258 (52%), ENSG00000167766 (52%), ENSG00000109861 (52%), ENSG00000197102 (47%), ENSG00000100592 (40%)        inter-chromosomal               leucine rich repeat (in FLII) interacting protein 2 [Source:HGNC Symbol;Acc:6703]        HtrA serine peptidase 4 [Source:HGNC Symbol;Acc:26909]  aagaccgattttctgcagaagatgaagctttgagtaacattgccagagagGGACATGGATGGAGCTGGAAACCATTATCCTCAGCAAACTAACGCAGGAA     5.71    4.53    14.98   0.8173  0.5807  0.5      0.886136913714164
LOX     NACA    5       121413795       -       12      57119135        -       ENSG00000113083 ENSG00000196531 4       1       17.51           inter-chromosomal                lysyl oxidase [Source:HGNC Symbol;Acc:6664]     nascent polypeptide-associated complex alpha subunit [Source:HGNC Symbol;Acc:7629]       ggtctgaatcccacccttggcattgcttggtggagactgagatacccgtgCTCCGCTCCAGGCTTCCTTCTGCAACAGGCGTGGGTCACGCTCTCGCTCG    3.47    54.47   761.59   0.7238  0.6486  0.25    0.701552057244671
LOX     ACTB    5       121413388       -       7       5570233 -       ENSG00000113083 ENSG00000075624 3       2       125.75          inter-chromosomal                lysyl oxidase [Source:HGNC Symbol;Acc:6664]     actin, beta [Source:HGNC Symbol;Acc:132]        gcctccgcccagcagccccgcactccgatcctgctgatccgcgacaaccgCACCGCCGAGACCGCGTCCGCCCCGCGAGCACAGAGCCTCGCCTTTGCCG     3.47    47.37   566.28  0.9817  0.8933  0.666666666666667       0.953766310689514
SLC38A1 RAD21   12      46582743        -       8       117859323       -       ENSG00000111371 ENSG00000164754 3       2       9.18            inter-chromosomal                solute carrier family 38, member 1 [Source:HGNC Symbol;Acc:13447]       RAD21 homolog (S. pombe) [Source:HGNC Symbol;Acc:9811]   tgactgggcctgctcatcgagtagtgacgaaggccactgaaacccgccgaGAAAAAAAAACAAACAAACACTTCTTTCCATCAGTAACACTGGCAATCTT    35.7    130.32  26.05   0.8829   0.5573  0.666666666666667       0.972783886952759
PTMA    NCL     2       232577226       +       2       232326275       -       ENSG00000187514 ENSG00000115053 7       99      61.11           intra-chromosomal                prothymosin, alpha [Source:HGNC Symbol;Acc:9623]        nucleolin [Source:HGNC Symbol;Acc:7667] atgaggaagctgagtcagctacgggcaagcgggcagctgaagatgatgagGATGACGATGACGATGAGGAAGATGACTCTGAAGAAGAAGCTATGGAGAC     1.7     19.74   591.74  0.4333  0.3367  0.0707070707070707       0.381102208803414
PTMA    GPR157  2       232577803       +       1       9162486 -       ENSG00000187514 ENSG00000180758 5       5       41.93   More than 30 homologies found: ENSG00000141696 (53%), ENSG00000158122 (65%), ENSG00000187024 (50%), ENSG00000153721 (50%), ENSG00000056972 (50%), ENSG00000087206 (65%), ENSG00000163827 (65%), ENSG00000262273 (49%), ENSG00000187021 (49%), ENSG00000125843 (49%), ENSG00000240720 (50%), ENSG00000135926 (50%), ENSG00000001630 (50%), ENSG00000198740 (50%), ENSG00000103351 (50%), ENSG00000171467 (65%), ENSG00000166938 (65%), ENSG00000105866 (50%), ENSG00000168172 (50%), ENSG00000176834 (65%), ENSG00000068878 (62%), ENSG00000175984 (50%), ENSG00000076513 (50%), ENSG00000271000 (50%), ENSG00000263100 (50%), ENSG00000262923 (50%), ENSG00000262832 (50%), ENSG00000262717 (50%), ENSG00000262680 (50%), ENSG00000262661 (50%)      inter-chromosomalprothymosin, alpha [Source:HGNC Symbol;Acc:9623]        G protein-coupled receptor 157 [Source:HGNC Symbol;Acc:23687]   ggaggaaaaaagaaccaaaacttccaaggccctgctttttttcttaaaagTTTTTTTTTTTTTTTTTGAGATGGAGTCTCACTCTATTGCCCAGGCTGGA     1.7     9.95    1535.5  0.8994  0.9447  1       0.950613037333006
PTMA    STC2    2       232576679       +       5       172755321       -       ENSG00000187514 ENSG00000113739 3       22      176.57          inter-chromosomal                prothymosin, alpha [Source:HGNC Symbol;Acc:9623]        stanniocalcin 2 [Source:HGNC Symbol;Acc:11374]  aatgaggtagacgaagaagaggaagaaggtggggaggaagaggaggaggaGGAAAAGGCGAGCAAAAAGGAAGAGTGGGAGGAGGAGGGGAAGCGGCGAA     1.7     58.16   254.96  0.3669  0.2694  0.136363636363636        0.529423884397523
```


