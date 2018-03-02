#
# Example of a WDL workflow fro quantifying transcripts in 
# single read BAM files. 


# Fetch the reference FASTA formatted file
task getFasta {
   String RefUrl 

   command <<<
      curl ${RefUrl} > ref.fasta
   >>>

   output {
      File RefFasta = "ref.fasta"
      }
}

# Index a FASTA formatted file of target sequences
task indexFasta {
   File RefFasta

   command <<<
      kallisto index -i "kal.idx" ${RefFasta}
   >>>

   runtime {
      docker : "nareshr/kallisto:v0.43"
      memory: "15G"
      disks: "local-disk 500 HDD"
      cpu: "4"
      }

   output {
      File RefIndex = "kal.idx"
      }
   }

# Convert BAM to fastq
task convertBam {
   File ReadsBam

   command <<<
      bedtools bamtofastq -i ${ReadsBam} -fq reads.fastq
   >>>

   runtime {
      docker : "biocontainers/bedtools"
      memory: "15G"
      disks: "local-disk 500 HDD"
      cpu: "4"
      }

   output {
      File ReadsFastq = "reads.fastq"
   }
}

# Runs the quantification algorithm
task quant {
   File ReadsFastq
   File RefIndex

   command <<<
      kallisto quant -i ${RefIndex} -o kalout --single -l 180 -s 20 -t 8 ${ReadsFastq}
   >>>

   runtime {
      docker : "nareshr/kallisto:v0.43"
      memory: "15G"
      disks: "local-disk 500 HDD"
      cpu: "4"
      }

   output {
      File Abundance = "kalout/abundance.tsv"
      File RunInfo   = "kalout/run_info.json"
   }
}

# Actual workflow
workflow kallisto {
   Array[File] readsBam

   call getFasta {
      input: RefUrl="ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz"
   }
   call indexFasta {
      input: RefFasta=getFasta.RefFasta
      }

   scatter (bam in readsBam) {
      call convertBam {
         input: ReadsBam=bam
         }
      call quant {
         input: 
            ReadsFastq=convertBam.ReadsFastq,
            RefIndex=indexFasta.RefIndex
         }
      }

   output {
      File RefIndex = indexFasta.RefIndex
      Array[File] Abundance = quant.Abundance
      Array[File] RunInfo = quant.RunInfo
      }
   }
   
