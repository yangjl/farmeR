#' \code{Run GATK job on farm}
#'
#' GATK Best Practices: recommended workflows for variant discovery analysis.
#'
#' see more detail about GATK:
#' \url{https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1}
#'
#' @param fq An input data.frame for fastq files. Must contains fq1, fq2 and out.
#' @param kitpath The absolute or relative path of the fermi.kit directory that can invoke the pipeline.
#' @param genome The full path of genome with bwa indexed reference fasta file.
#' @param s Approximate genome size, default=3g.
#' @param t Number of threads, default=16.
#' @param l Primary read length, default=100.
#' @param arrayjobs A character specify the number of array you try to run, i.e. 1-100.
#' @param jobid The job name show up in your sq NAME column.
#' @param email Your email address that farm will email to once the job was done/failed.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' fq <- data.frame(fq1=c("f_1.fq", "t_1.fq"), fq2=c("f_1.fq", "t_2.fq"), out=c("t1", "t2"))
#' run_fermikit(fq, kitpath="/home/jolyang/bin/fermikit/",
#' genome="/home/jolyang/dbcenter/AGP/AGPv2", s='3g', t=16, l=100, arrayjobs="1-2",
#' jobid="fermi", email=NULL)
#'
#' @export
run_GATK <- function(fq,
                         kitpath="/home/jolyang/bin/fermikit",
                         s='3g', t=16, l=100,
                         arrayjobs="1-2",
                         jobid="fermi",
                         email=NULL){

  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:nrow(fq)){

    shid <- paste0("slurm-script/run_bwamem_", i, ".sh")

    #Generate a SAM file containing aligned reads
    #http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
    rg <- "@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1"

    cat("### GATK pipeline created by farmeR",
        paste("###", format(Sys.time(), "%a %b %d %X %Y")),
        file=shid, sep="\n", append=FALSE)

    if(bwa){
      cat(paste("###Generate a SAM file containing aligned reads"),
          paste("bwa mem -M -R", rg, "-p", ref.fa, fq1, fq2, ">", sam),

          "java -jar picard.jar SortSam \\",
          paste0("    INPUT=", aligned_sam, " \\"),
          paste0("    OUTPUT=", sorted_bam, " \\"),
                 "    SORT_ORDER=coordinate",



          file=shid, sep="\n", append=TRUE)
    }
    if(markDup){
      cat("",
          "java -jar picard.jar MarkDuplicates \\",
          paste0("INPUT=", sorted_bam, " \\"),
          paste0("OUTPUT=", dedup_bam, " \\"),
          paste0("METRICS_FILE=metrics.txt"),
          paste0(""),
          "java -jar picard.jar BuildBamIndex \\",
          paste0("INPUT=dedup_reads.bam"),

          file=shid, sep="\n", append=TRUE)

    }

    if(realignInDels){
      cat(paste0("java -Xmx4g -Djava.io.tmpdir=/path/to/tmpdir \\"),
          paste0("-jar /path/to/GenomeAnalysisTK.jar \\"),
          paste0("-I <lane-level.bam> \\"),
          paste0("-R <ref.fasta> \\"),
          paste0("-T IndelRealigner \\"),
          paste0("-targetIntervals <intervalListFromStep1Above.intervals> \\"),
          paste0("-o <realignedBam.bam> \\"),
          paste0("-known /path/to/indel_calls.vcf"),
          paste0("--consensusDeterminationModel KNOWNS_ONLY \\"),
          paste0("-LOD 0.4"),
          file=shid, sep="\n", append=FALSE)
    }
    if(recalBases){
      cat(paste0("java -Xmx4g -Djava.io.tmpdir=/path/to/tmpdir \\"),
          paste0("-jar /path/to/GenomeAnalysisTK.jar \\"),
          paste0("-I <lane-level.bam> \\"),
          paste0("-R <ref.fasta> \\"),
          paste0("-T IndelRealigner \\"),
          paste0("-targetIntervals <intervalListFromStep1Above.intervals> \\"),
          paste0("-o <realignedBam.bam> \\"),
          paste0("-known /path/to/indel_calls.vcf \\"),
          paste0("--consensusDeterminationModel KNOWNS_ONLY \\"),
          paste0("-LOD 0.4"),
          file=shid, sep="\n", append=FALSE)
    }

    cat("java	–jar	GenomeAnalysisTK.jar	–T	HaplotypeCaller	\\",
        paste0("–R	human.fasta	\\"),
        paste0("–I	sample1.bam	\\"),
        paste0("–o	sample1.g.vcf	\\"),
        paste0("–L	exome_targets.intervals	\\"),
        paste0("–ERC	GVCF	\\"),

        file=shid, sep="\n", append=FALSE)



  }

  shcode <- paste("sh slurm-script/run_bwamem_$SLURM_ARRAY_TASK_ID.sh", sep="\n")

  set_array_job(shid="slurm-script/run_bwamem_array.sh",
                shcode=shcode, arrayjobs=arrayjobs,
                wd=NULL, jobid=jobid, email=email)
}





"java -Xmx64g -jar ~/Programs/GenomeAnalysisTK.jar
-R /reference.fasta
-T GenotypeGVCFs
-o woth300-rawsnps.vcf
-nt 24
--variant sample1.g.vcf --variant sample2.g.vcf ... ...
--variant sample296.g.vcf"

