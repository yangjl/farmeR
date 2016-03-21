#' \code{Run BWA-mem job on farm}
#'
#' Fermikit is a de novo assembly based variant calling pipeline for Illumina short reads
#'
#' see more detail about fermikit by Li, Heng:
#' \url{https://github.com/lh3/fermikit}
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
run_bwamem <- function(fq,
                         kitpath="/home/jolyang/bin/fermikit",
                         s='3g', t=16, l=100,
                         arrayjobs="1-2",
                         jobid="fermi",
                         email=NULL){

  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:nrow(fq)){

    shid <- paste0("slurm-script/run_bwamem_", i, ".sh")


    #-t INT     number of threads [1]
    #-M         mark shorter split hits as secondary (for Picard/GATK compatibility)
    #-T INT     minimum score to output [30]
    # genome ftp://ftp.ensemblgenomes.org/pub/plants/release-22/fasta/zea_mays/dna/README
    cmd1 <- paste0("bwa mem -t 18 -T 5 ../Zea_mays.AGPv2.14.dna/Zea_mays.AGPv2.14.dna.toplevel.fa
                   ../fastq/$fastq1 ../fastq/$fastq2 | samtools view -bSh - >$outfile.bam")
    #cmd2 <- paste0("samtools index sorted.$output.bam")
    #cmd3 <- paste0("samtools sort -m 10G -@ 2 $bam sorted.$output)
    cat(cmd1, file=shid, sep="\n", append=FALSE)
  }

  shcode <- paste("sh slurm-script/run_bwamem_$SLURM_ARRAY_TASK_ID.sh", sep="\n")

  set_array_job(shid="slurm-script/run_bwamem_array.sh",
                shcode=shcode, arrayjobs=arrayjobs,
                wd=NULL, jobid=jobid, email=email)
}

