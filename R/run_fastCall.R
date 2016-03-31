#' \code{Run fastCall on farm}
#'
#' Super fast variant caller for whole genome shotgun (WGS) sequencing data. It works for diploid species,
#' including both inbreds and outcrossers.
#'
#' see more details:
#' \url{https://github.com/Fei-Lu/FastCall}
#'
#'
#' Prerequisites:
#' module load java/1.8
#' Samtools
#'
#'
#' @param ref.fa The full path of genome with bwa indexed reference fasta file.
#' @param fastCallpwd The absolute path of GenomeAnalysisTK.jar.
#' @param bamdir The full path of the bam files.
#' @param baminfofile Taxa info for bam files, each line contains 1 to N bam files for the taxa.
#' @param chr Chr number, i.e. 1, default=NULL, ten chrs will run as 10 array jobs.
#' @param email Your email address that farm will email to once the jobs were done/failed.
#' @param outdir The full path of the output files. Note log and shell codes will also put here.
#' @param runinfo Parameters specify the array job partition information.
#' A vector of c(FALSE, "bigmemh", "1"): 1) run or not, default=FALSE
#' 2) -p partition name, default=bigmemh and 3) --cpus, default=1.
#' It will pass to \code{set_array_job}.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' bam <- data.frame(taxa=c("t1", "t2"), bam=c("t1_1.bam", "t2.bam"), bam=c("t1_2.bam", ""))
#' write.table(bam, "test/taxaBamMap.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE )
#'
#' run_fastCall(ref.fa = "~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
#'              fastCallpwd = "$HOME/bin/fastCall", bamdir = "$HOME/test",
#'              baminfofile = "test/taxaBamMap.txt", chr = NULL, outdir = "test",
#'              email = NULL,
#'              runinfo = c(FALSE, "bigmemh", 1))
#'
#' @export
run_fastCall <- function(
  ref.fa="~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
  fastCallpwd="$HOME/bin/fastCall",
  bamdir="/", baminfofile="taxaBamMap.txt", chr=NULL,
  outdir="",
  email=NULL, runinfo = c(FALSE, "bigmemh", 1)   ){


  ### determine memory based on partition
  runinfo <- get_runinfo(runinfo)

  #### create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)

  if(is.null(chr)) {
    chr <- 1:10
    arrayjobs <- "1-10"
  }else{
    arrayjobs <- as.character(chr)
  }

  for(i in chr){
    parid <- paste0(outdir, "/parameter_chr", i, ".txt")
    cat("### Parameter1: Reference genome file with an index file (.fai).",
        "### The reference should be in FastA format. Chromosomes are labled as 1-based numbers (1,2,3,4,5...).",
        ref.fa,
        "### Parameter2: Bam directory, where your bam files are",
        bamdir,
        "### Parameter3: Taxa bam information file",
        "### including the info about what bams are included for each taxon",
        baminfofile,
        "#Parameter4:	Chromosome on which genotyping will be performed (e.g. 1)",
        i,
        "#Parameter5:	VCF output directory",
        outdir,
        file=parid, sep="\n", append=FALSE)

    shid <- paste0(outdir, "/run_fastCall_chr", i, ".sh")
    cat(paste0("cd ", fastCallpwd),
        paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g -jar FastCall.jar ",
               parid,
               " > ", outdir, "/log_chr", i, ".txt"),
        file=shid, sep="\n", append=FALSE)
  }

  cmd1 <- "module load java"
  cmd2 <- paste0("sh ", outdir, "/run_fastCall_$SLURM_ARRAY_TASK_ID.sh")
  shcode <- c(cmd1, cmd2)
  set_array_job(shid="slurm-script/run_fastCall_array.sh",
                shcode=shcode, arrayjobs=arrayjobs,
                wd=NULL, jobid="fastCall", email=email, runinfo=runinfo)
  #  sbatch -p bigmemh --mem 32784 --ntasks=4  slurm-script/run_gatk_array.sh
}




