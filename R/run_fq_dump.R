#' \code{Run Aspera Connect to download from SRA.}
#'
#' Downloading SRA using 'ascp' utility or Aspera Connect.
#'
#' see more detail about SRA with Aspera downloading:
#' \url{http://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.downloading_sra_data_using}
#'
#' @param sra An input data.frame for SRA ids. Must contains column: SRR.
#' @param maxspeed The max speed for downloading.
#' @param outdir The output directory.
#' @param arrayjobs A character specify the number of array you try to run, i.e. 1-100.
#' @param jobid The job name show up in your sq NAME column.
#' @param email Your email address that farm will email to once the job was done/failed.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' sra <- data.frame(SRR=c("ERR877647", "ERR877648"),
#' SRX=c( "ERX957210", "ERX957211"),
#' pid=c( "1_Base1_Bbreve-sc-2188486", "P1_ECvsrS_1-sc-2201977"))
#' run_aspera(sra, maxspeed="200m", outdir=".", arrayjobs="1-2", jobid="aspera", email=NULL)
#'
#' @export
run_fq_dump <- function(pwd="/group/jrigrp4/BS_teo20/WGBS",
                        slurmsh="slurm-script/dump_WGBS.sh",
                        email=NULL){

  files <- list.files(path=pwd, pattern="sra$")
  mysh <- paste("cd", pwd)
  for(i in 1:length(files)){
    out1 <- paste("fastq-dump --split-spot --split-3 -A", files[i])
    mysh <- c(mysh, out)
  }


  set_farm_job(slurmsh=slurmsh,
             codesh=mysh,
             wd=NULL, jobid="dump", email=email)
}
