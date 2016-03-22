#' \code{Run fastq-dump.}
#'
#' fastq-dump to dump SRA file.
#'
#' see more detail about SRA with Aspera downloading:
#' \url{http://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.downloading_sra_data_using}
#'
#' @param filepath The absolute path of the SRA files.
#' @param slurmsh File name of the output shell command.
#' @param rmsra Remove the original SRA file after dumpping.
#' @param email Your email address that farm will email to once the job was done/failed.
#'
#' @return return a single shell script to run.
#'
#' @examples
#' run_fq_dump(filepath="/group/jrigrp4/BS_teo20/WGBS",
#' slurmsh="slurm-script/dump_WGBS.sh",
#' rmsra=TRUE, email=NULL)
#'
#' @export
run_fq_dump <- function(filepath="/group/jrigrp4/BS_teo20/WGBS",
                        slurmsh="slurm-script/dump_WGBS.sh",
                        rmsra=TRUE, email=NULL){

  files <- list.files(path=filepath, pattern="sra$")
  mysh <- paste("cd", filepath)
  for(i in 1:length(files)){
    out <- paste("fastq-dump --split-spot --split-3 -A", files[i])
    mysh <- c(mysh, out)
  }
  if(rmsra){
    mysh <- c(mysh, "rm *sra")
  }
  ### set up a single node job
  dir.create("slurm-script", showWarnings = FALSE)
  set_farm_job(slurmsh=slurmsh, shcode=mysh,
             wd=NULL, jobid="dump", email=email)
}
