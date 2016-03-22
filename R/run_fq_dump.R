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
#' @param gzip GZIP the fastq files.
#' @param run  Parameters control the array job partition.
#' A vector of c(TRUE, "bigmemh", "8196", "1"): 1) run or not, 2) -p partition name, 3) --mem, adn 4) --ntasks.
#'
#' @return return a single shell script to run.
#'
#' @examples
#' ## run a single node job:
#' run_fq_dump(filepath="/group/jrigrp4/BS_teo20/WGBS",
#'             slurmsh="slurm-script/dump_WGBS.sh", rmsra=TRUE, email=NULL)
#'
#' ##  run array job:
#' run_fq_dump2(filepath="test", rmsra=TRUE, gzip=TRUE, email=NULL, run=c(TRUE, "bigmemh", "8196", "1"))
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

#' @rdname run_fq_dump
#' @export
run_fq_dump2 <- function(filepath="/group/jrigrp4/BS_teo20/WGBS",
                         rmsra=TRUE, gzip=TRUE, email=NULL,
                         run=c(TRUE, "bigmemh", "8196", "1")){

  files <- list.files(path=filepath, pattern="sra$")
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:length(files)){

    shid <- paste0("slurm-script/run_dump_", i, ".sh")
    cmd1 <- paste0("cd ", filepath)
    cmd2 <- paste0("fastq-dump --split-spot --split-3 -A ", files[i])
    cmd <- c(cmd1, cmd2)
    if(rmsra){
      cmd3 <- paste0("rm ", files[i])
      cmd <- c(cmd, cmd3)
    }
    if(gzip){
      cmd4 <- paste0("gzip ", paste0(files[i], "_1.fastq"))
      cmd5 <- paste0("gzip ", paste0(files[i], "_2.fastq"))
      cmd <- c(cmd, cmd4, cmd5)
    }
    cat(cmd, file=shid, sep="\n", append=FALSE)
  }

  shcode <- paste("sh slurm-script/run_dump_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  myshid <- "slurm-script/run_dump_array.sh"

  set_array_job(shid=myshid,
                shcode=shcode, arrayjobs=paste("1", length(files), sep="-"),
                wd=NULL, jobid="dump", email=email)

  if(run[1]){
    runcode <- paste0("sbatch -p ", run[2], " --mem ", run[3], " --ntasks=", run[4], " ", myshid)
    print(runcode)
    system(runcode)
  }
}
