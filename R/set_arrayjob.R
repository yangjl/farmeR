<<<<<<< HEAD
#' Set up array job on farm
#'
#' 'Farm' is a computer cluster running slurm system.
#' Note, 'bigmem' mem=8000/cpu, "hi/med/low' mem=25000/cpu, "serial" mem=1500/cpu.
#'
#' @param shid Relative or absolute path and file name of your shell code, i.e. 'largedata/GenSel/CL_test.sh'.
#' @param shcode The commands inside your sh file.
#' @param arrayjobs A character specify the number of array you try to run, i.e. '1-100'.
#' @param wd Working directory, default=NULL => using your current directory.
#' @param jobid The job name show up in your sq "NAME' column.
#' @param email Your email address that farm will email to once the job was done/failed.
#'
#' @return a shell file.
#' @examples
#' set_array_job <- function(shid="largedata/GenSel/CL_test.sh",
#' shcode="sh largedata/myscript.sh", arrayjobs="1-700",wd=NULL, jobid="myjob", email=NULL)
#'
=======
#' \code{Set up array job on farm}
#'
#' Farm is a computer cluster running slurm system.
#' Note, bigmem mem=8000 per cpu, hi med low mem=25000 per cpu, serial mem=1500 per cpu.
#'
#' @param shid Relative or absolute path and file name of your shell code, i.e. CL_test.sh.
#' @param shcode The commands inside your sh file.
#' @param arrayjobs A character specify the number of array you try to run, i.e. 1-100.
#' @param wd Working directory, default=NULL. It will use your current directory.
#' @param jobid The job name show up in your sq NAME column.
#' @param email Your email address that farm will email to once the job was done or failed.
#'
#' @return return a shell file.
#'
#' @examples
#' for(i in 1:10){
#'     shid <- paste0("slurm-script/run_", i, ".sh")
#'     command <- paste0("bedtools getfasta -name -tab -fi roast.chrom.", i, ".msa.in")
#'     cat(command, file=shid, sep="\n", append=FALSE)
#' }
#' shcode <- paste("module load bismark/0.14.3", "sh slurm-script/run_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
#'
#'
#' @export

set_array_job <- function(shid="largedata/GenSel/CL_test.sh",
                         shcode="sh largedata/myscript.sh",
                         arrayjobs="1-700",
                         wd=NULL, jobid="myjob", email=NULL){

    #message(sprintf("###>>> cp from Introgression, tailored for pvpDiallel"))

    ##### setup working directory
    if(is.null(wd)){
       wd <- getwd()
    }
    dir.create("slurm-log", showWarnings = FALSE)
    sbath <- paste0(wd, "/slurm-log/")
    sbatho <- paste0(sbath, "testout-%j.txt")
    sbathe <- paste0(sbath, "err-%j.txt")

    #### parameters pass to slurm script
    cat(paste("#!/bin/bash -l"),
        #-D sets your project directory.
        #-o sets where standard output (of your batch script) goes.
        #-e sets where standard error (of your batch script) goes.
        #-J sets the job name.
        paste("#SBATCH -D", wd, sep=" "),
        paste("#SBATCH -o", sbatho, sep=" "),
        paste("#SBATCH -e", sbathe, sep=" "),
        paste("#SBATCH -J", jobid, sep=" "),
        paste0("#SBATCH --array=", arrayjobs),
        paste0("#SBATCH --mail-user=", email),
        paste("#SBATCH --mail-type=END"),
        paste("#SBATCH --mail-type=FAIL #email if fails"),


        "set -e",
        "set -u",
        "",
        #"module load gmap/2014-05-15",
        file=shid, sep="\n", append=FALSE);

    #### attach some sh scripts
    cat(shcode, file=shid, sep="\n", append=TRUE)
    message(paste("###>>> In this path: cd ", wd, sep=""), "\n",
<<<<<<< HEAD
            paste("###>>> RUN: sbatch -p bigmemh", shid),
            "")

}

#' Set up one farm job
#'
#' 'Farm' is a computer cluster running slurm system.
#' Note, 'bigmem' mem=8000/cpu, "hi/med/low' mem=25000/cpu, "serial" mem=1500/cpu.
#'
#' @param slurmsh Relative or absolute path and file name of your shell code, i.e. 'largedata/GenSel/CL_test.sh'.
#' @param shcode The commands inside your sh file.
#' @param wd Working directory, default=NULL => using your current directory.
#' @param jobid The job name show up in your sq "NAME' column.
#' @param email Your email address that farm will email to once the job was done/failed.
#'
#' @return a shell file.
#' @examples
#'
#'
=======
            paste("###>>> RUN: sbatch -p bigmemh", shid))

}

#' \code{Set up one farm job}
#'
#' Farm is a computer cluster running slurm system.
#' Note, bigmem mem=8000/cpu, hi/med/low mem=25000/cpu, serial mem=1500/cpu.
#'
#' @param slurmsh Relative or absolute path and file name of your shell code, i.e. largedata/GenSel/CL_test.sh.
#' @param shcode The commands inside your sh file.
#' @param wd Working directory, default=NULL. It will use your current directory.
#' @param jobid The job name show up in your sq NAME column.
#' @param email Your email address that farm will email to once the job was done/failed.
#'
#' @return a shell file.
#'
#' @examples
#'
#' @export
>>>>>>> e2dfb9d60d5ce3dc09c5a85b3d3b2193373d0989
set_farm_job <- function(slurmsh="largedata/GenSel/CL_test.sh",
                         shcode="sh largedata/myscript.sh",
                         wd=NULL, jobid="myjob", email=NULL){
  ##### setup working directory
  if(is.null(wd)){
    wd <- getwd()
  }
  dir.create("slurm-log", showWarnings = FALSE)
  sbath <- paste0(wd, "/slurm-log/")
  sbatho <- paste0(sbath, "testout-%j.txt")
  sbathe <- paste0(sbath, "err-%j.txt")
  sbathJ <- jobid

  #### parameters pass to slurm script
  cat(paste("#!/bin/bash"),
      #-D sets your project directory.
      #-o sets where standard output (of your batch script) goes.
      #-e sets where standard error (of your batch script) goes.
      #-J sets the job name.
      paste("#SBATCH -D", wd, sep=" "),
      paste("#SBATCH -o", sbatho, sep=" "),
      paste("#SBATCH -e", sbathe, sep=" "),
      paste("#SBATCH -J", sbathJ, sep=" "),
      paste0("#SBATCH --mail-user=", email),
      paste("#SBATCH --mail-type=END"),
      paste("#SBATCH --mail-type=FAIL #email if fails"),


      "set -e",
      "set -u",
      "",
      #"module load gmap/2014-05-15",
      file=slurmsh, sep="\n", append=FALSE);

  #### attach some sh scripts
  cat(shcode, file=slurmsh, sep="\n", append=TRUE)

  message(paste("###>>> In this path: cd ", wd, sep=""), "\n",
<<<<<<< HEAD
          paste("###>>> RUN: sbatch -p bigmemh --ntasks=1", slurmsh),
          "")
=======
          paste("###>>> RUN: sbatch -p bigmemh --ntasks=1", slurmsh))
>>>>>>> e2dfb9d60d5ce3dc09c5a85b3d3b2193373d0989

}
