#' \code{Run iget and iput to get/put data from/to iPlant.}
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
#'
#' @export
run_iput <- function(files, jobs1cpu=1,
                     localpwd="largedata/gatk_vcf",
                     email=NULL, runinfo = c(FALSE, "bigmemh", 1)){

  ####
  for(j in 1:ceiling(length(files)/jobs1cpu)){

    srow <- jobs1cpu*(j-1) + 1
    erow <- jobs1cpu*j

    cmd <- c()
    for(k in srow:erow){
      sh1 <- paste0("iput -X ", files[k])
      cmd <- c(cmd, sh1)
    }

    shid <- paste0(localpwd, "/iput_", j, ".sh")
    ### header of the shell code
    cat("### iput data from farm to iPlant",
        paste("###", format(Sys.time(), "%a %b %d %X %Y")),
        paste(""),
        cmd,
        file=shid, sep="\n", append=FALSE)
  }
  message("### run following:")
  message("### for example: sh input_j.sh")

  #shcode <- "sh slurm-script/iput_$SLURM_ARRAY_TASK_ID.sh"
  #set_array_job(shid="slurm-script/run_iput_array.sh",
  #              shcode=shcode, arrayjobs=paste0("1-", ceiling(length(files)/jobs1cpu)),
  #              wd=NULL, jobid="iput", email=email, runinfo=runinfo)
  #  sbatch -p bigmemh --mem 32784 --ntasks=4  slurm-script/run_gatk_array.sh

}
