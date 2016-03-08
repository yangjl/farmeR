#' \code{Run array Bismark job on farm}
#'
#' Run bowtie2/2.2.5
#' Run bismark/0.14.3
#'
#' Allow one mismatch 'n 1'
#'
#'
#' see more detail about Bismark:
#' \url{http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf}
#'
#' @param input_df An input data.frame object. Must contains fq1, fq2 and out.
#' @param genome The folder of genome prepared by bismark.
#' @param cpu Number of CPU to use per job.
#' @param outdir Folder for output.
#' @param arrayjobs A character specify the number of array you try to run, i.e. 1-100.
#' @param jobid The job name show up in your sq NAME column.
#' @param email Your email address that farm will email to once the job was done/failed.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' input_df <- data.frame(fq1=c("f_1.fq", "t_1.fq"), fq2=c("f_1.fq", "t_2.fq"), out=c("t1", "t2"))
#' runa_bismark(input_df, genome="/home/jolyang/dbcenter/AGP/AGPv2",
#' cpu=4, outdir="/group/jrigrp4/BS_teo20/WGBS/BSM", arrayjobs="1-5", jobid="bs1-5",
#' email=NULL)
#' @export
runa_bismark <- function(input_df,
                         genome="/home/jolyang/dbcenter/AGP/AGPv2",
                         cpu=4,
                         outdir="/group/jrigrp4/BS_teo20/WGBS/BSM",
                         arrayjobs="1-5",
                         jobid="bs1-5",
                         email="yangjl0930@gmail.com"){

  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)


  for(i in 1:nrow(input_df)){

    shid <- paste0("slurm-script/run_bismark_", i, ".sh")
    #out <- gsub(".*/", "", out)
    #outfile <- paste0(outdir, "/", out)
    cmd <- paste("bismark --bowtie2 -n 1", genome, "-p", cpu,
                 "-1", input_df$fq1[i],  "-2", input_df$fq2[i],
                 "--output_dir", outdir,  "--basename", input_df$out[i])

    cat(cmd, file=shid, sep="\n", append=FALSE)
  }

  shcode <- paste("module load bismark/0.14.3", "module load bowtie2/2.2.5",
                  "sh slurm-script/run_bismark_$SLURM_ARRAY_TASK_ID.sh", sep="\n")

  set_array_job(shid="slurm-script/run_bismark_array.sh",
                shcode=shcode, arrayjobs=arrayjobs,
                wd=NULL, jobid=jobid, email=email)
}

