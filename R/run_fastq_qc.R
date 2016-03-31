#' \code{Run fastq QC.}
#'
#' Run this 'seqtk fqchk -q 20' for all your fastq files.
#'
#' dependency:
#' seqtk version 1.0-r82-dirty
#'
#' @param df Data.frame with fq and out columns.
#' fq is the path of your fastq files and out is the path of your output file.
#' @param email Your email address that farm will email to once the job was done/failed.
#' @param runinfo Parameters specify the array job partition information.
#' A vector of c(FALSE, "bigmemh", "1"): 1) run or not, default=FALSE
#' 2) -p partition name, default=bigmemh and 3) --cpus, default=1.
#' It will pass to \code{set_array_job}.
#'
#' @return return A batch of shell scripts.
#'
#' @examples
#' fqs <- c("$HOME/dbcenter/Ecoli/fastq/SRR2921970.sra_1.fastq.gz",
#'          "$HOME/dbcenter/Ecoli/fastq/SRR2921970.sra_2.fastq.gz")
#' df <- data.frame(fq=fqs, out=fqs)
#' df$out <- gsub("sra_.*", "qc", df$out)
#' run_fastq_qc(df, email=NULL, runinfo = c(FALSE, "bigmemh", 1))
#'
#' @export
run_fastq_qc <- function(df, email=NULL, runinfo = c(FALSE, "bigmemh", 1)){

  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:nrow(df)){
    shid <- paste0("slurm-script/run_fqqc_", i, ".sh")
    cmd <- paste0("seqtk fqchk -q 20 ", df$fq[i], " > ", df$out[i])
    cat(cmd, file=shid, sep="\n", append=FALSE)
  }

  shcode <- paste("sh slurm-script/run_fqqc_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  set_array_job(shid="slurm-script/run_fqqc_array.sh",
                shcode=shcode, arrayjobs=paste("1", nrow(df), sep="-"),
                wd=NULL, jobid="fqQC", email=email, runinfo=runinfo)
}
