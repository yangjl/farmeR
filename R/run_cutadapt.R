#' \code{Run cutadapt}
#'
#' install in .local/bin using 'pip install --user --upgrade cutadapt'.
#'
#' For more information:
#' https://cutadapt.readthedocs.io/en/stable/index.html
#' version: v1.10
#'
#' @param inputdf Data.frame with fq1, fq2, out1 and out2 columns. Both input and output can be gzipped (.gz).
#' @param ad1 Adapter 3' method, see cutadapt documentation for more details. Same for ad2.
#'
#' @param q quality_cutoff Default=20. Note:use -q0 to get the distribution of all quality values
#'
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
#' inputdf <- data.frame(fq1=fqs[1], fq2=fqs[2], out1=fqs[1], out2=fqs[2])
#' df$out <- gsub("sra_.*", "qc", df$out)
#' run_fastq_qc(df, email=NULL, runinfo = c(FALSE, "bigmemh", 1))
#'
#' @export
run_cutadapt <- function(inputdf, ad1="AGATCGGAAGAGC", ad2="AGATCGGAAGAGC", q=20,
                         email=NULL, runinfo = c(FALSE, "bigmemh", 1)){

  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:nrow(inputdf)){
    shid <- paste0("slurm-script/run_cutadp_", i, ".sh")
    # AGATCGGAAGAGCGGTTCAGCAGGAATGCCG
    # AGATCGGAAGAGCACACGTCTGAACTCCAG
    cmd <- set_trim_PE(inputdf, i, ad1, ad2, quality_cutoff=q)
    cat(cmd, file=shid, sep="\n", append=FALSE)
  }

  set_array_job(shid="slurm-script/run_cutadp.sh",
                shcode="sh slurm-script/run_cutadp_$SLURM_ARRAY_TASK_ID.sh",
                arrayjobs=paste("1", nrow(inputdf), sep="-"),
                wd=NULL, jobid="cutadp", email=email, runinfo=runinfo)
}

#' @export
set_trim_PE <- function(inputdf, i, ad1, ad2, quality_cutoff){
  cmd <- paste("cutadapt -a", ad1, "-A", ad2,
               "-o", inputdf$out1[i], "-p", inputdf$out2[i],
               "-q", quality_cutoff,
               inputdf$fq1[i], inputdf$fq2[i])
  return(cmd)
}



