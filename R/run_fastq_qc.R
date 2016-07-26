#' \code{Run fastq QC.}
#'
#' Run this 'seqtk fqchk -q 20' for all your fastq files.
#'
#' dependency:
#' seqtk version 1.0-r82-dirty
#'
#' @param df Data.frame with fq and out columns.
#' fq is the path of your fastq files and out is the path of your output file.
#' @param method Method for QC, "seqtk" or "FastQC" only.
#'
#' @param q Default=20. Note:use -q0 to get the distribution of all quality values
#' @param genomesize Maize genome size, default=2500000000.
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
#' df <- data.frame(fq=fqs, out=fqs)
#' df$out <- gsub("sra_.*", "qc", df$out)
#' run_fastq_qc(df, email=NULL, runinfo = c(FALSE, "bigmemh", 1))
#'
#' @export
run_fastq_qc <- function(df, method= "seqtk", q=20, email=NULL, runinfo = c(FALSE, "bigmemh", 1)){

  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)

  for(i in 1:nrow(df)){
    shid <- paste0("slurm-script/run_fqqc_", i, ".sh")

    if(method == "seqtk"){
      cmd <- paste0("seqtk fqchk -q 20 ", df$fq[i], " > ", df$out[i])
    }
    if(method == "FastQC"){
      cmd <- paste("fastqc --extract -f fastq", "-t", runinfo[3], df$fq[i])
    }
    cat(cmd, file=shid, sep="\n", append=FALSE)
  }

  if(method == "FastQC"){
    shcode <- c("module load fastqc/0.11.5", "sh slurm-script/run_fqqc_$SLURM_ARRAY_TASK_ID.sh")
  }else if(method == "seqtk"){
    shcode <- "sh slurm-script/run_fqqc_$SLURM_ARRAY_TASK_ID.sh"
  }
  set_array_job(shid="slurm-script/run_fqqc_array.sh",
                shcode=shcode, arrayjobs=paste("1", nrow(df), sep="-"),
                wd=NULL, jobid="fqQC", email=email, runinfo=runinfo)
}

#' @export
get_qc <- function(files, genomesize=2500000000){
  df <- data.frame()
  for(i in 1:length(files)){
    qc <- read.delim(files[i], skip = 2, header=FALSE)
    names(qc) <- c("pos", "bases", "A", "C", "G", "T", "N", "avgQ", "errQ", "plow", "phigh")
    #bases <- nrow(qc) -1
    tem <- data.frame(qc[1, ], bp=nrow(qc) -1, file=files[i])
    df <- rbind(df, tem)
  }
  df$depth <- df$bases/genomesize
  return(df)

}
