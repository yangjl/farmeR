#' \code{Run bs2vcf job on farm}
#'
#' convert BSMAP format to VCF format
#' Find help about 'bs2vcf -h'.
#' bs2vcf -p /group/jrigrp4/BS_teo20/BSMAP_output -i test_methratio.txt -o test_methratio.vcf
#'
#' @param inputdf An input data.frame for fastq files. Must contains bsmap and out.
#'
#' @param email Your email address that farm will email to once the jobs were done/failed.
#' @param runinfo Parameters specify the array job partition information.
#' A vector of c(FALSE, "bigmemh", "1"): 1) run or not, default=FALSE
#' 2) -p partition name, default=bigmemh and 3) --cpus, default=1.
#' It will pass to \code{set_array_job}.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' inputdf <- data.frame(bsmap="test_ratio.txt", out="test_ratio.vcf")
#' run_bs2vcf( inputdf, email=NULL, runinfo = c(FALSE, "bigmemh", 1))
#'
#' @export
run_bs2vcf <- function(
  inputdf,
  email=NULL, runinfo = c(FALSE, "bigmemh", 1)
){

  # create dir if not exist
  runinfo <- get_runinfo(runinfo)
  dir.create("slurm-script", showWarnings = FALSE)

  for(i in 1:nrow(inputdf)){
    shid <- paste0("slurm-script/run_bs2vcf_", i, ".sh")
    ### header of the shell code
    cmd1 <- paste("bs2vcf -i", inputdf$bsmap[i], "-o", inputdf$out[i])

    #bgzip test_methratio.vcf; tabix -p vcf test_methratio.vcf.gz
    cmd2 <- paste0("bgzip ", inputdf$out[i], "; tabix -p vcf ", inputdf$out[i], ".gz")

    cat("### ba2vcf pipeline created by farmeR",
        paste("###", format(Sys.time(), "%a %b %d %X %Y")),
        paste(""),
        cmd1,
        cmd2,
        file=shid, sep="\n", append=FALSE)
  }

  shcode <- paste("sh slurm-script/run_bs2vcf_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  set_array_job(shid="slurm-script/run_bs2vcf_array.sh",
                shcode=shcode, arrayjobs=paste("1", nrow(inputdf), sep="-"),
                wd=NULL, jobid="bs2vcf", email=email, runinfo=runinfo)
}
