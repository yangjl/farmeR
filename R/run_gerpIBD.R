#' \code{Run gerpIBD on farm}
#'
#' gerpIBD: compute sum of GERP score in a given IBD region across genotypes.
#'
#' see more detail about gerpIBD
#' \url{https://github.com/yangjl/zmSNPtools}
#'
#'
#' @param inputdf An input data.frame. Type "gerpIBD -h" to understand the meaning of the input parameters.
#' @param email Your email address that farm will email to once the jobs were done/failed.
#' @param runinfo Parameters specify the array job partition information.
#' A vector of c(FALSE, "bigmemh", "1"): 1) run or not, default=FALSE
#' 2) -p partition name, default=bigmemh and 3) --cpus, default=1.
#' It will pass to \code{set_array_job}.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' traits <- tolower(c("ASI", "DTP", "DTS", "EHT", "GY", "PHT", "TW"))
#' inputdf <- data.frame(trait=rep(traits, each=44),
#'    mode=rep(rep(c("add", "dom", "k", "k5"), each=11), times=7 ),
#'    nperm=rep(0:10, times=28),
#'    d="largedata/IBD/allsnps_11m_IBD.bed",
#'    s="largedata/SNP/allsnps_11m.dsf5",
#'    g="largedata/SNP/geno_b0_cs/gerpv2_b0_cs0.csv",
#'    f="largedata/snpeff/perse/asi_k.txt",
#'    out="res")
#'
#' run_gerpIBD(inputdf[1:10,], email=NULL, runinfo = c(FALSE, "bigmemh", 1) )
#'
#' @export
run_gerpIBD <- function(inputdf, email=NULL, runinfo = c(FALSE, "bigmemh", 1) ){
  runinfo <- get_runinfo(runinfo)
  #### create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:nrow(inputdf)){

    shid <- paste0("slurm-script/run_gerpibd_", i, ".sh")
    ### header of the shell code
    cat("### gerpIBD pipeline created by farmeR",
        paste("###", format(Sys.time(), "%a %b %d %X %Y")),
        paste(""),
        paste("gerpIBD -d", inputdf$d[i], "-s", inputdf$s[i], "-g", inputdf$g[i], "-f", inputdf$f[i],
              "-o", inputdf$out[i], "-l", inputdf$l[i], "-t", inputdf$t[i]),
        file=shid, sep="\n", append=FALSE)
  }

  shcode <- paste("sh slurm-script/run_gerpibd_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  set_array_job(shid="slurm-script/run_gerpibd_array.sh",
                shcode=shcode, arrayjobs=paste("1", nrow(inputdf), sep="-"),
                wd=NULL, jobid="gerpibd", email=email, runinfo=runinfo)
}
