#' \code{Run array job of R codes}
#'
#'
#' @param inputdf An input data.frame, with columns of file and out.
#' @param outdir The dir of shell files.
#' @param cmdno Number of commands to excute in each array.
#' @param rcodes The abosulte path of your R codes to run.
#' @param arrayshid The sbatch id.
#' @param email Your email address that farm will email to once the jobs were done/failed.
#' @param cmdno Number of commands per CPU, i.e. number of rows per inputdf.
#' @param runinfo Parameters specify the array job partition information.
#' A vector of c(FALSE, "bigmemh", "1"): 1) run or not, default=FALSE
#' 2) -p partition name, default=bigmemh and 3) --cpus, default=1.
#' It will pass to \code{set_array_job}.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#'
#'
#' @export
run_Rcodes <- function(
  inputdf, outdir, cmdno=100,
  rcodes = "lib/C_format.R",
  arrayshid = "slurm-script/run_bcf_query_array.sh",
  email=NULL, runinfo = c(FALSE, "bigmemh", 1)
){

  runinfo <- get_runinfo(runinfo)
  #### create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  dir.create(outdir, showWarnings = FALSE)

  tot <- ceiling(nrow(inputdf)/cmdno)
  for(j in 1:tot){

    shid <- paste0(outdir, "/run_rcode_", j, ".sh")

    ##chr:start-end
    #sh1 <- paste("cd", outdir)
    sh <- paste0('R --no-save --args ', j, ' < ', rcodes)

    cat(paste("### run Rcode", Sys.time(), sep=" "),
        sh,
        file=shid, sep="\n", append=FALSE)
  }

  shcode <- paste0("sh ", outdir, "/run_rcode_$SLURM_ARRAY_TASK_ID.sh")
  set_array_job(shid=arrayshid, shcode=shcode, arrayjobs=paste("1", tot, sep="-"),
                wd=NULL, jobid="rcode", email=email, runinfo=runinfo)
}

