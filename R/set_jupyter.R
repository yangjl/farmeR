#' \code{Set up Jupyter Nodebook}
#'
#'
#' @param port The port you are going to listen.
#'
#' @return nothing.
#'
#' @examples
#'
#' ## find ssh process ID and kill it in nodes
#' ps -ef | grep ssh | grep "jolyang"
#' kill -9 pid
#'
#' ## I don’t know what local port is available for me to use.
#' Linux: netstat -a | grep "999"
#' Mac OS X: lsof -i -P | grep LISTEN
#'
#' set_jupyter(port=9998)
#'
#' @export
set_jupyter <- function(port=9999){
  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  shid <- "slurm-script/run_jupyter.sh"

  cmd1 <- paste0("ssh -N -f -R ", port, ":localhost:", port, " $SLURM_SUBMIT_HOST")
  cmd2 <- paste0("jupyter notebook --port=", port)
  cat(c(cmd1, cmd2), file=shid, sep="\n", append=FALSE)
  message("###>>> in your local machine, run the following:")
  message(sprintf("ssh -f -N -L localhost:%s:localhost:%s farm", port, port))

  Sys.sleep(3)
  system("sh slurm-script/run_jupyter.sh")
}



