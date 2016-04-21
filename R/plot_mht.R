#' \code{manhattan plot}
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
#' run_Rcodes(inputdf=data.frame(file=1:11, out=10), outdir="slurm-script", cmdno=10,
#'            rcodes = "lib/C_format.R", arrayshid = "slurm-script/run_rcode_array.sh",
#'            email=NULL, runinfo = c(FALSE, "bigmemh", 1))
#'
#' @export
plot_mht <- function(res=res, cex=.9, pch=16, col=rep(c("slateblue", "cyan4"), 5),
                         GAP=5e+06, yaxis=NULL,
                         col2plot="ModelFreq", ... ){

  res <- newpos(res, GAP = GAP)
  chrtick <- chrline_tick(GAP = GAP)

  #### setup the cavon
  if(is.null(yaxis)){
    plot(x=-1000, y=-1000,  type="p", xaxt="n", xlab="",
         xlim=c(0, max(chrtick$chrlines)), ylim=c(0, max(res[, col2plot], na.rm=TRUE)*1.3 ),
         ...)
  }else{
    plot(x=-1000, y=-1000,  type="p", xaxt="n", yaxt="n", xlab="",
         xlim=c(0, max(chrtick$chrlines)),
         ...)
    axis(side=2, at=yaxis, labels=yaxis)
  }
  axis(side=1, at=chrtick$ticks, labels=c("chr1", "chr2", "chr3", "chr4", "chr5",
                                          "chr6", "chr7", "chr8", "chr9", "chr10"))
  abline(v=chrtick$chrlines, col="grey")

  for(i in 1:10){
    points(x=subset(res, chr==i)$newpos, y=res[res$chr==i, col2plot],
           pch = pch, col=col[i], cex=cex);
  }
}


newpos <- function (dataframe, GAP = 5e+06, version = "v3")
{
  d <- dataframe
  if (!("chr" %in% names(d) & "pos" %in% names(d))){
    stop("Make sure your data frame contains columns chr and pos")
  }

  if(version == "v2"){
    cl <- read.csv("~/Documents/Github/zmSNPtools/Rcodes/chr_length_B73v2.csv")
  }
  if(version == "v3"){
    cl <- read.csv("~/Documents/Github/zmSNPtools/Rcodes/chr_length_B73v3.csv")
  }

  cl$accumpos <- cl$BP
  cl <- cl[order(cl$CHR), ]
  d$newpos <- d$pos;
  for (i in 2:10) {
    cl[cl$CHR == i, ]$accumpos <- cl[cl$CHR == (i - 1), ]$accumpos + cl[cl$CHR == i, ]$accumpos + GAP
    d[d$chr == i, ]$newpos <- d[d$chr == i, ]$pos + cl[cl$CHR == (i - 1), ]$accumpos + GAP
  }
  return(d)
}

chrline_tick <- function(GAP=5e+06, version = "v3"){
  #xscale:
  if(version == "v2"){
    cl <- read.csv("~/Documents/Github/zmSNPtools/Rcodes/chr_length_B73v2.csv")
  }
  if(version == "v3"){
    cl <- read.csv("~/Documents/Github/zmSNPtools/Rcodes/chr_length_B73v3.csv")
  }

  names(cl) <- c("chr", "snp", "pos")
  cl <- newpos(cl, GAP=GAP)

  cl$ticks <- cl$pos[1]/2
  cl$chrlines <- cl$pos[1]+GAP/2
  for(i in 2:10){
    cl$ticks[i] <- cl$newpos[i-1] + (cl$newpos[i]-cl$newpos[i-1])/2;
    cl$chrlines[i] <- cl$newpos[i]+ GAP/2;
  }
  return(cl)
}
