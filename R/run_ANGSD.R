#' \code{Run ANGSD on farm}
#'
#' angsd
#'
#' see more detail about ANGSD
#' \url{http://popgen.dk/angsd/index.php/Main_Page}
#'
#' Scripts learned from Tim's paper \url{https://github.com/timbeissinger/Maize-Teo-Scripts}
#'
#'
#' @param inputdf An input data.frame.
#' @param minQ Minimum base quality score.
#' @param minMapQ [int]=0 Minimum mapQ quality.
#'
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
#' traits <- tolower(c("ASI", "DTP", "DTS", "EHT", "GY", "PHT", "TW"))
#' inputdf <- data.frame(pi=0.995,
#'    geno="/Users/yangjl/Documents/GWAS2_KRN/SNP/merged/geno_chr.newbin",
#'    trainpheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
#'    testpheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
#'    chainLength=1000, burnin=100, varGenotypic=1.4, varResidual=2,
#'    out="test_out")
#'
#' run_GenSel4(inputdf, inpdir="slurm-script", email=NULL, runinfo = c(FALSE, "bigmemh", 1) )
#'
#' @export
run_ANGSD <- function(
  inputdf, inpdir="largedata/", cmdno=1,
  shid = "slurm-script/run_gensel_array.sh",
  email=NULL, runinfo = c(FALSE, "bigmemh", 1)
){

  runinfo <- get_runinfo(runinfo)
  #### create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  dir.create(inpdir, showWarnings = FALSE)

  for(i in 1:nrow(inputdf)){
    inpid <- paste0(inpdir, "/", inputdf$out[i], ".inp")
    ### output the inp file:
    GS_cv_inp(
      inp= inpid, pi=inputdf$pi[i], geno=inputdf$geno[i],
      trainpheno=inputdf$trainpheno[i], testpheno=inputdf$testpheno[i],
      chainLength=inputdf$chainLength[i], burnin=inputdf$burnin[i],
      varGenotypic=inputdf$varGenotypic[i], varResidual=inputdf$varResidual[i]
    )
  }

  ### setup shell id
  set_GS(inputdf, cmdno, inpdir)

  shcode <- paste0("sh ", inpdir, "/run_gensel_$SLURM_ARRAY_TASK_ID.sh")
  set_array_job(shid=shid, shcode=shcode, arrayjobs=paste("1", nrow(inputdf)/cmdno, sep="-"),
                wd=NULL, jobid="gensel", email=email, runinfo=runinfo)
}

############################ GenSel for cross-validation
set_angsd <- function(
  shfile,
  bamlist, outfile, ref, cpu, nInd,
  minMapQ, minQ,
  glikehood, indF, anc

){

  cat(paste("### ANGSD input file written at", Sys.time(), sep=" "),

      # Now redo angsd sample allele frequency calculation by conditioning
      # on the sites that occur in both populations.

      ### input bam file arguments:
      paste("angsd -bam", bamlist),
      paste("-out", outfile),
      paste("-uniqueOnly 0 -minMapQ", minMapQ, "-minQ", minQ,  "-baq 1", "-ref", ref),

      ### get MAF for fixed major and minor alleles
      paste("-GL", glikehood, "-P", cpu),
      paste("-doMajorMinor 1 -doMaf 1 -doSaf 2 "),
      paste("-indF", indF, "-anc", anc),

      ### filters
      paste("-nInd", nInd, "-minInd", minInd),


      file=shfile, sep=" ", append=FALSE
  )
}

