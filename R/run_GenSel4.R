#' \code{Run GenSel on farm}
#'
#' GenSel4R
#'
#' see more detail about gerpIBD
#' \url{https://github.com/yangjl/zmSNPtools}
#'
#'
#' @param inputdf An input data.frame.
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
#' inputdf <- data.frame(pi=0.995,
#'    geno="/Users/yangjl/Documents/GWAS2_KRN/SNP/merged/geno_chr.newbin",
#'    trainpheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
#'    testpheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
#'    chainLength=1000, burnin=100, varGenotypic=1.4, varResidual=2)
#'
#' run_GenSel4(inputdf, inpdir="slurm-script", email=NULL, runinfo = c(FALSE, "bigmemh", 1) )
#'
#' @export
run_GenSel4 <- function(
  inputdf, inpdir="largedata/",
  email=NULL, runinfo = c(FALSE, "bigmemh", 1)
){

  runinfo <- get_runinfo(runinfo)
  #### create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  dir.create(inpdir, showWarnings = FALSE)

  for(i in 1:nrow(inputdf)){
    shid <- paste0(inpdir, "/run_gensel_", i, ".inp")
    ### output the inp file:
    GS_cv_inp(
      inp= shid, pi=inputdf$pi[i], geno=inputdf$geno[i],
      trainpheno=inputdf$trainpheno[i], testpheno=inputdf$testpheno[i],
      chainLength=inputdf$chainLength[i], burnin=inputdf$burnin[i],
      varGenotypic=inputdf$varGenotypic[i], varResidual=inputdf$varResidual[i]
    )
  }

  shcode <- paste0("GenSel4R ", inpdir, "/run_gensel", "_$SLURM_ARRAY_TASK_ID.inp > ",
                   inpdir, "/run_gensel", "_$SLURM_ARRAY_TASK_ID.log")
  set_array_job(shid="slurm-script/run_gensel_array.sh",
                shcode=shcode, arrayjobs=paste("1", nrow(inputdf), sep="-"),
                wd=NULL, jobid="gensel", email=email, runinfo=runinfo)
}

############################ GenSel for cross-validation
GS_cv_inp <- function(
  inp, pi,geno, trainpheno,
  testpheno, chainLength, burnin, varGenotypic, varResidual
){

  cat(paste("// gensel input file written", Sys.time(), sep=" "),

      "analysisType Bayes",
      "bayesType BayesC",
      paste("chainLength", chainLength, sep=" "),
      paste("burnin", burnin=burnin, sep=" "),
      paste("probFixed", pi, sep=" "),

      paste("varGenotypic",  varGenotypic, sep=" "),
      paste("varResidual",  varResidual, sep=" "),
      "nuRes 10",
      "degreesFreedomEffectVar 4",
      "outputFreq 100",
      "seed 1234",
      "mcmcSamples yes",
      "plotPosteriors no",
      "FindScale no",
      "modelSequence no",
      "isCategorical no",

      "",
      "// markerFileName",
      paste("markerFileName", geno, sep=" "),
      "",
      "// train phenotypeFileName",
      paste("trainPhenotypeFileName", trainpheno, sep=" "),
      "",
      "// test phenotypeFileName",
      paste("testPhenotypeFileName", testpheno, sep=" "),

      #"// includeFileName",
      #paste("includeFileName", inmarker, sep=" "),

      file=inp, sep="\n", append=FALSE
  )
}

