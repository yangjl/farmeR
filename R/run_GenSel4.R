#' \code{Run GenSel on farm}
#'
#' GenSel4R
#'
#' see more detail about gerpIBD
#' \url{https://github.com/yangjl/zmSNPtools}
#'
#'
#' @param inputdf An input data.frame.
#' @param cv Cross-validation experiment [TRUE/FALSE, default=FALSE].
#'           If TRUE, inputdf must contain `trainpheno` and `testpheno`; if FALSE, only need `pheno`.
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
run_GenSel4 <- function(
  inputdf, cv=FALSE, inpdir="largedata/", cmdno=1,
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

    if(cv){
      GS_cv_inp(
        inp= inpid, pi=inputdf$pi[i], geno=inputdf$geno[i],
        trainpheno=inputdf$trainpheno[i], testpheno=inputdf$testpheno[i],
        chainLength=inputdf$chainLength[i], burnin=inputdf$burnin[i],
        varGenotypic=inputdf$varGenotypic[i], varResidual=inputdf$varResidual[i]
      )
    }else{
      GS_regular_inp(
        inp= inpid, pi=inputdf$pi[i], geno=inputdf$geno[i],
        pheno=inputdf$pheno[i],
        chainLength=inputdf$chainLength[i], burnin=inputdf$burnin[i],
        varGenotypic=inputdf$varGenotypic[i], varResidual=inputdf$varResidual[i]
      )
    }

  }

  ### setup shell id
  set_GS(inputdf, cmdno, inpdir)

  shcode <- paste0("sh ", inpdir, "/run_gensel_$SLURM_ARRAY_TASK_ID.sh")
  set_array_job(shid=shid, shcode=shcode, arrayjobs=paste("1", nrow(inputdf)/cmdno, sep="-"),
                wd=NULL, jobid="gensel", email=email, runinfo=runinfo)
}

############################ GenSel for cross-validation
GS_regular_inp <- function(
  inp, pi,geno, pheno,
  chainLength, burnin, varGenotypic, varResidual
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
      paste("phenotypeFileName", pheno, sep=" "),

      #"// includeFileName",
      #paste("includeFileName", inmarker, sep=" "),

      file=inp, sep="\n", append=FALSE
  )
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

set_GS <- function(inputdf, cmdno, inpdir){
  ####
  for(j in 1:(nrow(inputdf)/cmdno)){

    srow <- cmdno*(j-1) + 1
    erow <- cmdno*j
    shid <- paste0(inpdir, "/run_gensel_", j, ".sh")

    sh1 <- paste0("GenSel4R ", inpdir, "/", inputdf$out[srow:erow], ".inp > ",
                     inpdir, "/", inputdf$out[srow:erow], ".log")
    sh2 <- c(paste0("rm ", inpdir, "/", inputdf$out[srow:erow], ".mcmcSamples1"),
             paste0("rm ", inpdir, "/", inputdf$out[srow:erow], ".mrkRes1"),
             paste0("rm ", inpdir, "/", inputdf$out[srow:erow], ".inp"),
             paste0("rm ", inpdir, "/", inputdf$out[srow:erow], ".out1"),
             paste0("rm ", inpdir, "/", inputdf$out[srow:erow], ".cgrRes1"))
    cat(paste("### run GenSel", Sys.time(), sep=" "),
        c(sh1, sh2),
        file=shid, sep="\n", append=FALSE)
  }
}
