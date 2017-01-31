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
#' @param r Chr regions.
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
#' run_ANGSD(type="GL", shfile="slurm-script/run_angsd_cmd.sh",
#' bamlist="list.txt", outfile="out", ref="ref.fa", nInd=20, minInd=10,
#' indF="indF.txt", anc="ref.fa",
#' email=NULL, runinfo = c(FALSE, "bigmemh", 1)
#'
#'
#' @export
run_ANGSD <- function(
  type="GL",
  shfile="slurm-script/run_angsd_cmd.sh",
  bamlist="list.txt", outfile="out", ref="ref.fa", nInd=20, minInd=10,
  indF="indF.txt", anc="ref.fa", r,
  email=NULL, runinfo = c(FALSE, "bigmemh", 1)
){

  runinfo <- get_runinfo(runinfo)
  #### create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)


  ### set the first several lines of the bash file
  cat("#!/usr/bin/env bash",
      "set -e",
      "set -o pipefail",
      "",
      paste("### codes wrote by farmeR at", Sys.time(), sep=" "),
      file=shfile, sep="\n", append=FALSE)

  ### get gentoype likelihoods
  if(type=="GL"){
    #$ANGSD/angsd -doGlf 3 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-4 -out testF.HWE
    set_angsd_gl(
      shfile, r,
      bamlist, outfile, cpu=runinfo[3],
      ref, anc,
      minMapQ=30,minQ=20,
      glikehood=1,
      doMajorMinor=1, doMaf=1,
      doGlf=3, SNP_pval=2e-6
    )
    message(sprintf("### get gentoype likelihoods!"))
  }

  if(type == "indF"){
    set_angsd_gl(
      shfile, r,
      bamlist, outfile, cpu=runinfo[3],
      ref, anc,
      minMapQ=30,minQ=20,
      glikehood=1,
      doMajorMinor=1, doMaf=1,
      doGlf=3, SNP_pval=2e-6
    )
    set_angsd_indF(outfile, nInd)
    message(sprintf("### get gentoype likelihoods and then indF!"))
  }

  ### setup shell id
  if(type=="theta"){
    set_angsd_theta(
      shfile,
      bamlist, outfile, ref, cpu=runinfo[3], nInd, minInd,
      minMapQ=30, minQ=20,
      glikehood=1, anc,
      doMajorMinor=1, doMaf=1, doSaf=1, uniqueOnly=0, baq=1
    )
  }

  cmd <- paste0("sh ", shfile)
  set_farm_job(slurmsh = "slurm-script/run_angsd.sh",
               shcode = cmd, wd = NULL, jobid = "angsd", email=email,
               runinfo=runinfo)
}


#######
set_angsd_theta <- function(
  shfile,
  bamlist, outfile, ref, cpu, nInd, minInd,
  minMapQ, minQ,
  glikehood, anc,
  doMajorMinor=1, doMaf=1, doSaf=2,
  uniqueOnly=0, baq=1
){

  cat(paste("### ANGSD input file written at", Sys.time(), sep=" "),

      # Now redo angsd sample allele frequency calculation by conditioning
      # on the sites that occur in both populations.
      # command1="-bam "$pop1List" -out "$outputdir"/"$pop1"_WholeGenome
      # -doMajorMinor 1 -doMaf 1 -indF "$pop1F" -doSaf 2 -uniqueOnly 0
      # -anc "$anc" -minMapQ $minMapQ -minQ 20 -nInd $nIndPop1 -minInd $minIndPop1
      # -baq 1 -ref "$ref" -GL $glikehood -P $cpu -rf $regionfile"


      ### input bam file arguments:
      paste("angsd -bam", bamlist,
      "-out", outfile,
      "-uniqueOnly", uniqueOnly, "-minMapQ", minMapQ, "-minQ", minQ,  "-baq", baq, "-ref", ref,

      ### get MAF for fixed major and minor alleles
      "-GL", glikehood, "-P", cpu,
      "-doMajorMinor", doMajorMinor,  "-doMaf", doMaf, "-doSaf", doSaf,
      "-anc", anc,

      ### filters
      "-nInd", nInd, "-minInd", minInd),
      file=shfile, sep="\n", append=FALSE)

  cat(### input bam file arguments:
      "",
      paste0("realSFS ", outfile, ".saf.idx", " -P ", cpu, " > ", outfile, ".sfs"),
      "",
      paste("angsd -bam", bamlist,
            "-out", outfile,
            "-pest", paste0(outfile, ".sfs"),
            "-uniqueOnly", uniqueOnly, "-minMapQ", minMapQ, "-minQ", minQ,  "-baq", baq, "-ref", ref,

            ### get MAF for fixed major and minor alleles
            "-GL", glikehood, "-P", cpu,
            "-doMajorMinor", doMajorMinor,  "-doMaf", doMaf, "-doSaf", doSaf,
            "-anc", anc,

            ### filters
            "-nInd", nInd, "-minInd", minInd),
      "",
      paste("thetaStat make_bed", paste0(outfile, ".thetas.gz")),
      "",
      paste("#Estimate for every Chromosome/scaffold"),
      paste("thetaStat do_stat", paste0(outfile, ".thetas.gz"), "-nChr", 2*nInd),
      "",
      paste("#Do a sliding window analysis based on the output from the make_bed command."),
      paste("thetaStat do_stat", paste0(outfile, ".thetas.gz"), "-nChr", 2*nInd, "-win 50000 -step 10000  -outnames ",
            paste0(outfile, ".thetasWindow.gz")),

      file=shfile, sep="\n", append=T)
}

#######
set_angsd_gl <- function(
  shfile, r,
  bamlist, outfile, cpu,
  ref, anc,
  minMapQ,minQ,
  glikehood,
  doMajorMinor, doMaf,
  doGlf, SNP_pval
){
  if(is.null(r)){
    cmd <- paste("angsd -bam", bamlist,
                 "-out", outfile,
                 "-nThreads", cpu,

                 "-doGlf", doGlf,
                 "-GL", glikehood,
                 "-ref", ref,
                 "-anc", anc,

                 ### get MAF for fixed major and minor alleles
                 "-doMajorMinor", doMajorMinor,
                 "-doMaf", doMaf,
                 "-minMapQ", minMapQ,
                 "-minQ", minQ,

                 ### test for polymorphic sites with the likelihood ratio test
                 "-SNP_pval", SNP_pval)
  }else{
    cmd <- paste("angsd -bam", bamlist,
                 "-out", outfile,
                 "-nThreads", cpu,
                 "-r", r,

                 "-doGlf", doGlf,
                 "-GL", glikehood,
                 "-ref", ref,
                 "-anc", anc,

                 ### get MAF for fixed major and minor alleles
                 "-doMajorMinor", doMajorMinor,
                 "-doMaf", doMaf,
                 "-minMapQ", minMapQ,
                 "-minQ", minQ,

                 ### test for polymorphic sites with the likelihood ratio test
                 "-SNP_pval", SNP_pval)
  }

  cat(cmd, file=shfile, sep="\n", append=TRUE)

}


#######
set_angsd_indF <- function(
  outfile, nInd
){
  glf.gz <- paste0(outfile, ".glf.gz")
  approx_indF <- paste0(outfile, ".approx_indF")
  indF <- paste0(outfile, ".indF")

  #N_SITES=$((`zcat testF.HWE.mafs.gz | wc -l`-1))
  cmd1 <- paste("N_SITES=$((`zcat", glf.gz, "| wc -l`-1))")

  #zcat testF.HWE.glf.gz | ../ngsF --n_ind 20 --n_sites $N_SITES --glf -
  #--min_epsilon 0.001 --out testF.approx_indF --approx_EM --seed 12345 --init_values r
  cmd2 <- paste("zcat", glf.gz, "| ngsF", "--n_ind", nInd,
                "--n_sites $N_SITES --glf -",
                "--out", approx_indF,
                "--min_epsilon 0.001 --approx_EM --seed 12345 --init_values r")

  #zcat testF.HWE.glf.gz | ../ngsF --n_ind 20 --n_sites $N_SITES --glf -
  #--min_epsilon 0.001 --out testF.indF --init_values testF.approx_indF.pars
  cmd3 <- paste("zcat", glf.gz, "| ngsF", "--n_ind", nInd,
                "--n_sites $N_SITES --glf -",
                "--out", indF,
                "--min_epsilon 0.001 --init_values", paste0(approx_indF, ".pars"))

  cat(cmd1, cmd2, cmd3,
      file=shfile, sep="\n", append=TRUE)

}
