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
  indF="indF.txt", anc="ref.fa",
  email=NULL, runinfo = c(FALSE, "bigmemh", 1)
){

  runinfo <- get_runinfo(runinfo)
  #### create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)

  ### setup shell id
  if(type=="SFS"){
    set_angsd_sfs(
      shfile,
      bamlist, outfile, ref, cpu=runinfo[3], nInd, minInd,
      minMapQ=30, minQ=20,
      glikehood=1, indF, anc,
      doMajorMinor=1, doMaf=1, doSaf=2, uniqueOnly=0, baq=1
    )
  }

  if(type=="GL"){
    set_angsd_gl(
      shfile, bamlist, outfile, cpu=runinfo[3],
      minMapQ=30,
      glikehood=1,
      doMajorMinor=1, doMaf=1,
      doGlf=1, SNP_pval=2e-6
    )
  }

  cmd <- paste0("sh ", shfile)
  set_farm_job(slurmsh = "slurm-script/run_angsd.sh",
               shcode = cmd, wd = NULL, jobid = "angsd", email=email,
               runinfo=runinfo)
}


#######
set_angsd_sfs <- function(
  shfile,
  bamlist, outfile, ref, cpu, nInd, minInd,
  minMapQ, minQ,
  glikehood, indF, anc,
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
      "-indF", indF, "-anc", anc,

      ### filters
      "-nInd", nInd, "-minInd", minInd),
      file=shfile, sep="\n", append=FALSE)
}

#######
set_angsd_gl <- function(
  shfile,
  bamlist, outfile, cpu,
  minMapQ,
  glikehood,
  doMajorMinor=1, doMaf=1,
  doGlf, SNP_pval
){

  cat(paste("### ANGSD input file written at", Sys.time(), sep=" "),
      # $angsdir/angsd -bam $snpCallList -GL $glikehood -out $output/maize_snps_july -P $cpu
      # -doMaf 1 -indF $popF -doMajorMinor 1 -doGeno 5 -doPost 1 -postCutoff 0.95
      # -minMapQ $minMapQ -minQ 20 -minInd $minInd -rf $regionfile  -SNP_pval $SNP_pval

      ### input bam file arguments:
      paste("angsd -bam", bamlist,"-out", outfile, "-nThreads", cpu,
            ### BAM filters
            "-minMapQ", minMapQ,

            ### get MAF for fixed major and minor alleles
            "-GL", glikehood,
            "-doMajorMinor", doMajorMinor,"-doMaf", doMaf, "-doGlf", doGlf,

            ### test for polymorphic sites with the likelihood ratio test
            "-SNP_pval", SNP_pval
            ),
      file=shfile, sep="\n", append=FALSE)
}


