#' \code{Run GATK Joint Genotype on farm}
#'
#' GATK Best Practices: recommended workflows for variant discovery analysis.
#'
#' see more detail about GATK:
#' \url{https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1}
#'
#' local programs:
#' bwa Version: 0.7.5a-r405
#' picard-tools-2.1.1
#' GenomeAnalysisTK-3.5/
#'
#' @param fq An input data.frame for fastq files. Must contains fq1, fq2 and out.
#' @param kitpath The absolute or relative path of the fermi.kit directory that can invoke the pipeline.
#' @param genome The full path of genome with bwa indexed reference fasta file.
#' @param s Approximate genome size, default=3g.
#' @param t Number of threads, default=16.
#' @param l Primary read length, default=100.
#' @param arrayjobs A character specify the number of array you try to run, i.e. 1-100.
#' @param jobid The job name show up in your sq NAME column.
#' @param email Your email address that farm will email to once the job was done/failed.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' fq <- data.frame(fq1="fq_1.fq", fq2="f1_2.fq", out="mysample",
#'                  group="g1", sample="s1",PL="illumina", LB="lib1", PU="unit1")
#' run_fermikit(fq, kitpath="/home/jolyang/bin/fermikit/",
#' genome="/home/jolyang/dbcenter/AGP/AGPv2", s='3g', t=16, l=100, arrayjobs="1-2",
#' jobid="fermi", email=NULL)
#'
#' @export
run_alphaimpute <- function(
  pedfile="data/Parentage_for_imputeR.csv",
  geno,
  outputdir="largedata/teo_alpha",
  EditingParameters="95.0,2.0,98.0,AllSnpOut",
  RestartOption=1,
  email=NULL,
  runinfo = c(FALSE, "bigmemh", 4) ){

  ### prepare pedigree file
  ### It should be separated with space or comma with for missing parents coded as 0.

  p <- read.csv(pedfile)
  ped <- gsub(".*/", "", pedfile)
  ped <- gsub("\\..*", "_alpha.csv", ped)
  ped <- paste(outputdir, "/", ped)
  write.table(p, ped, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)

  message(sprintf("###>>> Load [ %s SNPs ] for [ %s plants]", nrow(geno), ncol(geno)-3))
  snpinfo <- geno[, 1:3]
  snpinfo$chr <- as.numeric(as.character(gsub("S|_.*", "", snpinfo$snpid)))
  snpinfo$pos <- as.numeric(as.character(gsub(".*_", "", snpinfo$snpid)))
  #snpinfo <- snpinfo[order(snpinfo$chr, snpinfo$pos),]

  dir.create("slurm-script", showWarnings = FALSE)
  for(chri in 1:10){
    sub <- subset(snpinfo, chr == chri)
    subgeno <- merge(sub, geno[, -2:-3])
    subgeno <- subgeno[order(subgeno$chr, subgeno$pos),]
    subinfo <- subgeno[, 1:5]
    subg <- subgeno[, -1:-5]
    tgeno <- t(subg)

    chrgeno <- paste0(outputdir, "/geno_", chri, ".txt")
    write.table(tgeno, chrgeno, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE )

    specid <- paste0("slurm-script/AlphaImputeSpec.txt")
    set_alphaimpute(ped, chrgeno, numsnp=ncol(tgeno), EditingParameters, RestartOption, runinfo, shid=specid)

    shid <- paste0("slurm-script/run_alpha_", chri, ".sh")
    cat("### AlphaImpute created by farmeR",
        paste("###", format(Sys.time(), "%a %b %d %X %Y")),
        paste("AlphaImputev1.3.2", specid),
        file=shid, sep="\n", append=FALSE)

  }


  shcode <- paste("sh slurm-script/run_alpha_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  set_array_job(shid="slurm-script/run_alpha_array.sh",
                shcode=shcode, arrayjobs="1",
                wd=NULL, jobid="alpha", email=NULL, runinfo=runinfo)
  #  sbatch -p bigmemh --mem 32784 --ntasks=4  slurm-script/run_gatk_array.sh
}

system("sbatch -p bigmemh --mem 32784 --ntasks=4 slurm-script/run_alpha_array.sh")

#' @rdname run_alphaimpute
set_alphaimpute <- function(ped, chrgeno, numsnp, EditingParameters,
                            RestartOption, runinfo, shid){

  cat(paste0("PedigreeFile                         ,", ped),
      paste0("GenotypeFile                         ,", chrgeno),
      paste0("SexChrom                             ,No"),
      paste0("NumberSnp                            ,", numsnp),
      paste0("InternalEdit                         ,Yes"),
      #EditingParameters=95.0,2.0,98.0,AllSnpOut
      paste0("EditingParameters                    ,", EditingParameters),
      #2-30 NumberOfPairsOfPhasingRounds?
      paste0("NumberPhasingRuns                    ,10"),
      paste0("CoreAndTailLengths                   ,200,300,400,500,600,250,325,410,290,700"),
      paste0("CoreLengths                          ,100,200,300,400,500,150,225,310,190,600"),
      paste0("PedigreeFreePhasing                  ,No"),
      paste0("GenotypeError                        ,1.00"),
      paste0("NumberOfProcessorsAvailable          ,", as.numeric(runinfo[3])-1),
      paste0("InternalIterations                   ,3"),
      paste0("PreprocessDataOnly                   ,No"),
      paste0("PhasingOnly                          ,No"),
      paste0("ConservativeHaplotypeLibraryUse      ,No"),
      # Those individuals with an imputation accuracy
      # above WellPhasedThreshold will be outputted in the WellPhasedIndividuals.txt file
      paste0("WellPhasedThreshold                  ,99.0"),
      paste0("UserDefinedAlphaPhaseAnimalsFile     ,None"),
      paste0("PrePhasedFile                        ,None"),
      paste0("BypassGeneProb                       ,No"),
      paste0("RestartOption                        ,", RestartOption),
      paste0("HMMOption                            ,Yes"),
      paste0("HmmParameters                        ,200,5,20,", runinfo[3], ",-123456789"),
      paste0("TrueGenotypeFile                     ,None"),
      "",
      file=shid, sep="\n", append=FALSE)
  message("###>>> Setup AlphaImputeSpec sheet!")
}




