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
#' @param gvcf A vector of g.vcf files.
#' @param outvcf File name of the output vcf files, default="mysamples.vcf".
#' @param ref.fa The full path of genome with bwa indexed reference fasta file.
#' @param gatkpwd The absolute path of GenomeAnalysisTK.jar.
#' @param includeNonVariantSites Include loci found to be non-variant after genotyping (for GenotypeGVCFs).
#' @param hardfilter Whether to filter variants. see detail about how to apply hard filters to a call set.
#'        \url{https://www.broadinstitute.org/gatk/guide/article?id=2806}
#' @param snpflt Parameters to apply the filter to the SNP call set.
#' @param indelflt Parameters to apply the filter to the Indel call set.
#' @param email Your email address that farm will email to once the job was done/failed.
#' @param runinfo Parameters specify the array job partition information.
#' A vector of c(FALSE, "bigmemh", "1"): 1) run or not, default=FALSE
#' 2) -p partition name, default=bigmemh and 3) --cpus, default=1.
#' It will pass to \code{set_array_job}.
#' @return return a batch of shell scripts.
#'
#' @examples
#' gvcf <- c("1.vcf", "2.vcf")
#' outvcf <- "out.vcf"
#' run_GATK_JointGenotype(
#' gvcf,
#' outvcf="mysamples.vcf",
#' ref.fa="~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
#' gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
#' includeNonVariantSites=FALSE,
#' hardfilter=TRUE,
#' snpflt="\"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\"",
#' indelflt="\"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\"",
#' email=NULL,
#' runinfo = c(FALSE, "bigmemh", 1) )
#'
#' @export
run_GATK_JointGenotype <- function(
  gvcf,
  outvcf,
  ref.fa="~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
  gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
  includeNonVariantSites=FALSE,
  hardfilter=TRUE,
  snpflt="\"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\"",
  indelflt="\"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\"",
  email=NULL,
  runinfo = c(FALSE, "bigmemh", 1) ){

  ### prepare parameters for running:
  runinfo <- get_runinfo(runinfo)
  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  shid <- paste0("slurm-script/run_jointgenotype_", 1, ".sh")

  ### header of the shell code
  cat("### GATK pipeline created by farmeR",
      paste("###", format(Sys.time(), "%a %b %d %X %Y")),
      paste(""),
      file=shid, sep="\n", append=FALSE)

  set_jointgenotype(gvcf, outvcf, gatkpwd, ref.fa, includeNonVariantSites, runinfo, shid)

  if(hardfilter) set_hardfilter(outvcf, gatkpwd, snpflt, indelflt, ref.fa, runinfo, shid)

  shcode <- paste("sh slurm-script/run_jointgenotype_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  set_array_job(shid="slurm-script/run_jointgenotype_array.sh",
                shcode=shcode, arrayjobs="1",
                wd=NULL, jobid="jgeno", email=email, runinfo=runinfo)
  #  sbatch -p bigmemh --mem 32784 --ntasks=4  slurm-script/run_gatk_array.sh
}


set_jointgenotype <- function(gvcf, outvcf, gatkpwd, ref.fa, includeNonVariantSites, runinfo, shid){

  cat("### Performs joint genotyping on all samples together",
      paste0("java -Xmx", floor(as.numeric(runinfo[4])/1024), "g ", "-jar ", gatkpwd, " –T GenotypeGVCFs \\"),
      paste0("–R ", ref.fa, " \\"),
      file=shid, sep="\n", append=TRUE)

  for(i in 1:length(gvcf)){
    cat(paste0("–V ", gvcf[i],	" \\"),
        file=shid, sep="\n", append=TRUE)
  }

  if(includeNonVariantSites){
    cat("--includeNonVariantSites \\",
        file=shid, sep="\n", append=TRUE)
  }
  cat(paste0("–o ", outvcf),
      "",
      file=shid, sep="\n", append=TRUE)

  message("###>>> set up Variants calling using GATK GenotypeGVCFs!")
}


set_VQSR <- function(gvcf, outvcf, gatkpwd, ref.fa, includeNonVariantSites, runinfo, shid){
  ### Recalibrate variant quality scores = run VQSR
  ## link https://www.broadinstitute.org/gatk/guide/article?id=2805
  message("###>>> set up Recalibrate variant quality scores = run GATK VQSR!")
}


#' @rdname run_GATK_JointGenotype
#' @export
set_hardfilter <- function(outvcf, gatkpwd, snpflt, indelflt, ref.fa, runinfo, shid){

  raw_snps.vcf <- gsub("vcf$", "raw_snps.vcf", outvcf)
  filtered_snps.vcf <- gsub("vcf$", "filtered_snps.vcf", outvcf)
  raw_indels.vcf <- gsub("vcf$", "raw_indels.vcf", outvcf)
  filtered_indels.vcf <- gsub("vcf$", "filtered_indels.vcf", outvcf)

  cat("### Apply hard filters to a call set",
      "### 1. Extract the SNPs from the call set",
      paste0("java -Xmx", floor(as.numeric(runinfo[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T SelectVariants \\"),
      paste0("-selectType SNP \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-V ", outvcf, " \\"),
      paste0("-o ", raw_snps.vcf),
      "",
      file=shid, sep="\n", append=TRUE)

  cat("### 2. Apply the filter to the SNP call set",
      paste0("java -Xmx", floor(as.numeric(runinfo[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T VariantFiltration \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-V ", raw_snps.vcf, " \\"),
      paste0("--filterExpression ", snpflt, " \\"),
      paste0("--filterName \"my_snp_filter\" \\"),
      paste0("-o ", filtered_snps.vcf),
      "",
      file=shid, sep="\n", append=TRUE)

  cat("### 3. Extract the Indels from the call set",
      paste0("java -Xmx", floor(as.numeric(runinfo[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T SelectVariants \\"),
      paste0("-selectType INDEL \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-V ", outvcf, " \\"),
      paste0("-o ", raw_indels.vcf),
      "",
      file=shid, sep="\n", append=TRUE)

  cat("### 4. Apply the filter to the Indel call set",
      paste0("java -Xmx", floor(as.numeric(runinfo[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T VariantFiltration \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-V ", raw_indels.vcf, " \\"),
      paste0("--filterExpression ", indelflt, " \\"),
      paste0("--filterName \"my_indel_filter\" \\"),
      paste0("-o ", filtered_indels.vcf),
      "",
      file=shid, sep="\n", append=TRUE)

  message("###>>> Apply hard filters to a call set!")
}



