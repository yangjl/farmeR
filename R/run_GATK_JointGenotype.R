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
run_GATK_JointGenotype <- function(
  inputdf, ref.fa="~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
  gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",

  inputbam=FALSE, indels.vcf="indels.vcf",
  dbsnp.vcf="dbsnp.vcf",
  snpflt="\"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\"",
  indelflt="\"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\"",
  email=NULL,
  run = c(TRUE, "bigmemh", 4) ){

  fq <- inputdf

  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:nrow(fq)){

    shid <- paste0("slurm-script/run_gatk_", i, ".sh")

    ### header of the shell code
    cat("### GATK pipeline created by farmeR",
        paste("###", format(Sys.time(), "%a %b %d %X %Y")),
        paste(""),
        file=shid, sep="\n", append=FALSE)

    ### alignment and sorting using picard-tools
    if(bwa){
      set_bwa(fq, run, picardpwd, shid)
    }

    #### mark duplicates
    if(markDup){
      set_markDup(fq, picardpwd, run, shid)
    }

    ### Perform local realignment around indels
    if(realignInDels){
      set_realignInDels(fq, inputbam, indels.vcf, ref.fa, gatkpwd, run, shid)
    }

    ### Recalibrate Bases
    if(recalBases){
      set_recalBases(fq, inputbam, indels.vcf, dbsnp.vcf, ref.fa, gatkpwd, run, shid)
    }

    ### Variant Discovery using HaplotypeCaller
    vcaller(fq, inputbam, ref.fa, gatkpwd, run, shid)
  }

  shcode <- paste("sh slurm-script/run_gatk_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  set_array_job(shid="slurm-script/run_gatk_array.sh",
                shcode=shcode, arrayjobs="1",
                wd=NULL, jobid="gatk", email=NULL)
  #  sbatch -p bigmemh --mem 32784 --ntasks=4  slurm-script/run_gatk_array.sh
}


set_jointgenotype <- function(gvcf, outvcf, gatkpwd, ref.fa, run, shid){

  cat("### Performs	joint	genotyping	on	all	samples	together",
      paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g ", "-jar ", gatkpwd, " –T	GenotypeGVCFs\\"),
      paste0("–R ", ref.fa, " \\"),
      file=shid, sep="\n", append=TRUE)

  for(i in 1:length(gvcf)){
    cat(paste0("–V ", gvcf[i],	" \\"),
        file=shid, sep="\n", append=TRUE)
  }
  cat(paste0("–o ", outvcf),
      file=shid, sep="\n", append=TRUE)

  message("###>>> set up Variants calling using GATK HaplotypeCaller!")
}


set_VQSR <- function(gvcf, outvcf, gatkpwd, ref.fa, run, shid){
  ### Recalibrate variant quality scores = run VQSR
  ## link https://www.broadinstitute.org/gatk/guide/article?id=2805
  message("###>>> set up Recalibrate variant quality scores = run GATK VQSR!")
}


#' @rdname run_GATK_JointGenotype
#' @export
set_hardfilter <- function(outvcf, gatkpwd, snpflt, indelflt, ref.fa, run, shid){

  raw_snps.vcf <- gsub("vcf$", "raw_snps.vcf", outvcf)
  filtered_snps.vcf <- gsub("vcf$", "filtered_snps.vcf", outvcf)
  raw_indels.vcf <- gsub("vcf$", "raw_indels.vcf", outvcf)
  filtered_indels.vcf <- gsub("vcf$", "filtered_indels.vcf", outvcf)

  cat("### Apply hard filters to a call set",
      "### 1. Extract the SNPs from the call set",
      paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T SelectVariants \\"),
      paste0("-selectType SNP \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-V ", outvcf, " \\"),
      paste0("-o ", raw_snps.vcf),
      "",
      file=shid, sep="\n", append=TRUE)

  cat("### 2. Apply the filter to the SNP call set",
      paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T VariantFiltration \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-V ", raw_snps.vcf, " \\"),
      paste0("--filterExpression ", snpflt, " \\"),
      paste0("--filterName \"my_snp_filter\" \\"),
      paste0("-o ", filtered_snps.vcf),
      "",
      file=shid, sep="\n", append=TRUE)

  cat("### 3. Extract the Indels from the call set",
      paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T SelectVariants \\"),
      paste0("-selectType INDEL \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-V ", outvcf, " \\"),
      paste0("-o ", raw_indels.vcf),
      "",
      file=shid, sep="\n", append=TRUE)

  cat("### 4. Apply the filter to the Indel call set",
      paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
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



