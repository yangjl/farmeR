#' \code{Run GATK job on farm}
#'
#' GATK Best Practices: recommended workflows for variant discovery analysis.
#'
#' see more detail about GATK:
#' \url{https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1}
#'
#' idxing:
#' bwa index Zea_mays.AGPv2.14.dna.toplevel.fa
#'
#' module load java/1.8
#' module load bwa/0.7.9a
#'
#' local programs:
#' bwa Version: 0.7.5a-r405
#' picard-tools-2.1.1
#' GenomeAnalysisTK-3.5/
#'
#' @param inputdf An input data.frame for fastq files. Must contains fq1, fq2, out (and/or bam).
#' If inputdf contained bam, bwa alignment will be escaped.
#' Additional columns: group (group id), sample (sample id), PL (platform, i.e. illumina),
#' LB (library id), PU (unit, i.e. unit1). These strings (or info) will pass to BWA mem through -R.
#'
#' @param ref.fa The full path of genome with bwa indexed reference fasta file.
#' @param gatkpwd The absolute path of GenomeAnalysisTK.jar.
#' @param picardpwd The absolute path of picard.jar.
#'
#' @param markDup Mark Duplicates, default=TRUE.
#' @param realignInDels Realign Indels, default=FALSE. IF TRUE, a golden indel.vcf file should be provided.
#' @param indels.vcf The full path of indels.vcf.
#' @param recalBases Recalibrate Bases, default=FALSE. IF TRUE, a golden snps.vcf file should be provided.
#' @param dbsnp.vcf The full path of dbsnp.vcf.
#' @param email Your email address that farm will email to once the jobs were done/failed.
#' @param runinfo Parameters specify the array job partition information.
#' A vector of c(FALSE, "bigmemh", "1"): 1) run or not, default=FALSE
#' 2) -p partition name, default=bigmemh and 3) --cpus, default=1.
#' It will pass to \code{set_array_job}.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' inputdf <- data.frame(fq1="fq_1.fq", fq2="f1_2.fq", out="mysample",
#'                  group="g1", sample="s1", PL="illumina", LB="lib1", PU="unit1")
#'
#' run_GATK(inputdf,
#'          ref.fa="~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
#'          gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
#'          picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
#'          markDup=TRUE,
#'          realignInDels=FALSE, indels.vcf="indels.vcf",
#'          recalBases=FALSE, dbsnp.vcf="dbsnp.vcf",
#'          email=NULL, runinfo = c(FALSE, "bigmemh", 1))
#'
#' @export
run_GATK <- function(inputdf,
                     ref.fa="~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
                     gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
                     picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
                     markDup=TRUE,
                     realignInDels=FALSE, indels.vcf="indels.vcf",
                     recalBases=FALSE, dbsnp.vcf="dbsnp.vcf",
                     email=NULL, runinfo = c(FALSE, "bigmemh", 1)){

  ##### prepare parameters:
  fq <- inputdf
  ### determine memory based on partition
  run <- get_runinfo(runinfo)

  #### create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:nrow(fq)){

    shid <- paste0("slurm-script/run_gatk_", i, ".sh")
    ### header of the shell code
    cat("### GATK pipeline created by farmeR",
        paste("###", format(Sys.time(), "%a %b %d %X %Y")),
        paste(""),
        file=shid, sep="\n", append=FALSE)

    if(sum(names(fq) %in% "bam") > 0){
      inputbam <- fq$bam[i]
    }else{
      ### alignment and sorting using picard-tools
      inputbam <- set_bwa(fq, run, picardpwd, i, ref.fa, shid)
    }

    #### mark duplicates
    if(markDup) inputbam <- set_markDup(fq, picardpwd, inputbam, i, run, shid)

    ### Perform local realignment around indels
    if(realignInDels) inputbam <- set_realignInDels(fq, inputbam, i, indels.vcf, ref.fa, gatkpwd, run, shid)

    ### Recalibrate Bases
    if(recalBases) inputbam <- set_recalBases(fq, inputbam, i, indels.vcf, dbsnp.vcf, ref.fa, gatkpwd, run, shid)

    ### Variant Discovery using HaplotypeCaller
    vcaller(fq, inputbam, i, ref.fa, gatkpwd, run, shid)
  }

  shcode <- paste("module load java/1.8", "module load bwa/0.7.9a",
                  "sh slurm-script/run_gatk_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  set_array_job(shid="slurm-script/run_gatk_array.sh",
                shcode=shcode, arrayjobs=paste("1", nrow(inputdf), sep="-"),
                wd=NULL, jobid="gatk", email=email, runinfo=runinfo)
  #  sbatch -p bigmemh --mem 32784 --ntasks=4  slurm-script/run_gatk_array.sh
}


##########
set_bwa <- function(fq, run, picardpwd, i, ref.fa, shid){
    #Generate a SAM file containing aligned reads
    #http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
    rg <- paste0("\'@RG\\tID:", fq$group[i], "\\tSM:", fq$sample[i],
                 "\\tPL:", fq$PL[i], "\\tLB:", fq$LB[i], "\\tPU:", fq$PU[i], "\'")

    aligned_sam <- paste0(fq$out[i], ".aln.sam")
    sorted_bam <- paste0(fq$out[i], ".sorted.bam")
    #Generate a SAM file containing aligned reads
    #http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates
    cat(paste("### Generate a SAM file containing aligned reads"),
        paste("bwa mem -M -R", rg, " -t", run[3], "-p", ref.fa, fq$fq1[i], fq$fq2[i], ">", aligned_sam),
        paste(""),

        ### http://broadinstitute.github.io/picard/
        paste0("java -Xmx", floor(as.numeric(run[4])/1024),
               "g -jar ", picardpwd, " SortSam \\"),
        paste0("    INPUT=", aligned_sam, " \\"),
        paste0("    OUTPUT=", sorted_bam, " \\"),
        "    SORT_ORDER=coordinate",
        #paste("rm", aligned_sam),
        paste(""),
        file=shid, sep="\n", append=TRUE)
    message("###>>> set up BWA mem and then sort to bam using picard-tools!")
    return(sorted_bam)
}

##########
set_markDup <- function(fq, picardpwd, inputbam, i, run, shid){
  ### http://broadinstitute.github.io/picard/
  sorted_bam <- inputbam
  dedup_bam <- gsub("bam$", "dedup.bam", sorted_bam)
  metrics <- paste0(fq$out[i], "_metrics.txt")

  cat(paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g ",
             "-jar ", picardpwd, " MarkDuplicates \\"),
      paste0("INPUT=", sorted_bam, " \\"),
      paste0("OUTPUT=", dedup_bam, " \\"),
      paste0("METRICS_FILE=", metrics),
      #paste("rm", sorted_bam),
      paste0(""),
      paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g ",
             "-jar $HOME/bin/picard-tools-2.1.1/picard.jar BuildBamIndex \\"),
      paste0("INPUT=", dedup_bam),
      "",
      file=shid, sep="\n", append=TRUE)
  message("###>>> set up Mark Duplicates using picard-tools!")
  return(dedup_bam)
}

##########
set_realignInDels <- function(fq, inputbam, i, indels.vcf, ref.fa, gatkpwd, run, shid){
  dir.create("$HOME/tmp", showWarnings = FALSE)

  ### input and output files
  realigned_bam <- gsub("bam", "indelrealigned.bam", inputbam)
  intervals <- paste0(fq$out[i], ".forIndelRealigner.intervals")

  cat("### Define intervals to target for local realignment",
      paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T RealignerTargetCreator \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-I ", inputbam, " \\"),
      paste0("--known ", indels.vcf, " \\"),
      paste0("-o ", intervals),
      paste(""),
      file=shid, sep="\n", append=TRUE)

  cat("### Perform local realignment around indels",
      paste0("### link: https://www.broadinstitute.org/gatk/guide/article?id=2800"),
      paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g -Djava.io.tmpdir=$HOME/tmp \\"),
      paste0("-jar ", gatkpwd, " \\"),
      paste0("-I ", inputbam, " \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-T IndelRealigner \\"),
      paste0("-targetIntervals ", intervals, " \\"),
      paste0("-o ", realigned_bam, " \\"),
      paste0("-known ", indels.vcf),
      paste0("--consensusDeterminationModel KNOWNS_ONLY \\"),
      paste0("-LOD 0.4"),
      paste(""),
      file=shid, sep="\n", append=TRUE)
  message("###>>> set up Realign InDels using GATK!")
  return(realigned_bam)

}

##########
set_recalBases <- function(fq, inputbam, i, indels.vcf, dbsnp.vcf, ref.fa, gatkpwd, run, shid){

  ### input and output files
  recal_table <- paste0(fq$out[i], ".recal_data.table")
  post_recal_table <- paste0(fq$out[i], ".post_recal_data.table")
  plotpdf <- paste0(fq$out[i], ".recalibration_plots.pdf")
  recal_bam <- gsub("bam$", "recal.bam", inputbam)

  #1. Analyze patterns of covariation in the sequence dataset
  cat("### Recalibrate base quality scores = run BQSR",
      "### link: https://www.broadinstitute.org/gatk/guide/article?id=2801",
      paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T BaseRecalibrator \\"),
      paste0("-R ", ref.fa, "\\"),
      paste0("-I ", inputbam, " \\"),
      paste0("-knownSites ", dbsnp.vcf, " \\"),
      paste0("-knownSites ", indels.vcf, " \\"),
      paste0("-o ", recal_table),
      "",
      file=shid, sep="\n", append=TRUE)

  #2. Do a second pass to analyze covariation remaining after recalibration
  cat(paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T BaseRecalibrator \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-I ", inputbam, " \\"),
      paste0("-knownSites ", dbsnp.vcf, " \\"),
      paste0("-knownSites ", indels.vcf, " \\"),
      paste0("-BQSR ", recal_table, " \\"),
      paste0("-o ", post_recal_table),
      "",
      file=shid, sep="\n", append=TRUE)

  #3. Generate before/after plots
  cat(paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T AnalyzeCovariates \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-before ", recal_table, " \\"),
      paste0("-after ", post_recal_table, " \\"),
      paste0("-plots ", plotpdf),
      "",
      file=shid, sep="\n", append=TRUE)

  #4. Apply the recalibration to your sequence data
  cat(paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T PrintReads \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-I ", inputbam, " \\"),
      paste0("-BQSR ", recal_table, " \\"),
      paste0("-o ", recal_bam),
      "",
      file=shid, sep="\n", append=TRUE)
  message("###>>> set up Recalibrate Bases using GATK!")
  return(recal_bam)
}

vcaller <- function(fq, inputbam, i, ref.fa, gatkpwd, run, shid){
  ### input and output files
  raw_variants.vcf <- paste0(fq$out[i], ".raw_variants.vcf")

  cat("### Call variants in your sequence data",
      paste0("java -Xmx", floor(as.numeric(run[4])/1024), "g ", "-jar ", gatkpwd, " \\"),
      paste0("-T HaplotypeCaller \\"),
      paste0("-R ", ref.fa, " \\"),
      paste0("-I ", inputbam, " \\"),
      paste0("--genotyping_mode DISCOVERY \\"),
      paste0("-stand_emit_conf 10 \\"),
      paste0("-stand_call_conf 30 \\"),
      paste0("-o ", raw_variants.vcf),
      file=shid, sep="\n", append=TRUE)
  message("###>>> set up Variants calling using GATK HaplotypeCaller!")

}





