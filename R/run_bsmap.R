#' \code{Run bsmap on farm}
#'
#' BSMAP pipeline to count C and CT bases across the genome.
#' In the pipeline, first, we will align adapter trimmed reads against reference genome.
#' Then, we will sort the bam file (bam file will be kept at the end) and remove PCR duplicates.
#' After marking duplicates, we will use bamtools to keep properly paired reads. And then use bamUtil to clip
#' overlapping reads. And finally, extract DNA methylation ratio at each cysotine site.
#'
#' see more detail about BSMAP for Methylation:
#' \url{https://sites.google.com/a/brown.edu/bioinformatics-in-biomed/bsmap-for-methylation}
#'
#' dependency:
#' bsmap-2.90
#' Note: the following parameters were used in bsmap "-v 5 -r 0 -q 20 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG".
#' bamUtil: https://github.com/statgen/bamUtil
#' module load bamtools/2.2.3
#' module load java/1.8
#'
#' @param inputdf An input data.frame for fastq files. Must contains fq1, fq2, out (and/or bam).
#' If inputdf contained bam, bwa alignment will be escaped.
#' Additional columns: group (group id), sample (sample id), PL (platform, i.e. illumina),
#' LB (library id), PU (unit, i.e. unit1). These strings (or info) will pass to BWA mem through -R.
#'
#' @param ref.fa The full path of genome with bwa indexed reference fasta file.
#' @param picardpwd The absolute path of picard.jar.
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
#' inputdf <- data.frame(fq1="$HOME/dbcenter/Ecoli/fastq/SRR2921970.sra_1.fastq.gz",
#'                       fq2="$HOME/dbcenter/Ecoli/fastq/SRR2921970.sra_2.fastq.gz",
#'                       out="$HOME/dbcenter/Ecoli/fastq/SRR2921970")
#'
#' runa_bsmap(inputdf, ref.fa="~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
#' picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
#' email=NULL, runinfo = c(FALSE, "bigmemh", 1))
#'
#' @export
runa_bsmap <- function(
  inputdf,
  ref.fa="~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
  picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
  email=NULL, runinfo = c(FALSE, "bigmemh", 1)
){

  # create dir if not exist
  runinfo <- get_runinfo(runinfo)
  dir.create("slurm-script", showWarnings = FALSE)

  for(i in 1:nrow(inputdf)){
    shid <- paste0("slurm-script/run_bsmap_", i, ".sh")
    ### header of the shell code
    cat("### BSMAP pipeline created by farmeR",
        paste("###", format(Sys.time(), "%a %b %d %X %Y")),
        paste(""),
        file=shid, sep="\n", append=FALSE)

    fq1 <- inputdf$fq1[i]
    fq2 <- inputdf$fq2[i]
    outbam <- paste0(inputdf$out[i], ".bam")
    bsalign(fq1, fq2, outbam, ref.fa, picardpwd, shid, runinfo)
  }

  shcode <- paste("module load java/1.8", "module load bamtools/2.2.3",
                  "sh slurm-script/run_bsmap_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  set_array_job(shid="slurm-script/run_bsmap_array.sh",
                shcode=shcode, arrayjobs=paste("1", nrow(inputdf), sep="-"),
                wd=NULL, jobid="bsmap", email=email, runinfo=runinfo)
}

#########
bsalign <- function(fq1, fq2, outbam, ref.fa, picardpwd, shid, runinfo){

  sorted.bam <- gsub("bam", "sorted.bam", outbam)
  markDup.bam <- gsub("bam", "markDup.bam", outbam)
  metrics <- gsub("bam", "metrics.txt", sorted.bam)
  markDup_pairs.bam <- gsub("bam", "markDup_pairs.bam", outbam)
  markDup_pairs_clipOverlap.bam <- gsub("bam", "markDup_pairs_clipOverlap.bam", outbam)
  methratio.txt <- gsub("bam", "methratio.txt", outbam)

  ### bsmap -a SMbs_16_158_S1_L001_R1_001_val_1.fq.gz -b SMbs_16_158_S1_L001_R2_001_val_2.fq.gz
  ### -d ZmB73_RefGen_v2.fa -o SMbs_DCL5_CRISPR.bam -v 5 -r 0 -p 8 -q 20 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG
  cat("### 1: align adapter trimmed datasets to B73 genome",
      paste("bsmap -a", fq1, "-b", fq2, "-d", ref.fa, "-o", outbam, "-p", runinfo[3],
            "-v 5 -r 0 -q 20 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG"),
      "",
      "### 2a, reomve PCR duplicates, must be sorted by coordinate",
      # java -jar /picard-tools-1.102/SortSam.jar INPUT=SMbs_DCL5_CRISPR.bam
      # OUTPUT=SMbs_DCL5_CRISPR_sorted.bam SORT_ORDER=coordinate
      paste0("java -Xmx", floor(as.numeric(runinfo[4])/1024), "g ",
             "-jar ", picardpwd, " SortSam \\"),
      paste0("INPUT=", outbam, " \\"),
      paste0("OUTPUT=", sorted.bam, " \\"),
      "SORT_ORDER=coordinate",
      "",
      "### 2a, reomve PCR duplicates, must be sorted by coordinate",
      # java -jar /picard-tools-1.102/SortSam.jar INPUT=SMbs_DCL5_CRISPR.bam
      # OUTPUT=SMbs_DCL5_CRISPR_sorted.bam SORT_ORDER=coordinate
      paste0("java -Xmx", floor(as.numeric(runinfo[4])/1024), "g ",
             "-jar ", picardpwd, " MarkDuplicates \\"),
      paste0("INPUT=", sorted.bam, " \\"),
      paste0("OUTPUT=", markDup.bam, " \\"),
      paste0("METRICS_FILE=", metrics, " \\"),
      "REMOVE_DUPLICATES=true",
      "",
      "### 2b, keep properly paired reads",
      paste("bamtools filter -isProperPair true -in", markDup.bam, "-out", markDup_pairs.bam ),
      "",
      "### 2c, clip overlapping reads",
      paste("bam clipOverlap --in", markDup_pairs.bam,
            "--out", markDup_pairs_clipOverlap.bam, "--stats"),
      "",
      "### 3, extract DNA methylation at each cysotine",
      paste("methratio.py -o", methratio.txt, "-d", ref.fa, "-u -z -r", markDup_pairs_clipOverlap.bam),
      paste("rm", outbam),
      paste("rm", markDup.bam),
      paste("rm", markDup_pairs.bam),
      paste("rm", markDup_pairs_clipOverlap.bam),

      file=shid, sep="\n", append=TRUE)
}





