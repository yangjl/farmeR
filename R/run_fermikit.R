#' \code{Run Fermikit job on farm}
#'
#' Fermikit is a de novo assembly based variant calling pipeline for Illumina short reads
#'
#' see more detail about fermikit by Li, Heng:
#' \url{https://github.com/lh3/fermikit}
#'
#' @param fq An input data.frame for fastq files. Must contains fq1, fq2 and out.
#' @param kitpath The absolute or relative path of the fermi.kit directory that can invoke the pipeline.
#' @param ref.fa The full path of genome with bwa indexed reference fasta file.
#' @param bamdir Path of the bam files.
#' @param s Approximate genome size, default=2.3g.
#' @param l Primary read length, default=100.
#' @param email Your email address that farm will email to once the job was done/failed.
#' @param runinfo Parameters specify the array job partition information.
#' A vector of c(FALSE, "bigmemh", "1"): 1) run or not, default=FALSE
#' 2) -p partition name, default=bigmemh and 3) --cpus, default=1.
#' It will pass to \code{set_array_job}.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' fq <- data.frame(fq1=c("f_1.fq", "t_1.fq"), fq2=c("f_1.fq", "t_2.fq"), out=c("t1", "t2"))
#' run_fermikit(fq, kitpath="/home/jolyang/bin/fermikit/",
#'              ref.fa="path/to/fastq.fa", s='2.3g', l=100,
#'              email=NULL, runinfo = c(FALSE, "bigmemh", 1)
#'
#' run_fermikit(bamdir="", kitpath="$HOME/bin/fermikit",
#'              ref.fa="",
#'              email=NULL, runinfo = c(FALSE, "bigmemh", 1))
#'
#' @export
run_fermikit <- function(
  fq,
  kitpath="$HOME/bin/fermikit",
  ref.fa="",
  s='2.3g', l=100,
  email=NULL, runinfo = c(FALSE, "bigmemh", 1)
){

  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:nrow(fq)){

    shid <- paste0("slurm-script/run_fermikit_", i, ".sh")
    #If you have multiple FASTQ files and want to trim adapters before assembly:
    #fermi.kit/fermi2.pl unitig -s3g -t16 -l100 -p prefix \
    #"fermi.kit/seqtk mergepe r1.fq r2.fq | fermi.kit/trimadap-mt -p4" > prefix.mak
    mak <- paste0(fq$out[i], ".mak")
    mag.gz <- paste0(fq$out[i], ".mag.gz")
    gz <- paste0(fq$out[i], "*fq.gz")
    pre <- paste0(fq$out[i], ".pre.gz")
    fmd <- paste0(fq$out[i], ".flt.fmd")
    sam <- paste0(fq$out[i], ".unsrt.sam.gz")
    cmd1 <- paste0(kitpath, "/fermi2.pl unitig -s", s, " -t", runinfo[3], " -l", l, " -p ", fq$out[i], " \\", "\n",
                  "\"", kitpath,"/seqtk mergepe ", fq$fq1[i], " ", fq$fq2[i], " | ", " \\\n",
                  kitpath, "/trimadap-mt -p", runinfo[3], "\" > ", mak)

    cmd2 <- paste0("make -f ", mak)
    cmd3 <- paste0(kitpath,"/run-calling -t", runinfo[3], " ", ref.fa, " ", mag.gz, " | sh")
    cmd4 <- paste0("rm ", gz)
    cmd5 <- paste0("rm ", fmd)
    cmd6 <- paste0("rm ", pre)
    cmd7 <- paste0("rm ", sam)
    cat(c(cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7), file=shid, sep="\n", append=FALSE)
  }

  message(sprintf("###>>> mergepe, trimadap-mt and then fermi unitig !"))
  shcode <- paste("sh slurm-script/run_fermikit_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  set_array_job(shid="slurm-script/run_fermikit_array.sh",
                shcode=shcode, arrayjobs=paste("1", nrow(fq), sep="-"),
                wd=NULL, jobid="fermikit", email=email, runinfo=runinfo)
  #  sbatch -p bigmemh --mem 32784 --ntasks=4  slurm-script/run_gatk_array.sh
}

#' @rdname run_fermikit
#' @examples

#'
#' @export
run_fermikit_vcfcall <- function(
  bamdir="",
  kitpath="$HOME/bin/fermikit",
  ref.fa="",
  email=NULL, runinfo = c(FALSE, "bigmemh", 1)
){

  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  #It is also possible to call SNPs and short INDELs
  #from multiple BAMs at the same time and produce a multi-sample VCF:
  #fermi.kit/htsbox pileup -cuf ref.fa pre1.srt.bam pre2.srt.bam > out.raw.vcf
  #fermi.kit/k8 fermi.kit/hapdip.js vcfsum -f out.raw.vcf > out.flt.vcf

  cmd1 <- paste0("cd ", bamdir)
  cmd2 <- paste0(kitpath, "/htsbox pileup -cuf ", ref.fa, " *srt.bam", " > out.raw.vcf")
  cmd3 <- paste0(kitpath,"/k8 ", kitpath,"/hapdip.js vcfsum -f out.raw.vcf > out.flt.vcf")

  message(sprintf("###>>> call SNPs and short INDELs from multiple BAMs !"))
  shcode <- c(cmd1, cmd2, cmd3)
  set_array_job(shid="slurm-script/run_fermikit_array.sh",
                shcode=shcode, arrayjobs="1",
                wd=NULL, jobid="fermikit", email=email, runinfo=runinfo)
  #  sbatch -p bigmemh --mem 32784 --ntasks=4  slurm-script/run_gatk_array.sh

}
