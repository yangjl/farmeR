#' \code{Run Fermikit job on farm}
#'
#' Fermikit is a de novo assembly based variant calling pipeline for Illumina short reads
#'
#' see more detail about fermikit by Li, Heng:
#' \url{https://github.com/lh3/fermikit}
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
#' fq <- data.frame(fq1=c("f_1.fq", "t_1.fq"), fq2=c("f_1.fq", "t_2.fq"), out=c("t1", "t2"))
#' run_fermikit(fq, kitpath="/home/jolyang/bin/fermikit/",
#' genome="/home/jolyang/dbcenter/AGP/AGPv2", s='3g', t=16, l=100, arrayjobs="1-2",
#' jobid="fermi", email=NULL)
#'
#' @export
run_fermikit <- function(fq,
                         kitpath="/home/jolyang/bin/fermikit/",
                         s='3g', t=16, l=100,
                         arrayjobs="1-2",
                         jobid="fermi",
                         email=NULL){

  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:nrow(fq)){

    shid <- paste0("slurm-script/run_fermikit_", i, ".sh")
    #out <- gsub(".*/", "", out)
    #outfile <- paste0(outdir, "/", out)

    #If you have multiple FASTQ files and want to trim adapters before assembly:
    #fermi.kit/fermi2.pl unitig -s3g -t16 -l100 -p prefix \
    #"fermi.kit/seqtk mergepe r1.fq r2.fq | fermi.kit/trimadap-mt -p4" > prefix.mak
    cmd1 <- paste0(kitpath, "/fermi2.pl unitig -s", s, " -t", t, " -l", l, " -p ", fq$out[i], " \\", "\n",
                  "\"", kitpath,"/seqtk mergepe ", fq$fq1[i], " ", fq$fq2[i], " | ", " \\\n",
                  kitpath, "/trimadap-mt -p",t, "\" > ", fq$out[i])

    cmd2 <- paste0("make -f prefix.mak")
    cat(c(cmd1, cmd2), file=shid, sep="\n", append=FALSE)
  }

  shcode <- paste("sh slurm-script/run_fermikit_$SLURM_ARRAY_TASK_ID.sh", sep="\n")

  set_array_job(shid="slurm-script/run_fermikit_array.sh",
                shcode=shcode, arrayjobs=arrayjobs,
                wd=NULL, jobid=jobid, email=email)
}

