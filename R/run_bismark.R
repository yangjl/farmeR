#' \code{Run array Bismark job on farm}
#'
#' Run bowtie2/2.2.5
#' Run bismark/0.14.3
#'
#' Allow one mismatch 'n 1'
#'
#'
#' see more detail about Bismark:
#' \url{http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf}
#'
#' (I) Running bismark_genome_preparation
#' module load bismark/0.14.3
#' module load bowtie2/2.2.5
#' bismark_genome_preparation --bowtie2 /home/jolyang/dbcenter/AGP/AGPv2/
#'
#' (II) Running bismark
#'
#' #uses 0-based genomic start and 1-based end coordinates.
#' bismark_methylation_extractor -s --bedGraph --counts --buffer_size 10G --CX
#' --cytosine_report --genome_folder /home/jolyang/dbcenter/AGP/AGPv2 test_pe.bam
#' #<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
#'
#' @param inputdf An input data.frame object. Must contains fq1, fq2 and outbase, bam (optional).
#' @param genome The folder of genome prepared by bismark.
#' @param N Number of miss match.
#' @param align Whether to conduct alignment, default=TRUE.
#' @param outdir Folder for output.
#' @param email Your email address that farm will email to once the job was done/failed.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' input_df <- data.frame(fq1=c("f_1.fq", "t_1.fq"), fq2=c("f_1.fq", "t_2.fq"), out=c("t1", "t2"))
#' runa_bismark(input_df, genome="/home/jolyang/dbcenter/AGP/AGPv2",
#' cpu=4, outdir="/group/jrigrp4/BS_teo20/WGBS/BSM", arrayjobs="1-5", jobid="bs1-5",
#' email=NULL)
#' @export
run_bismark <- function(inputdf,
                        genome="/home/jolyang/dbcenter/AGP/AGPv2",
                        outdir="/group/jrigrp4/BS_teo20/WGBS/BSM",
                        N=1,
                        align=TRUE, email=NULL, runinfo = c(FALSE, "bigmemh", 1)){

  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  ### determine memory based on partition
  runinfo <- get_runinfo(runinfo)

  for(i in 1:nrow(inputdf)){
    if(sum(names(inputdf) %in% "bam") ==1){
      bamfile <- inputdf$bam[i]
    }else{
      bamfile <- paste0(outdir, "/", inputdf$outbase[i], "_pe.bam")
    }

    shid <- paste0("slurm-script/run_bismark_", i, ".sh")
    cmd1 <- paste("bismark --bowtie2 -N", N, genome, "-p", runinfo[3],
                  "-1", inputdf$fq1[i],  "-2", inputdf$fq2[i],
                  "--output_dir", outdir,  "--basename", inputdf$outbase[i])
    cmd2 <- paste("bismark_methylation_extractor -p --bedGraph --counts --buffer_size 30%",
                  "-o", outdir,
                  "--CX --cytosine_report --genome_folder", genome, bamfile)

    if(align){
      cmd <- c(cmd1, cmd2)
    }else{
      cmd <- cmd2
    }
    cat(cmd, file=shid, sep="\n", append=FALSE)
  }

  shcode <- paste("module load bismark/0.14.3",
                  "module load bowtie2/2.2.5",
                  "sh slurm-script/run_bismark_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  set_array_job(shid="slurm-script/run_bismark_array.sh",
                shcode=shcode, arrayjobs=paste("1", nrow(inputdf), sep="-"),
                wd=NULL, jobid="bismark", email=email, runinfo=runinfo)
}

