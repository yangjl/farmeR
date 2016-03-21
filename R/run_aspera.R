#' \code{Run Aspera Connect to download from SRA.}
#'
#' Downloading SRA using 'ascp' utility or Aspera Connect.
#'
#' see more detail about SRA with Aspera downloading:
#' \url{http://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.downloading_sra_data_using}
#'
#' @param sra An input data.frame for SRA ids. Must contains column: SRR.
#' @param maxspeed The max speed for downloading.
#' @param outdir The output directory.
#' @param arrayjobs A character specify the number of array you try to run, i.e. 1-100.
#' @param jobid The job name show up in your sq NAME column.
#' @param email Your email address that farm will email to once the job was done/failed.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' sra <- data.frame(SRR=c("ERR877647", "ERR877648"),
#' SRX=c( "ERX957210", "ERX957211"),
#' pid=c( "1_Base1_Bbreve-sc-2188486", "P1_ECvsrS_1-sc-2201977"))
#' run_aspera(sra, maxspeed="200m", outdir=".", arrayjobs="1-2", jobid="aspera", email=NULL)
#'
#' @export
run_aspera <- function(sra,
                       maxspeed="200m",
                       outdir=".",
                       arrayjobs="1-2",
                       jobid="aspera",
                       email=NULL){

  # create dir if not exist
  dir.create("slurm-script", showWarnings = FALSE)
  for(i in 1:nrow(sra)){
    shid <- paste0("slurm-script/run_aspera_", i, ".sh")
    # ascp root: vog.hin.mln.ibcn.ptf@ptfnona:
    #Remainder of path:
    #    /sra/sra-instant/reads/ByRun/sra/{SRR|ERR|DRR}/<first 6 characters of accession>/<accession>/<accession>.sra
    #Where
    #{SRR|ERR|DRR} should be either ‘SRR’, ‘ERR’, or ‘DRR’ and should match the prefix of the target .sra file
    #ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -QT -l 100m
    # anonftpftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR161/SRR1610960/SRR1610960.sra
    sraid <- sra$SRR[i]
    id1 <- substr(sraid, start=1, stop=3)
    id2 <- substr(sraid, start=1, stop=6)
    # genome ftp://ftp.ensemblgenomes.org/pub/plants/release-22/fasta/zea_mays/dna/README
    #http://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.downloading_sra_data_using

    cmd <- paste0("ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -QT -l ", maxspeed,
                  " anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/",
                  id1, "/", id2, "/", sraid,"/", sraid, ".sra ", outdir)
    cat(cmd, file=shid, sep="\n", append=FALSE)
  }

  shcode <- paste("sh slurm-script/run_aspera_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
  set_array_job(shid="slurm-script/run_aspera_array.sh",
                shcode=shcode, arrayjobs=arrayjobs,
                wd=NULL, jobid=jobid, email=email)
}
