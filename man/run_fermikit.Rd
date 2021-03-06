% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_fermikit.R
\name{run_fermikit}
\alias{run_fermikit}
\alias{run_fermikit_vcfcall}
\title{\code{Run Fermikit job on farm}}
\usage{
run_fermikit(fq, kitpath = "$HOME/bin/fermikit", ref.fa = "", s = "2.3g",
  l = 100, email = NULL, runinfo = c(FALSE, "bigmemh", 1))

run_fermikit_vcfcall(bamdir = "", kitpath = "$HOME/bin/fermikit",
  ref.fa = "", email = NULL, runinfo = c(FALSE, "bigmemh", 1))
}
\arguments{
\item{fq}{An input data.frame for fastq files. Must contains fq1, fq2 and out.}

\item{kitpath}{The absolute or relative path of the fermi.kit directory that can invoke the pipeline.}

\item{ref.fa}{The full path of genome with bwa indexed reference fasta file.}

\item{s}{Approximate genome size, default=2.3g.}

\item{l}{Primary read length, default=100.}

\item{email}{Your email address that farm will email to once the job was done/failed.}

\item{runinfo}{Parameters specify the array job partition information.
A vector of c(FALSE, "bigmemh", "1"): 1) run or not, default=FALSE
2) -p partition name, default=bigmemh and 3) --cpus, default=1.
It will pass to \code{set_array_job}.}

\item{bamdir}{Path of the bam files.}
}
\value{
return a batch of shell scripts.
}
\description{
Fermikit is a de novo assembly based variant calling pipeline for Illumina short reads
}
\details{
see more detail about fermikit by Li, Heng:
\url{https://github.com/lh3/fermikit}
}
\examples{
fq <- data.frame(fq1=c("f_1.fq", "t_1.fq"), fq2=c("f_1.fq", "t_2.fq"), out=c("t1", "t2"))
run_fermikit(fq, kitpath="/home/jolyang/bin/fermikit/",
             ref.fa="path/to/fastq.fa", s='2.3g', l=100,
             email=NULL, runinfo = c(FALSE, "bigmemh", 1)

run_fermikit(bamdir="", kitpath="$HOME/bin/fermikit",
             ref.fa="",
             email=NULL, runinfo = c(FALSE, "bigmemh", 1))


}

