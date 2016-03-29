# farmeR

This is an R packages to generate pipelines for maize Genomic and Genetic analysis.

## Install

Install [devtools](https://github.com/hadley/devtools) first, and then use `devtools` to install `imputeR` from github.

```R
devtools::install_github("hadley/devtools")
library(devtools)
install_github("yangjl/farmeR")
library(farmeR)
```

OR
```R
git clone 
library(devtools)
load_all(farmeR)
```

List all the functions in the package and find help.

```R
ls(getNamespace("farmeR"), all.names=TRUE)
?run_GATK
```
## Dependencies:

1. GenomeAnalysisTK-3.5
2. picard-tools-2.1.1
3. bwa 0.7.5a

## RUN GATK for varinat calling in two steps:

1. alignment, mark duplicates, realign Indel, recal bases to variant calling for PE fq files in parallel.
```R
inputdf <- data.frame(fq1="fq_1.fq", fq2="f1_2.fq", out="mysample",
                 group="g1", sample="s1", PL="illumina", LB="lib1", PU="unit1")

run_GATK(inputdf,
         ref.fa="$HOME/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
         gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
         picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
         markDup=TRUE,
         realignInDels=FALSE, indels.vcf="indels.vcf",
         recalBases=FALSE, dbsnp.vcf="dbsnp.vcf",
         email=NULL, runinfo = c(FALSE, "bigmemh", 1))
```

2. Joint genotype calling and VCF filtering.
```R
gvcf <- c("1.vcf", "2.vcf")
outvcf <- "out.vcf"
run_GATK_JointGenotype(gvcf, outvcf,
    ref.fa="~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
    gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
    hardfilter=TRUE,
    snpflt="\"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\"",
    indelflt="\"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\"",
    email=NULL,
    runinfo = c(FALSE, "bigmemh", 1) )
```
## Documentation

coming soon.



