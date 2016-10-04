[![Build Status](https://travis-ci.org/dmlc/xgboost.svg?branch=master)](https://travis-ci.org/dmlc/xgboost)
[![Documentation Status](https://readthedocs.org/projects/xgboost/badge/?version=latest)](https://xgboost.readthedocs.org)
[![GitHub license](http://dmlc.github.io/img/apache2.svg)](./LICENSE)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/xgboost)](http://cran.r-project.org/web/packages/xgboost)
[![PyPI version](https://badge.fury.io/py/xgboost.svg)](https://pypi.python.org/pypi/xgboost/)
[![Gitter chat for developers at https://gitter.im/dmlc/xgboost](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/dmlc/xgboost?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[Documentation](https://xgboost.readthedocs.org) |
[Resources](demo/README.md) |
[Installation](https://xgboost.readthedocs.org/en/latest/build.html) |
[Release Notes](NEWS.md) |
[RoadMap](https://github.com/dmlc/xgboost/issues/873)

This is an R packages to generate genomic and bioinformatic pipelines and submit jobs on HPC system running slurm.

## Install

Install [devtools](https://github.com/hadley/devtools) first, and then use `devtools` to install `imputeR` from github.

```R
#install.packages(devtools)
devtools::install_github("yangjl/farmeR")
library(farmeR)
```

List all the functions in the package and find help.

```R
ls(getNamespace("farmeR"), all.names=TRUE)
?run_GATK
```

## RUN GATK for varinat calling in two steps:

### Dependencies:
1. GenomeAnalysisTK-3.5
2. picard-tools-2.1.1
3. bwa 0.7.5a

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



