
<!-- README.md is generated from README.Rmd. Please edit that file -->

# copynumber R package with support for extra assembly

This package is a fork of Bioconductor R package
[‘copynumber’](https://bioconductor.org/packages/release/bioc/html/copynumber.html)
with minor modification for supporting extra genome assemblies
(e.g. hg38, mm10). The idea from
<https://github.com/aroneklund/copynumber> is adopted to generate ‘hg38’
object, the process is recorded in `data-raw/`.

This package may be useful for running sequenza with `assembly = "hg38"`
and other software calling copy number segments with ‘copynumber’
package.

Vignette is removed in this package, please read official documentation
at
<https://bioconductor.org/packages/release/bioc/html/copynumber.html>.

**Hope someday the official package will adopt all necessary genome
assemblies.**

## Update

  - 2020-08-20
      - fork code from <https://github.com/igordot/copynumber/> to
        support mm10.
      - add code suggested by @nurmians to bypass bug 1 row `arm.data`
        in `multipcf()`.

## NOTE

The source code comes from copynumber v1.26.0, any package updates
please inform me by issue or email.

Contribution is welcome.

## Installation

Install this modified package from GitHub:

``` r
# You can install the orignal package from bioconductor fistly
# to install all dependencies
# BiocManager::install("copynumber")

# Then install this hg38-supported package
remotes::install_github("ShixiangWang/copynumber")
```

## Test copynumber

``` r
library(copynumber)

# Test pcf ----------------------------------------------------------------

#Load the lymphoma data set:
data(lymphoma)

#Take out a smaller subset of 3 samples (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:3)

#First winsorize data to handle outliers:
wins.lymph <- winsorize(sub.lymphoma)

#Run pcf (using small gamma because of low-density data):
pcf.segments <- pcf(data=wins.lymph,gamma=12,Y=sub.lymphoma, assembly = "hg38")


# Test aspcf --------------------------------------------------------------

#Load LogR and BAF data:
data(logR)
data(BAF)

#First winsorize logR to handle outliers:
wins.logR <- winsorize(logR)

#Run aspcf:
aspcf.segments <- aspcf(wins.logR,BAF, assembly = "hg38")
```

## Test sequenza

``` r
library(sequenza)
data.file = system.file("extdata", "example.seqz.txt.gz", package="sequenza", mustWork = TRUE)
test = sequenza.extract(data.file, assembly="hg38")
```
