## code to prepare `DATASET` dataset goes here
# Reference: https://github.com/aroneklund/copynumber

# Load the development package
devtools::load_all()

## download cytoband file for hg38
##   http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz
a <- data.table::fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz", header = FALSE)
a2 <- a[a$V1 %in% c("chrX", "chrY", paste0("chr", 1:22))]
a3 <- as.data.frame(a2)
a3$V1 <- factor(a3$V1)
a3$V4 <- factor(a3$V4)
a3$V5 <- factor(a3$V5)

hg38 <- a3

## the cytoband data for various genome builds are in this file
oldthings <- load(system.file("R", "sysdata.rda", package = "copynumber", mustWork = TRUE))

## we just add hg38 and leave the rest as it was
save(list = c(oldthings, "hg38"), file = system.file("R", "sysdata.rda", package = "copynumber", mustWork = TRUE))
