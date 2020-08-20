## code to add hg38 and mm10 genome builds


# cytoband data for various genome builds are in this file
original_sysdata <- load("R/sysdata.rda")

# hg38
cytoband_hg38 <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz"
cb <- data.table::fread(cytoband_hg38, header = FALSE, data.table = FALSE)
cb <- cb[cb$V1 %in% c(paste0("chr", 1:22), "chrX", "chrY"), ]
rownames(cb) <- seq(1:nrow(cb))
cb$V1 <- factor(cb$V1)
cb$V4 <- factor(cb$V4)
cb$V5 <- factor(cb$V5)
hg38 <- cb

# mm10
cytoband_mm10 <- "http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/cytoBand.txt.gz"
cb <- data.table::fread(cytoband_mm10, header = FALSE, data.table = FALSE)
cb <- cb[cb$V1 %in% c(paste0("chr", 1:19), "chrX", "chrY"), ]
rownames(cb) <- seq(1:nrow(cb))
cb$V1 <- factor(cb$V1)
cb$V4 <- factor(cb$V4)
cb$V5 <- factor(cb$V5)
mm10 <- cb

# add hg38 and mm10 to sysdata
save(list = c(original_sysdata, "hg38", "mm10"), file = "R/sysdata.rda")
