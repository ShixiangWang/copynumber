####################################################################
## Author: Gro Nilsen, Knut Liest?l and Ole Christian Lingj?rde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liest?l et al. (2012), BMC Genomics
####################################################################


# Function that detects and modifies outliers by winsorization

## Input:
### data : either a numeric matrix/data frame or the name of a tab-separated file which data can be read from. The matrix or file should have chromosome numbers in the first column, local probe positions in the second, and copy number data for one or more samples in subsequent columns, and the header of copy number column(s) should be a sample identifier
### pos.unit : the unit used to represent the probe positions. Allowed options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs). By default assumed to be "bp"
### arms : optional character vector of same length as the number of rows in \code{data} containing chromosome arms (denoted 'p' or 'q'). If not specified, arms are estimated using the function \code{getArms}
### tau : winsorization threshold, default is 2.5
### k : the window size to be applied in median filtering, default is 25
### gamma : penalty for each discontinuity in the pcf curve, default is 40
### iter : number of iterations in pcf winsorization, default is 1
### method : the winsorization method to be applied, must be one of "mad" (default) or "pcf"
### digits : the number of decimal places to be used for winsorized data values in the output. Default is 4
### return.outliers: logical, should a data frame containing outlier statuses be returned in addition to wins.data?
### save.res : logical value indicating whether winsorization results should be saved(TRUE) or returned (FALSE).
### file.names : optional character vector of length two giving the name of the files where the winsorized data and outliers, respectively, should be saved if \code{return=FALSE}
### verbose: logical value indicating whether or not to print message each time the analysis of a chromosome arm is finished

## Output:
### wins.data :a matrix with chromosomes numbers in the first column, probe positions in the second and the winsorized copy number values for the sample(s) in subsequent column(s)
### outliers : a matrix with chromosome numbers in the first column, probe positions in the second and outlier statuses for each sample in the subsequent column(s). The values +/- 1 indicate that the observation is an outlier, whereas the value 0 indicates that it is not


## Required by:
### plotGamma


## Requires:
### getArms
### numericArms
### numericChrom
### pcf (exactPcf)
### medianFilter

winsorize <- function(data, pos.unit = "bp", arms = NULL, method = "mad", tau = 2.5, k = 25, gamma = 40, iter = 1, assembly = "hg19", digits = 4, return.outliers = FALSE, save.res = FALSE, file.names = NULL, verbose = TRUE) {

  # Check pos.unit input:
  stopifnot(pos.unit %in% c("bp", "kbp", "mbp"))

  # Check method input:
  stopifnot(method %in% c("mad", "pcf"))

  # Check assembly input:
  if (!assembly %in% c("hg38", "hg19", "hg18", "hg17", "hg16", "mm7", "mm8", "mm9", "mm10")) {
    stop("assembly must be one of hg{16, 17, 18, 19, 38}, mm{7, 8, 9, 10}", call. = FALSE)
  }


  # Check data input: can either be a matrix/data frame or a filename
  stopifnot(class(data) %in% c("matrix", "data.frame", "character"))

  # Is data a file:
  isfile <- class(data) == "character"

  if (!isfile) {
    # data is matrix/data.frame:
    stopifnot(ncol(data) >= 3) # something is missing in input data
    # Extract information from data:
    chrom <- data[, 1]
    pos <- data[, 2]
    nSample <- ncol(data) - 2
    sample.names <- colnames(data)[-c(1:2)]
  } else {
    # data is a datafile which contains data
    f <- file(data, "r") # open file connection
    head <- scan(f, nlines = 1, what = "character", quiet = TRUE, sep = "\t") # Read header
    if (length(head) < 3) {
      stop("Data in file must have at least 3 columns", call. = FALSE)
    }
    sample.names <- head[-c(1:2)]
    nSample <- length(sample.names)

    # Read just the two first columns to get chrom and pos
    chrom.pos <- read.table(file = data, sep = "\t", header = TRUE, colClasses = c(rep(NA, 2), rep("NULL", nSample)), as.is = TRUE) # chromosomes could be character or numeric
    chrom <- chrom.pos[, 1]
    pos <- chrom.pos[, 2]
  }

  # Make sure chrom is not factor:
  if (is.factor(chrom)) {
    # If chrom is factor; convert to character
    chrom <- as.character(chrom)
  }

  # Make sure chromosomes are numeric (check for factor and replace X and Y by 23 and 24)
  num.chrom <- numericChrom(chrom)
  nProbe <- length(num.chrom)

  # Make sure position is numeric:
  if (!is.numeric(pos)) {
    stop("input in data column 2 (posistions) must be numeric", call. = FALSE)
  }

  # Get character arms:
  if (is.null(arms)) {
    arms <- getArms(num.chrom, pos, pos.unit, get(assembly))
  } else {
    if (length(arms) != nProbe) {
      stop("'arms' must be the same length as number of rows in data", call. = FALSE)
    }
  }
  # Translate to numeric arms:
  num.arms <- numericArms(num.chrom, arms)
  # Unique chromosome arms
  arm.list <- unique(num.arms)
  nArm <- length(arm.list)

  # Initialize
  if (!save.res) {
    wins.data <- matrix(nrow = 0, ncol = nSample)
    if (return.outliers) {
      wins.outliers <- matrix(nrow = 0, ncol = nSample)
    }
  } else {
    if (is.null(file.names)) {
      dir.res <- "Wins_res"
      if (!dir.res %in% dir()) {
        # Create folder in working directory where results are saved:
        dir.create(dir.res)
      }
      file.names <- c(paste(dir.res, "/", "wins.data.txt", sep = ""), paste(dir.res, "/", "wins.outliers.txt", sep = ""))
    } else {
      # Check that file.names is the correct length
      if (length(file.names) < 2) {
        stop("'file.names' must be of length 2", call. = FALSE)
      }
    }
  } # endif



  # Run Winsorization separately on the chromosomearms:
  for (c in 1:nArm) {
    probe.c <- which(num.arms == arm.list[c])

    # Result matrices for this arm
    wins.data.c <- matrix(nrow = length(probe.c), ncol = 0)
    if (return.outliers || save.res) {
      wins.outliers.c <- matrix(nrow = length(probe.c), ncol = 0)
    }

    if (!isfile) {
      arm.data <- data[probe.c, -c(1:2), drop = FALSE]
    } else {
      # Read data for this arm from file; since f is a opened connection, the reading will start on the next line which has not already been read
      # two first columns are skipped
      arm.data <- read.table(f, nrows = length(probe.c), sep = "\t", colClasses = c(rep("NULL", 2), rep("numeric", nSample)))
    }
    # Make sure data is numeric:
    if (any(!sapply(arm.data, is.numeric))) {
      # arm.data = sapply(arm.data,as.numeric)}
      stop("input in data columns 3 and onwards (copy numbers) must be numeric", call. = FALSE)
    }

    # Run winsorization independently on each sample:
    for (i in 1:nSample) {

      # Get sample data for this arm
      y <- arm.data[, i]

      # Only run winsorization on non-missing values
      na <- is.na(y)
      use.y <- y[!na]

      # Initialize:
      ywins <- rep(NA, length(y))
      outliers <- rep(NA, length(y))

      # Do winsorization
      wins <- switch(method,
        mad = madWins(use.y, tau = tau, k = k, digits = digits),
        pcf = pcfWins(use.y, tau = tau, k = k, gamma = gamma, iter = iter, digits = digits)
      )
      ywins[!na] <- wins$ywin
      outliers[!na] <- wins$outliers

      # Round to the desired number of digits:
      ywins <- round(ywins, digits = digits)

      # Add this sample's winsdata to arm results
      wins.data.c <- cbind(wins.data.c, ywins)
      if (return.outliers || save.res) {
        wins.outliers.c <- cbind(wins.outliers.c, outliers)
      }
    } # endfor

    if (!save.res) {
      # all winsdata are saved in one matrix
      wins.data <- rbind(wins.data, wins.data.c)
      if (return.outliers) {
        wins.outliers <- rbind(wins.outliers, wins.outliers.c)
      }
    } else {
      # Write arm-wins res to file:
      if (c == 1) {
        # open connection for writing to file
        wd <- file(file.names[1], "w")
        wo <- file(file.names[2], "w")
      }
      write.table(data.frame(chrom[probe.c], pos[probe.c], wins.data.c, stringsAsFactors = FALSE), file = wd, col.names = if (c == 1) c("chrom", "pos", sample.names) else FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

      write.table(data.frame(chrom[probe.c], pos[probe.c], wins.outliers.c, stringsAsFactors = FALSE), file = wo, col.names = if (c == 1) c("chrom", "pos", sample.names) else FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    }

    if (verbose) {
      chr <- unique(chrom[probe.c])
      a <- unique(arms[probe.c])
      cat(paste("winsorize finished for chromosome arm ", chr, a, sep = ""), "\n")
    }
  } # endfor

  if (isfile) {
    # Close connection
    close(f)
  }
  if (!save.res) {
    wins.data <- data.frame(chrom, pos, wins.data, stringsAsFactors = FALSE)
    colnames(wins.data) <- c("chrom", "pos", sample.names)
    if (return.outliers) {
      wins.outliers <- data.frame(chrom, pos, wins.outliers, stringsAsFactors = FALSE)
      colnames(wins.outliers) <- c("chrom", "pos", sample.names)
      return(list(wins.data = wins.data, wins.outliers = wins.outliers))
    } else {
      return(wins.data)
    }
  } else {
    # close connections:
    close(wd)
    close(wo)
    cat(paste("winsorized data were saved in file", file.names[1]), sep = "\n")
    cat(paste("outliers were saved in file", file.names[2]), sep = "\n")
    return(invisible(NULL))
  }
} # end winsorize



# Perform MAD winsorization:
madWins <- function(x, tau, k, digits) {
  xhat <- medianFilter(x, k)
  d <- x - xhat
  SD <- mad(d)
  z <- tau * SD
  xwin <- xhat + psi(d, z)
  outliers <- rep(0, length(x))
  outliers[round(x, digits) > round(xwin, digits)] <- 1 # Need to round here to make sure comparison works; may fail when different number of digits!!!
  outliers[round(x, digits) < round(xwin, digits)] <- -1
  return(list(ywin = xwin, sdev = SD, outliers = outliers))
}



# PCF winsorization:

pcfWins <- function(x, tau, k, gamma, iter, digits) {
  xhat <- medianFilter(x, k)
  for (j in 1:iter) {
    d <- x - xhat
    sdev <- mad(d)
    z <- tau * sdev
    xwin <- xhat + pmax(pmin(d, z), -z)
    if (length(x) < 400) {
      xhat <- exactPcf(y = xwin, gamma = gamma * sdev^2, kmin = 5, yest = TRUE)$yhat
    } else {
      xhat <- selectFastPcf(x = xwin, gamma = gamma * sdev^2, kmin = 5, yest = TRUE)$yhat
    }
  }
  sdev <- mad(xwin - xhat)
  z <- tau * sdev
  xwin <- xhat + psi(x - xhat, z)
  outliers <- rep(0, length(x))
  outliers[round(x, digits) > round(xwin, digits)] <- 1 # Need to round here to make sure comparison works; may fail when different number of digits!!!
  outliers[round(x, digits) < round(xwin, digits)] <- -1
  return(list(ywin = xwin, sdev = sdev, outliers = outliers))
}


psi <- function(x, z) {
  xwin <- x
  xwin[x < -z] <- -z
  xwin[x > z] <- z
  return(xwin)
}
