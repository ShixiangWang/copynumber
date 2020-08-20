
####################################################################
## Author: Gro Nilsen, Knut Liest?l and Ole Christian Lingj?rde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liest?l et al. (2012), BMC Genomics
####################################################################

## Required by:
### none


## Requires:
### getMad
### getArms
### numericArms
### numericChrom
### pullOutContent

## Main function for multipcf-analysis to be called by the user

multipcf <- function(data, pos.unit = "bp", arms = NULL, Y = NULL, gamma = 40, normalize = TRUE, w = 1, fast = TRUE, assembly = "hg19", digits = 4, return.est = FALSE, save.res = FALSE, file.names = NULL, verbose = TRUE) {

  # Check pos.unit input:
  if (!pos.unit %in% c("bp", "kbp", "mbp")) {
    stop("'pos.unit' must be one of bp, kbp and mbp", call. = FALSE)
  }

  # Check assembly input:
  if (!assembly %in% c("hg19", "hg18", "hg17", "hg16", "mm7", "mm8", "mm9", "hg38", "mm10")) {
    stop("assembly must be one of hg19, hg18, hg17 or hg16", call. = FALSE)
  }

  # Is data a file:
  isfile.data <- class(data) == "character"

  # Check data input:
  if (!isfile.data) {
    # Input could come from winsorize and thus be a list; check and possibly retrieve data frame wins.data
    data <- pullOutContent(data, what = "wins.data")
    stopifnot(ncol(data) >= 4) # something is missing in input data
    # Extract information from data:
    chrom <- data[, 1]
    position <- data[, 2]
    nSample <- ncol(data) - 2
    sampleid <- colnames(data)[-c(1:2)]
  } else {
    # data is a datafile which contains data
    f <- file(data, "r") # open file connection
    head <- scan(f, nlines = 1, what = "character", quiet = TRUE, sep = "\t") # Read header
    if (length(head) < 3) {
      stop("Data in file must have at least 3 columns", call. = FALSE)
    }
    sampleid <- head[-c(1:2)]
    nSample <- length(sampleid)

    # Read just the two first columns to get chrom and pos
    chrom.pos <- read.table(file = data, sep = "\t", header = TRUE, colClasses = c(rep(NA, 2), rep("NULL", nSample)), as.is = TRUE)
    chrom <- chrom.pos[, 1]
    position <- chrom.pos[, 2]
  }

  # Make sure chrom is not factor:
  if (is.factor(chrom)) {
    # If chrom is factor; convert to character
    chrom <- as.character(chrom)
  }
  # Make sure chromosomes are numeric (replace X and Y by 23 and 24)
  num.chrom <- numericChrom(chrom)
  nProbe <- length(num.chrom)

  # Make sure position is numeric:
  if (!is.numeric(position)) {
    stop("input in data column 2 (posistions) must be numeric", call. = FALSE)
  }

  # Get character arms:
  if (is.null(arms)) {
    arms <- getArms(num.chrom, position, pos.unit, get(assembly))
  } else {
    stopifnot(length(arms) == nProbe)
  }
  # Translate to numeric arms:
  num.arms <- numericArms(num.chrom, arms)

  # Unique arms:
  arm.list <- unique(num.arms)
  nArm <- length(arm.list)

  # Check that w is same length as number of samples:
  if (length(w) == 1) {
    w <- rep(w, nSample)
  } else if (length(w) != nSample) {
    stop("'w' must be a single number or a vector of same length as the number of samples in 'data'", call. = FALSE)
  }

  # Check Y input:
  if (!is.null(Y)) {
    stopifnot(class(Y) %in% c("matrix", "data.frame", "character"))
    isfile.Y <- class(Y) == "character"

    if (!isfile.Y) {
      ncol.Y <- ncol(Y)
      nrow.Y <- nrow(Y)
    } else {
      f.y <- file(Y, "r")
      ncol.Y <- length(scan(f.y, nlines = 1, what = "character", quiet = TRUE, sep = "\t"))
      nrow.Y <- nrow(read.table(file = Y, sep = "\t", header = TRUE, colClasses = c(NA, rep("NULL", ncol.Y - 1)), as.is = TRUE))
    }
    if (nrow.Y != nProbe || ncol.Y != nSample + 2) {
      stop("Input Y does not represent the same number of probes and samples as found in input data", call. = FALSE)
    }
  }

  # Create folders in working directory where results are saved:
  # Initialize
  seg.names <- c("chrom", "arm", "start.pos", "end.pos", "n.probes", sampleid)
  mpcf.names <- c("chrom", "pos", sampleid)
  segments <- data.frame(matrix(nrow = 0, ncol = nSample + 5))
  colnames(segments) <- seg.names
  if (return.est) {
    mpcf.est <- matrix(nrow = 0, ncol = nSample)
  }
  if (save.res) {
    if (is.null(file.names)) {
      # Create directory where results are to be saved
      dir.res <- "multipcf_results"
      if (!dir.res %in% dir()) {
        dir.create(dir.res)
      }
      file.names <- c(paste(dir.res, "/", "estimates.txt", sep = ""), paste(dir.res, "/", "segments.txt", sep = ""))
    } else {
      # Check that file.names is the correct length
      if (length(file.names) < 2) {
        stop("'file.names' must be of length 2", call. = FALSE)
      }
    }
  }

  # estimates must be returned from routines if return.est or save.res
  yest <- any(return.est, save.res)

  # If normalize=T, we will scale by the sample residual standard error. If the number of probes in the data set < 100K, the MAD sd-estimate is calculated using all obs for the sample:
  sd <- rep(1, nSample) # to have a value to check in if-test later; not used otherwise. sd is only used if normalize=TRUE and then these values are replaced by MAD-sd
  if (normalize && nProbe < 100000) {
    for (j in 1:nSample) {
      if (!isfile.data) {
        sample.data <- data[, j + 2]
      } else {
        cc <- rep("NULL", nSample + 2)
        cc[j + 2] <- "numeric"
        # only read data for the j'th sample
        sample.data <- read.table(file = data, sep = "\t", header = TRUE, colClasses = cc)[, 1]
      }
      sd[j] <- getMad(sample.data[!is.na(sample.data)], k = 25) # Take out missing values before calculating mad
    }
  }

  # Scale gamma according to the number of samples:
  gamma <- gamma * nSample

  # run multiPCF separately on each chromosomearm:
  for (c in 1:nArm) {
    probe.c <- which(num.arms == arm.list[c])
    pos.c <- position[probe.c]
    nProbe.c <- length(probe.c)

    # get data for this arm
    if (!isfile.data) {
      arm.data <- data[probe.c, -c(1:2), drop = FALSE]
    } else {
      # Read data for this arm from file; since f is a opened connection, the reading will start on the next line which has not already been read
      # two first columns skipped
      arm.data <- read.table(f, nrows = nProbe.c, sep = "\t", colClasses = c(rep("NULL", 2), rep("numeric", nSample)))
    }

    # Check that there are no missing values:
    if (any(is.na(arm.data))) {
      stop("multiPCF cannot be run because there are missing data values, see 'imputeMissing' for imputation of missing values")
    }

    # Make sure data is numeric:
    if (any(!sapply(arm.data, is.numeric))) {
      stop("input in data columns 3 and onwards (copy numbers) must be numeric", call. = FALSE)
    }

    # If normalize=T and nProbe>=100K, we calculate the MAD sd-estimate for each sample using only obs in this arm
    if (normalize && nProbe >= 100000) {
      sd <- apply(arm.data, 2, getMad)
    }

    # Check sd; cannot normalize if sd=0 or if sd=NA:
    if (any(sd == 0) || any(is.na(sd))) {
      # not run multipcf, return mean for each sample:
      m <- apply(arm.data, 2, mean)
      dim(m) <- c(length(m), 1)
      if (yest) {
        yhat <- sapply(m, rep, nrow(arm.data))
        mpcf <- list(pcf = t(yhat), nIntervals = 1, start0 = 1, length = nrow(arm.data), mean = m)
      } else {
        mpcf <- list(nIntervals = 1, start0 = 1, length = nrow(arm.data), mean = m)
      }
    } else {
      # normalize data data (sd=1 if normalize=FALSE)
      arm.data <- sweep(arm.data, 2, sd, "/")

      # weight data (default weights is 1)
      arm.data <- sweep(arm.data, 2, w, "*")

      # Run multipcf:
      if (!fast || nrow(arm.data) < 400) {
        mpcf <- doMultiPCF(as.matrix(t(arm.data)), gamma = gamma, yest = yest) # requires samples in rows, probes in columns
        # note: returns samples in rows, estimates in columns.
      } else {
        mpcf <- selectFastMultiPcf(as.matrix(arm.data), gamma = gamma, L = 15, yest = yest) # requires samples in columns, probes in rows
      }

      # "Unweight" estimates:
      mpcf$mean <- sweep(mpcf$mean, 1, w, "/")
      if (yest) {
        mpcf$pcf <- sweep(mpcf$pcf, 1, w, "/")
      }
      # "Un-normalize" estimates:
      mpcf$mean <- sweep(mpcf$mean, 1, sd, "*")
      if (yest) {
        mpcf$pcf <- sweep(mpcf$pcf, 1, sd, "*")
      }
    }

    # Information about segments:
    nSeg <- mpcf$nIntervals
    start0 <- mpcf$start0
    n.pos <- mpcf$length
    seg.mean <- t(mpcf$mean) # get samples in columns
    posStart <- pos.c[start0]
    posEnd <- c(pos.c[start0 - 1], pos.c[nProbe.c])

    # Chromosome number and character arm id:
    chr <- unique(chrom[probe.c])
    a <- unique(arms[probe.c])
    chrid <- rep(chr, times = nSeg)
    armid <- rep(a, times = nSeg)

    # May use mean of input data or the observed data specified in Y:
    if (!is.null(Y)) {
      # get Y for this arm
      if (!isfile.Y) {
        arm.Y <- Y[probe.c, -c(1:2), drop = FALSE]
      } else {
        arm.Y <- read.table(f.y, nrows = length(probe.c), sep = "\t", colClasses = c(rep("NULL", 2), rep("numeric", nSample)))
      }
      # Make sure Y is numeric:
      if (any(!sapply(arm.Y, is.numeric))) {
        stop("input in Y columns 3 and onwards (copy numbers) must be numeric", call. = FALSE)
      }
      # Use observed data to calculate segment mean (recommended)
      seg.mean <- matrix(NA, nrow = nSeg, ncol = nSample)
      for (s in 1:nSeg) {
        seg.Y <- as.matrix(arm.Y[start0[s]:(start0[s] + n.pos[s] - 1), ])
        seg.mean[s, ] <- apply(seg.Y, 2, mean, na.rm = TRUE)
      }
    }

    # Round
    if (yest) {
      yhat <- round(mpcf$pcf, digits = digits)
    }
    seg.mean <- round(seg.mean, digits = digits)

    # Data frame:
    segments.c <- data.frame(chrid, armid, posStart, posEnd, n.pos, seg.mean, stringsAsFactors = FALSE)
    colnames(segments.c) <- seg.names

    # Should results be written to files or returned to user:
    if (save.res) {
      if (c == 1) {
        # open connection for writing to file
        w1 <- file(file.names[1], "w")
        w2 <- file(file.names[2], "w")
      }

      # Write segments to file for this arm
      write.table(segments.c, file = w2, col.names = if (c == 1) seg.names else FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
      # Write estimated multiPCF-values file for this arm:
      write.table(data.frame(chrom[probe.c], pos.c, t(yhat), stringsAsFactors = FALSE), file = w1, col.names = if (c == 1) mpcf.names else FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    }

    # Append results for this arm:
    segments <- rbind(segments, segments.c)
    if (return.est) {
      mpcf.est <- rbind(mpcf.est, t(yhat))
    }

    if (verbose) {
      cat(paste("multipcf finished for chromosome arm ", chr, a, sep = ""), "\n")
    }
  } # endfor

  # Close connections
  if (isfile.data) {
    close(f)
  }
  if (!is.null(Y)) {
    if (isfile.Y) {
      close(f.y)
    }
  }

  if (save.res) {
    cat(paste("multipcf-estimates were saved in file", file.names[1]), sep = "\n")
    close(w1)
    cat(paste("segments were saved in file", file.names[2]), sep = "\n")
    close(w2)
  }

  # return results:
  if (return.est) {
    mpcf.est <- data.frame(chrom, position, mpcf.est, stringsAsFactors = FALSE)
    colnames(mpcf.est) <- mpcf.names
    return(list(estimates = mpcf.est, segments = segments))
  } else {
    return(segments)
  }
} # endfunction




# Run exact multipcf algorithm, to be called by multipcf (main function)
doMultiPCF <- function(y, gamma, yest) {

  ## y: input matrix of copy number estimates, samples in rows
  ## gamma: penalty for discontinuities
  ## yest: logical, should estimates be returned
  N <- length(y)
  nSamples <- nrow(y)
  nProbes <- ncol(y)
  ## initialisations
  if (yest) {
    yhat <- rep(0, N)
    dim(yhat) <- c(nSamples, nProbes)
  }
  bestCost <- rep(0, nProbes)
  bestSplit <- rep(0, nProbes + 1)
  bestAver <- rep(0, N)
  dim(bestAver) <- c(nSamples, nProbes)
  Sum <- rep(0, N)
  dim(Sum) <- c(nSamples, nProbes)
  Aver <- rep(0, N)
  dim(Aver) <- c(nSamples, nProbes)
  Nevner <- rep(0, N)
  dim(Nevner) <- c(nProbes, nSamples)
  Nevner <- t(Nevner + (nProbes:1))
  eachCost <- rep(0, N)
  dim(eachCost) <- c(nSamples, nProbes)
  Cost <- rep(0, nProbes)
  ## Filling of first elements
  y1 <- y[, 1]
  Sum[, 1] <- y1
  Aver[, 1] <- y1
  bestCost[1] <- 0
  bestSplit[1] <- 0
  bestAver[, 1] <- y1
  helper <- rep(1, nSamples)
  ## Solving for gradually longer arrays. Sum accumulates
  ## values for errors for righthand plateau downward from n;
  ## this error added to gamma and the stored cost in bestCost
  ## give the total cost. Cost stores the total cost for breaks
  ## at any position below n, and which.min finds the position
  ## with lowest cost (best split). Aver is the average of the
  ## righthand plateau.
  for (n in 2:nProbes) {
    Sum[, 1:n] <- Sum[, 1:n] + y[, n]
    Aver[, 1:n] <- Sum[, 1:n] / Nevner[, (nProbes - n + 1):nProbes]
    eachCost[, 1:n] <- -(Sum[, 1:n] * Aver[, 1:n])
    Cost[1:n] <- helper %*% eachCost[, 1:n]
    Cost[2:n] <- Cost[2:n] + bestCost[1:(n - 1)] + gamma
    Pos <- which.min(Cost[1:n])
    bestCost[n] <- Cost[Pos]
    bestAver[, n] <- Aver[, Pos]
    bestSplit[n] <- Pos - 1
  }
  ## The final solution is found iteratively from the sequence
  ## of split positions stored in bestSplit and the averages
  ## for each plateau stored in bestAver
  n <- nProbes
  antInt <- 0
  while (n > 0) {
    if (yest) {
      yhat[, (bestSplit[n] + 1):n] <- bestAver[, n]
    }
    n <- bestSplit[n]
    antInt <- antInt + 1
  }
  n <- nProbes
  lengde <- rep(0, antInt)
  start0 <- rep(0, antInt)
  verdi <- rep(0, antInt * nSamples)
  dim(verdi) <- c(nSamples, antInt)
  oldSplit <- n
  antall <- antInt
  while (n > 0) {
    start0[antall] <- bestSplit[n] + 1
    lengde[antall] <- oldSplit - bestSplit[n]
    verdi[, antall] <- bestAver[, n]
    n <- bestSplit[n]
    oldSplit <- n
    antall <- antall - 1
  }
  if (yest) {
    return(list(pcf = yhat, length = lengde, start0 = start0, mean = verdi, nIntervals = antInt))
  } else {
    return(list(length = lengde, start0 = start0, mean = verdi, nIntervals = antInt))
  }
}


## Choose fast multipcf version, called by multipcf (main function)
selectFastMultiPcf <- function(x, gamma, L, yest) {
  xLength <- nrow(x)
  if (xLength < 1000) {
    result <- runFastMultiPCF(x, gamma, L, 0.15, 0.15, yest)
  } else {
    if (xLength < 3000) {
      result <- runFastMultiPCF(x, gamma, L, 0.12, 0.10, yest)
    } else {
      if (xLength < 15000) {
        result <- runFastMultiPCF(x, gamma, L, 0.12, 0.05, yest)
      } else {
        result <- runMultiPcfSubset(x, gamma, L, 0.12, 0.05, yest)
      }
    }
  }
  return(result)
}


# Fast version 1, for moderately long sequences, called by selectFastMultiPcf
runFastMultiPCF <- function(x, gamma, L, frac1, frac2, yest) {
  mark <- rep(0, nrow(x))
  mark <- sawMarkM(x, L, frac1, frac2)
  dense <- compactMulti(t(x), mark)
  compPotts <- multiPCFcompact(dense$Nr, dense$Sum, gamma)
  if (yest) {
    potts <- expandMulti(nrow(x), ncol(x), compPotts$Lengde, compPotts$mean)
    return(list(
      pcf = potts, length = compPotts$Lengde, start0 = compPotts$sta,
      mean = compPotts$mean, nIntervals = compPotts$nIntervals
    ))
  } else {
    return(list(
      length = compPotts$Lengde, start0 = compPotts$sta,
      mean = compPotts$mean, nIntervals = compPotts$nIntervals
    ))
  }
}

# Fast version 2, for very long sequences, called by selectFastMultiPcf
runMultiPcfSubset <- function(x, gamma, L, frac1, frac2, yest) {
  SUBSIZE <- 5000 # length of subsets
  antGen <- nrow(x)
  mark <- sawMarkM(x, L, frac1, frac2)
  markInit <- c(mark[1:(SUBSIZE - 1)], TRUE)
  compX <- compactMulti(t(x[1:SUBSIZE, ]), markInit)
  mark2 <- rep(FALSE, antGen)
  mark2[1:SUBSIZE] <- markMultiPotts(compX$Nr, compX$Sum, gamma, SUBSIZE)
  mark2[4 * SUBSIZE / 5] <- TRUE
  start0 <- 4 * SUBSIZE / 5 + 1
  while (start0 + SUBSIZE < antGen) {
    slutt <- start0 + SUBSIZE - 1
    markSub <- c(mark2[1:(start0 - 1)], mark[start0:slutt])
    markSub[slutt] <- TRUE
    compX <- compactMulti(t(x[1:slutt, ]), markSub)
    mark2[1:slutt] <- markMultiPotts(compX$Nr, compX$Sum, gamma, slutt)
    start0 <- start0 + 4 * SUBSIZE / 5
    mark2[start0 - 1] <- TRUE
  }
  markSub <- c(mark2[1:(start0 - 1)], mark[start0:antGen])
  compX <- compactMulti(t(x), markSub)
  compPotts <- multiPCFcompact(compX$Nr, compX$Sum, gamma)
  if (yest) {
    potts <- expandMulti(nrow(x), ncol(x), compPotts$Lengde, compPotts$mean)
    return(list(
      pcf = potts, length = compPotts$Lengde, start0 = compPotts$sta,
      mean = compPotts$mean, nIntervals = compPotts$nIntervals
    ))
  } else {
    return(list(
      length = compPotts$Lengde, start0 = compPotts$sta,
      mean = compPotts$mean, nIntervals = compPotts$nIntervals
    ))
  }
}

# function that accumulates numbers of observations and sums between potential breakpoints
compactMulti <- function(y, mark) {
  antGen <- ncol(y)
  antSample <- nrow(y)
  antMark <- sum(mark)
  ant <- rep(0, antMark)
  sum <- rep(0, antMark * antSample)
  dim(sum) <- c(antSample, antMark)
  pos <- 1
  oldPos <- 0
  count <- 1
  delSum <- rep(0, antSample)
  while (pos <= antGen) {
    delSum <- 0
    while (mark[pos] < 1) {
      delSum <- delSum + y[, pos]
      pos <- pos + 1
    }
    ant[count] <- pos - oldPos
    sum[, count] <- delSum + y[, pos]
    oldPos <- pos
    pos <- pos + 1
    count <- count + 1
  }
  list(Nr = ant, Sum = sum)
}


# main calculations for fast multipcf-versions
multiPCFcompact <- function(nr, sum, gamma) {
  ## nr,sum : numbers and sums for one analysis unit,
  ## typically one chromosomal arm. Samples assumed to be in rows.
  ## gamma: penalty for discontinuities
  N <- length(nr)
  nSamples <- nrow(sum)
  ## initialisations
  yhat <- rep(0, N * nSamples)
  dim(yhat) <- c(nSamples, N)
  bestCost <- rep(0, N)
  bestSplit <- rep(0, N + 1)
  bestAver <- rep(0, N * nSamples)
  dim(bestAver) <- c(nSamples, N)
  Sum <- rep(0, N * nSamples)
  dim(Sum) <- c(nSamples, N)
  Nevner <- rep(0, N * nSamples)
  dim(Nevner) <- c(nSamples, N)
  eachCost <- rep(0, N * nSamples)
  dim(eachCost) <- c(nSamples, N)
  Cost <- rep(0, N)
  ## Filling of first elements
  Sum[, 1] <- sum[, 1]
  Nevner[, 1] <- nr[1]
  bestSplit[1] <- 0
  bestAver[, 1] <- sum[, 1] / nr[1]
  helper <- rep(1, nSamples)
  bestCost[1] <- helper %*% (-Sum[, 1] * bestAver[, 1])
  lengde <- rep(0, N)

  ## Solving for gradually longer arrays. Sum accumulates
  ## error values for righthand plateau downward from n;
  ## this error added to gamma and the stored cost in bestCost
  ## give the total cost. Cost stores the total cost for breaks
  ## at any position below n, and which.min finds the position
  ## with lowest cost (best split). Aver is the average of the
  ## righthand plateau.
  for (n in 2:N) {
    Sum[, 1:n] <- Sum[, 1:n] + sum[, n]
    Nevner[, 1:n] <- Nevner[, 1:n] + nr[n]
    eachCost[, 1:n] <- -(Sum[, 1:n]^2) / Nevner[, 1:n]
    Cost[1:n] <- helper %*% eachCost[, 1:n]
    Cost[2:n] <- Cost[2:n] + bestCost[1:(n - 1)] + gamma
    Pos <- which.min(Cost[1:n])
    cost <- Cost[Pos]
    aver <- Sum[, Pos] / Nevner[, Pos]
    bestCost[n] <- cost
    bestAver[, n] <- aver
    bestSplit[n] <- Pos - 1
  }

  ## The final solution is found iteratively from the sequence
  ## of split positions stored in bestSplit and the averages
  ## for each plateau stored in bestAver

  n <- N
  antInt <- 0
  while (n > 0) {
    yhat[, (bestSplit[n] + 1):n] <- bestAver[, n]
    antInt <- antInt + 1
    lengde[antInt] <- sum(nr[(bestSplit[n] + 1):n])
    n <- bestSplit[n]
  }
  lengdeRev <- lengde[antInt:1]
  init <- rep(0, antInt)
  init[1] <- 1
  if (antInt >= 2) {
    for (k in 2:antInt) {
      init[k] <- init[k - 1] + lengdeRev[k - 1]
    }
  }

  n <- N
  verdi <- rep(0, antInt * nSamples)
  dim(verdi) <- c(nSamples, antInt)
  bestSplit[n + 1] <- n
  antall <- antInt
  while (n > 0) {
    verdi[, antall] <- bestAver[, n]
    n <- bestSplit[n]
    antall <- antall - 1
  }

  list(Lengde = lengdeRev, sta = init, mean = verdi, nIntervals = antInt)
}

# helper function for fast version 2
markMultiPotts <- function(nr, sum, gamma, subsize) {
  ## nr,sum: numbers and sums for one analysis unit,
  ##  typically one chromosomal arm. Samples assumed to be in rows.
  ## gamma: penalty for discontinuities
  N <- length(nr)
  nSamples <- nrow(sum)
  markSub <- rep(FALSE, N)
  bestCost <- rep(0, N)
  bestSplit <- rep(0, N + 1)
  bestAver <- rep(0, N * nSamples)
  dim(bestAver) <- c(nSamples, N)
  Sum <- rep(0, N * nSamples)
  dim(Sum) <- c(nSamples, N)
  Nevner <- rep(0, N * nSamples)
  dim(Nevner) <- c(nSamples, N)
  eachCost <- rep(0, N * nSamples)
  dim(eachCost) <- c(nSamples, N)
  Cost <- rep(0, N)
  ## Filling of first elements
  Sum[, 1] <- sum[, 1]
  Nevner[, 1] <- nr[1]
  bestSplit[1] <- 0
  bestAver[, 1] <- sum[, 1] / nr[1]
  helper <- rep(1, nSamples)
  bestCost[1] <- helper %*% (-Sum[, 1] * bestAver[, 1])
  lengde <- rep(0, N)
  for (n in 2:N) {
    Sum[, 1:n] <- Sum[, 1:n] + sum[, n]
    Nevner[, 1:n] <- Nevner[, 1:n] + nr[n]
    eachCost[, 1:n] <- -(Sum[, 1:n]^2) / Nevner[, 1:n]
    Cost[1:n] <- helper %*% eachCost[, 1:n]
    Cost[2:n] <- Cost[2:n] + bestCost[1:(n - 1)] + gamma
    Pos <- which.min(Cost[1:n])
    cost <- Cost[Pos]
    aver <- Sum[, Pos] / Nevner[, Pos]
    bestCost[n] <- cost
    bestAver[, n] <- aver
    bestSplit[n] <- Pos - 1
    markSub[Pos - 1] <- TRUE
  }
  help <- findMarksMulti(markSub, nr, subsize)
  return(help)
}


# Function which finds potential breakpoints
findMarksMulti <- function(markSub, Nr, subsize) {
  ## markSub: marks in compressed scale
  ## NR: number of observations between potenstial breakpoints
  mark <- rep(FALSE, subsize) ## marks in original scale
  if (sum(markSub) < 1) {
    return(mark)
  } else {
    N <- length(markSub)
    ant <- seq(1:N)
    help <- ant[markSub]
    lengdeHelp <- length(help)
    help0 <- c(0, help[1:(lengdeHelp - 1)])
    lengde <- help - help0
    start0 <- 1
    oldStart <- 1
    startOrig <- 1
    for (i in 1:lengdeHelp) {
      start0 <- start0 + lengde[i]
      lengdeOrig <- sum(Nr[oldStart:(start0 - 1)])
      startOrig <- startOrig + lengdeOrig
      mark[startOrig - 1] <- TRUE
      oldStart <- start0
    }
    return(mark)
  }
}


## expand compact solution
expandMulti <- function(nProbes, nSamples, lengthInt, mean) {
  ## input: nr of probes, length of intervals,
  ## value in intervals; returns the expansion

  Potts <- rep(0, nProbes * nSamples)
  dim(Potts) <- c(nSamples, nProbes)
  lengthCompArr <- length(lengthInt)
  k <- 1
  for (i in 1:lengthCompArr) {
    for (j in 1:lengthInt[i]) {
      Potts[, k] <- mean[, i]
      k <- k + 1
    }
  }
  return(Potts)
}

## sawtooth-filter for multiPCF - marks potential breakpoints. Uses two
## sawtoothfilters, one lang (length L) and one short (fixed length 6)
sawMarkM <- function(x, L, frac1, frac2) {
  nrProbes <- nrow(x)
  nrSample <- ncol(x)
  mark <- rep(0, nrProbes)
  sawValue <- rep(0, nrProbes)
  filter <- rep(0, 2 * L)
  sawValue2 <- rep(0, nrProbes)
  filter2 <- rep(0, 6)
  for (k in 1:L) {
    filter[k] <- k / L
    filter[2 * L + 1 - k] <- -k / L
  }

  for (k in 1:3) {
    filter2[k] <- k / 3
    filter2[7 - k] <- -k / 3
  }

  for (l in 1:(nrProbes - 2 * L + 1)) {
    for (m in 1:nrSample) {
      diff <- crossprod(filter, x[l:(l + 2 * L - 1), m])
      sawValue[l + L - 1] <- sawValue[l + L - 1] + abs(diff)
    }
  }
  limit <- quantile(sawValue, (1 - frac1))
  for (l in 1:(nrProbes - 2 * L)) {
    if (sawValue[l + L - 1] > limit) {
      mark[l + L - 1] <- 1
    }
  }
  for (l in (L - 1):(nrProbes - L - 2)) {
    for (m in 1:nrSample) {
      diff2 <- crossprod(filter2, x[l:(l + 5), m])
      sawValue2[l + 2] <- sawValue2[l + 2] + abs(diff2)
    }
  }
  limit2 <- quantile(sawValue2, (1 - frac2))
  for (l in (L - 1):(nrProbes - L - 2)) {
    if (sawValue2[l + 2] > limit2) {
      mark[l + 2] <- 1
    }
  }
  for (l in 1:L) {
    mark[l] <- 1
    mark[nrProbes + 1 - l] <- 1
  }

  return(mark)
}
