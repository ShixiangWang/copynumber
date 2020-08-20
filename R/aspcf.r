
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
### findNN
### getArms
### getMad
### numericArms
### numericChrom
### pullOutContent

# Main function for allele-specific PCF to be called by user

aspcf <- function(logR, BAF, pos.unit = "bp", arms = NULL, kmin = 5, gamma = 40, baf.thres = c(0.1, 0.9), skew = 3, assembly = "hg19", digits = 4, return.est = FALSE, save.res = FALSE, file.names = NULL, verbose = TRUE) {

  # Check pos.unit input:
  if (!pos.unit %in% c("bp", "kbp", "mbp")) {
    stop("pos.unit must be one of bp, kbp and mbp", call. = FALSE)
  }

  # Check assembly input:
  if (!assembly %in% c("hg38", "hg19", "hg18", "hg17", "hg16", "mm7", "mm8", "mm9", "mm10")) {
    stop("assembly must be one of hg{16, 17, 18, 19, 38}, mm{7, 8, 9, 10}", call. = FALSE)
  }

  # Check if logR and BAF are files:
  isfile.logR <- class(logR) == "character"
  isfile.BAF <- class(BAF) == "character"

  # Check and extract logR-data input:
  if (!isfile.logR) {
    # Input could come from winsorize and thus be a list; check and possibly retrieve data frame wins.data
    logR <- pullOutContent(logR, what = "wins.data")
    stopifnot(ncol(logR) >= 3) # something is missing in input data
    # Extract information from logR:
    chrom <- logR[, 1]
    position <- logR[, 2]
    nSample <- ncol(logR) - 2
    sampleid <- colnames(logR)[-c(1:2)]
  } else {
    # logR is a datafile which contains logR-data
    f.logR <- file(logR, "r") # open file connection
    head <- scan(f.logR, nlines = 1, what = "character", quiet = TRUE, sep = "\t") # Read header
    if (length(head) < 3) {
      stop("Data in logR file must have at least 3 columns", call. = FALSE)
    }
    sampleid <- head[-c(1:2)]
    nSample <- length(sampleid)

    # Read just the two first columns to get chrom and pos
    chrom.pos <- read.table(file = logR, sep = "\t", header = TRUE, colClasses = c(rep(NA, 2), rep("NULL", nSample)), as.is = TRUE) # chromosomes could be character or numeric
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
    stop("input in logR column 2 (posistions) must be numeric", call. = FALSE)
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


  # Check BAF input:
  if (!isfile.BAF) {
    # Input could come from winsorize and thus be a list; check and possibly retrieve data frame wins.data
    BAF <- pullOutContent(BAF, what = "wins.data")
    ncol.BAF <- ncol(BAF)
    nrow.BAF <- nrow(BAF)
  } else {
    f.BAF <- file(BAF, "r")
    ncol.BAF <- length(scan(f.BAF, nlines = 1, what = "character", quiet = TRUE, sep = "\t"))
    nrow.BAF <- nrow(read.table(file = BAF, sep = "\t", header = TRUE, colClasses = c(NA, rep("NULL", ncol.BAF - 1)), as.is = TRUE))
  }
  if (nrow.BAF != nProbe || ncol.BAF != nSample + 2) {
    stop("Input in BAF does not represent the same number of probes and samples as found in input in logR", call. = FALSE)
  }

  # Initialize
  yhat.names <- c("chrom", "pos", sampleid)
  seg.names <- c("sampleID", "chrom", "arm", "start.pos", "end.pos", "n.probes", "logR.mean", "BAF.mean")
  segments <- data.frame(matrix(nrow = 0, ncol = 8))
  colnames(segments) <- seg.names

  if (return.est) {
    logR.yhat <- matrix(nrow = 0, ncol = nSample)
  }
  if (save.res) {
    if (is.null(file.names)) {
      # Create directory where results are to be saved
      dir.res <- "aspcf_results"
      if (!dir.res %in% dir()) {
        dir.create(dir.res)
      }
      file.names <- c(paste(dir.res, "/", "logR_estimates.txt", sep = ""), paste(dir.res, "/", "segments.txt", sep = ""))
    } else {
      # Check that file.names is the correct length
      if (length(file.names) < 2) {
        stop("'file.names' must be of length 2", call. = FALSE)
      }
    }
  }

  # estimates must be returned from routines if return.est or save.res
  yest <- any(return.est, save.res)

  # run ASPCF separately on each chromosomearm:
  for (c in 1:nArm) {
    probe.c <- which(num.arms == arm.list[c])
    pos.c <- position[probe.c]

    # Result matrices for this arm
    segments.c <- data.frame(matrix(nrow = 0, ncol = 8))
    if (yest) {
      logR.yhat.c <- matrix(nrow = length(probe.c), ncol = 0)
    }

    # Get data for this arm
    if (!isfile.logR) {
      arm.logR <- logR[probe.c, -c(1:2), drop = FALSE]
    } else {
      # Read data for this arm from file; since f is a opened connection, the reading will start on the next line which has not already been read
      # skip two first columns
      arm.logR <- read.table(f.logR, nrows = length(probe.c), sep = "\t", colClasses = c(rep("NULL", 2), rep("numeric", nSample)))
    }
    if (!isfile.BAF) {
      arm.BAF <- BAF[probe.c, -c(1:2), drop = FALSE]
    } else {
      # Read data for this arm from file; since f is a opened connection, the reading will start on the next line which has not already been read
      arm.BAF <- read.table(f.BAF, nrows = length(probe.c), sep = "\t", colClasses = c(rep("NULL", 2), rep("numeric", nSample)))
    }

    # Checking that there are no missing values in logR:
    if (any(is.na(arm.logR))) {
      stop("Missing values are not allowed in logR input", call. = FALSE)
    }

    # Make sure data is numeric:
    if (any(!sapply(arm.logR, is.numeric))) {
      stop("input in logR columns 3 and onwards must be numeric", call. = FALSE)
    }
    if (any(!sapply(arm.BAF, is.numeric))) {
      stop("input in BAF columns 3 and onwards must be numeric", call. = FALSE)
    }

    # Run ASPCF separately for each sample:
    for (i in 1:nSample) {
      sample.logR <- arm.logR[, i]
      sample.BAF <- arm.BAF[, i]

      # Initialize:
      yhat1 <- rep(NA, length(probe.c))
      yhat2 <- rep(NA, length(probe.c))

      # Filter out BAF-values below/above threshold:
      filt.baf <- sample.BAF
      filt.baf[filt.baf < baf.thres[1]] <- NA
      filt.baf[filt.baf > baf.thres[2]] <- NA

      obs <- !is.na(filt.baf)

      if (sum(obs) != 0) { # Making sure at least one BAF-value in this arm and this sample passes the threshold test!
        # Find nearest non-missing neighbour for obs. that have been filtered:
        nn <- rep(NA, length(probe.c))
        nn[obs] <- which(obs)
        if (any(!obs)) {
          nn[!obs] <- findNN(pos = pos.c, obs = obs)
          # aggregate logR values for probes with the same NN:
          use.logR <- aggregate(sample.logR, by = list(nn), FUN = mean)$x
        } else {
          use.logR <- sample.logR
        }
        # use BAF-values inside threshold:
        use.BAF <- filt.baf[obs]

        # Run ASPCF
        res <- fastAspcf(logR = use.logR, allB = use.BAF, kmin = kmin, gamma = gamma, skewed.SD = skew)

        yhat1[obs] <- res$yhat1
        yhat2[obs] <- res$yhat2

        # Interpolate over filtered positions:
        yhat1[!obs] <- yhat1[nn[!obs]]
        yhat2[!obs] <- yhat2[nn[!obs]]
      } else {
        yhat1 <- rep(mean(sample.logR), length(probe.c))
        yhat2 <- rep(NA, length(probe.c))
        message(paste("aspcf is not run for ", sampleid[i], " in chromosome arm ", unique(chrom[probe.c]), unique(arms[probe.c]), " because all of the BAF-values are outside the threshold values. Mean is returned for logR.", sep = ""))
      }
      # Rounding:
      yhat1 <- round(yhat1, digits = digits)
      yhat2 <- round(yhat2, digits = digits)

      # Create segmentation information:
      # Note that yhat1 and yhat2 will have the same breakpoints; hence only use one of them to find these
      wd <- which(diff(yhat1) != 0)
      seg.start <- c(1, wd + 1)
      seg.stop <- c(wd, length(probe.c))
      nSeg <- length(seg.start)

      # Create table with relevant segment-information
      seg.arm <- rep(unique(arms[probe.c]), nSeg)
      seg.chrom <- rep(unique(chrom[probe.c]), nSeg)
      pos.start <- pos.c[seg.start]
      pos.stop <- pos.c[seg.stop]
      n.pos <- seg.stop - seg.start + 1

      logR.mean <- yhat1[seg.start]
      BAF.mean <- yhat2[seg.start]

      # Data frame:
      seg <- data.frame(rep(sampleid[i], nSeg), seg.chrom, seg.arm, pos.start, pos.stop, n.pos, logR.mean, BAF.mean, stringsAsFactors = FALSE)
      colnames(seg) <- seg.names
      segments.c <- rbind(segments.c, seg)

      if (yest) {
        logR.yhat.c <- cbind(logR.yhat.c, yhat1)
      }
    } # endfor


    # Should results be written to files or returned to user:
    if (save.res) {
      if (c == 1) {
        # open connection for writing to file
        w1 <- file(file.names[1], "w")
        w2 <- file(file.names[2], "w")
      }
      # Write estimated PCF-values file for this arm:
      write.table(data.frame(chrom[probe.c], pos.c, logR.yhat.c, stringsAsFactors = FALSE), file = w1, col.names = if (c == 1) yhat.names else FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

      # Write segments to file for this arm
      write.table(segments.c, file = w2, col.names = if (c == 1) seg.names else FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    }

    # Append to results for other arms:
    segments <- rbind(segments, segments.c)
    if (return.est) {
      logR.yhat <- rbind(logR.yhat, logR.yhat.c)
    }

    if (verbose) {
      cat(paste("aspcf finished for chromosome arm ", seg.chrom[1], seg.arm[1], sep = ""), "\n")
    }
  } # endfor

  if (isfile.logR) {
    # Close connection
    close(f.logR)
  }
  if (isfile.BAF) {
    # Close connection
    close(f.BAF)
  }

  if (save.res) {
    close(w1)
    close(w2)

    cat(paste("logR-estimates were saved in file", file.names[1]), sep = "\n")
    cat(paste("segments were saved in file", file.names[2]), sep = "\n")
  }

  if (return.est) {
    logR_yhat <- data.frame(chrom, position, logR.yhat, stringsAsFactors = FALSE)
    return(list(logR_estimates = logR.yhat, segments = segments))
  } else {
    return(segments)
  }
} # endfunction





# fast ASPCF  version
fastAspcf <- function(logR, allB, kmin, gamma, skewed.SD) {
  N <- length(logR)
  w <- 1000 # windowsize
  d <- 100

  startw <- -d
  stopw <- w - d

  nseg <- 0
  var2 <- 0
  breakpts <- 0
  larger <- TRUE
  repeat{
    from <- max(c(1, startw))
    to <- min(c(stopw, N))
    logRpart <- logR[from:to]
    allBpart <- allB[from:to]
    allBflip <- allBpart
    allBflip[allBpart > 0.5] <- 1 - allBpart[allBpart > 0.5]

    sd1 <- getMad(logRpart)
    sd2 <- getMad(allBflip)

    # Must check that sd1 and sd2 are defined and != 0:
    sd.valid <- c(!is.na(sd1), !is.na(sd2), sd1 != 0, sd2 != 0)
    if (all(sd.valid)) {
      # run aspcfpart:
      part.res <- aspcfpart(logRpart = logRpart, allBflip = allBflip, a = startw, b = stopw, d = d, sd1 = sd1, sd2 = sd2, N = N, kmin = kmin, gamma = gamma)
      breakptspart <- part.res$breakpts
      # the 'larger' is (occasionally) necessary in the last window of the segmentation!
      larger <- breakptspart > breakpts[length(breakpts)]
      breakpts <- c(breakpts, breakptspart[larger])
      var2 <- var2 + sd2^2
      nseg <- nseg + 1
    }

    if (stopw < N + d) {
      startw <- min(stopw - 2 * d + 1, N - 2 * d)
      stopw <- startw + w
    } else {
      break
    }
  } # end repeat
  breakpts <- unique(c(breakpts, N))
  if (nseg == 0) {
    nseg <- 1
  } # just in case the sd-test never passes.
  sd2 <- sqrt(var2 / nseg)

  # On each segment calculate mean of unflipped B allele data
  frst <- breakpts[1:length(breakpts) - 1] + 1
  last <- breakpts[2:length(breakpts)]
  nseg <- length(frst)

  yhat1 <- rep(NA, N)
  yhat2 <- rep(NA, N)

  for (i in 1:nseg) {
    yhat1[frst[i]:last[i]] <- rep(mean(logR[frst[i]:last[i]]), last[i] - frst[i] + 1)
    yi2 <- allB[frst[i]:last[i]]
    # Center data around zero (by subtracting 0.5) and estimate mean
    if (length(yi2) == 0) {
      mu <- 0
    } else {
      mu <- mean(abs(yi2 - 0.5))
    }

    # Make a (slightly arbitrary) decision concerning branches
    # This may be improved by a test of equal variances
    if (sqrt(sd2^2 + mu^2) < skewed.SD * sd2) {
      mu <- 0
    }
    yhat2[frst[i]:last[i]] <- rep(mu + 0.5, last[i] - frst[i] + 1)
  }

  return(list(yhat1 = yhat1, yhat2 = yhat2))
} # end fastAspcf


# Get breakpts for a given window
aspcfpart <- function(logRpart, allBflip, a, b, d, sd1, sd2, N, kmin, gamma) {
  from <- max(c(1, a))
  usefrom <- max(c(1, a + d))
  useto <- min(c(N, b - d))

  N <- length(logRpart)
  y1 <- logRpart
  y2 <- allBflip

  # Check that vectors are long enough to run algorithm:
  if (N < 2 * kmin) {
    breakpts <- 0
    return(list(breakpts = breakpts))
  }

  # Find initSum, initKvad, initAve for segment y[1..kmin]
  initSum1 <- sum(y1[1:kmin])
  initKvad1 <- sum(y1[1:kmin]^2)
  initAve1 <- initSum1 / kmin
  initSum2 <- sum(y2[1:kmin])
  initKvad2 <- sum(y2[1:kmin]^2)
  initAve2 <- initSum2 / kmin

  # Define vector of best costs
  bestCost <- rep(0, N)
  cost1 <- (initKvad1 - initSum1 * initAve1) / sd1^2
  cost2 <- (initKvad2 - initSum2 * initAve2) / sd2^2
  bestCost[kmin] <- cost1 + cost2

  # Define vector of best splits
  bestSplit <- rep(0, N)

  # Define vector of best averages
  bestAver1 <- rep(0, N)
  bestAver2 <- rep(0, N)
  bestAver1[kmin] <- initAve1
  bestAver2[kmin] <- initAve2


  # Initialize
  Sum1 <- rep(0, N)
  Sum2 <- rep(0, N)
  Kvad1 <- rep(0, N)
  Kvad2 <- rep(0, N)
  Aver1 <- rep(0, N)
  Aver2 <- rep(0, N)
  Cost <- rep(0, N)

  # We have to treat the region y(1..2*kmin-1) separately, as it
  # cannot be split into two full segments
  kminP1 <- kmin + 1
  for (k in (kminP1):(2 * kmin - 1)) {
    Sum1[kminP1:k] <- Sum1[kminP1:k] + y1[k]
    Aver1[kminP1:k] <- Sum1[kminP1:k] / ((k - kmin):1)
    Kvad1[kminP1:k] <- Kvad1[kminP1:k] + y1[k]^2
    Sum2[kminP1:k] <- Sum2[kminP1:k] + y2[k]
    Aver2[kminP1:k] <- Sum2[kminP1:k] / ((k - kmin):1)
    Kvad2[kminP1:k] <- Kvad2[kminP1:k] + y2[k]^2


    bestAver1[k] <- (initSum1 + Sum1[kminP1]) / k
    bestAver2[k] <- (initSum2 + Sum2[kminP1]) / k
    cost1 <- ((initKvad1 + Kvad1[kminP1]) - k * bestAver1[k]^2) / sd1^2
    cost2 <- ((initKvad2 + Kvad2[kminP1]) - k * bestAver2[k]^2) / sd2^2

    bestCost[k] <- cost1 + cost2
  }


  for (n in (2 * kmin):N) {
    nMkminP1 <- n - kmin + 1

    Sum1[kminP1:n] <- Sum1[kminP1:n] + y1[n]
    Aver1[kminP1:n] <- Sum1[kminP1:n] / ((n - kmin):1)
    Kvad1[kminP1:n] <- Kvad1[kminP1:n] + (y1[n])^2

    cost1 <- (Kvad1[kminP1:nMkminP1] - Sum1[kminP1:nMkminP1] * Aver1[kminP1:nMkminP1]) / sd1^2

    Sum2[kminP1:n] <- Sum2[kminP1:n] + y2[n]
    Aver2[kminP1:n] <- Sum2[kminP1:n] / ((n - kmin):1)
    Kvad2[kminP1:n] <- Kvad2[kminP1:n] + (y2[n])^2
    cost2 <- (Kvad2[kminP1:nMkminP1] - Sum2[kminP1:nMkminP1] * Aver2[kminP1:nMkminP1]) / sd2^2

    Cost[kminP1:nMkminP1] <- bestCost[kmin:(n - kmin)] + cost1 + cost2

    Pos <- which.min(Cost[kminP1:nMkminP1]) + kmin
    cost <- Cost[Pos] + gamma

    aver1 <- Aver1[Pos]
    aver2 <- Aver2[Pos]
    totAver1 <- (Sum1[kminP1] + initSum1) / n
    totCost1 <- ((Kvad1[kminP1] + initKvad1) - n * totAver1 * totAver1) / sd1^2
    totAver2 <- (Sum2[kminP1] + initSum2) / n
    totCost2 <- ((Kvad2[kminP1] + initKvad2) - n * totAver2 * totAver2) / sd2^2
    totCost <- totCost1 + totCost2

    if (totCost < cost) {
      Pos <- 1
      cost <- totCost
      aver1 <- totAver1
      aver2 <- totAver2
    }
    bestCost[n] <- cost
    bestAver1[n] <- aver1
    bestAver2[n] <- aver2
    bestSplit[n] <- Pos - 1
  } # endfor


  # Trace back
  n <- N
  breakpts <- n
  while (n > 0) {
    breakpts <- c(bestSplit[n], breakpts)
    n <- bestSplit[n]
  } # endwhile

  breakpts <- breakpts + from - 1
  breakpts <- breakpts[breakpts >= usefrom & breakpts <= useto]

  return(list(breakpts = breakpts))
} # end aspcfpart
