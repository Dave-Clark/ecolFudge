# parallel rarefaction curve
quickRareCurve <- function (x, step = 1, sample, xlab = "Sample Size",
  ylab = "Species", label = TRUE, col, lty, max.cores = T, nCores = 1, ...)
{
    require(parallel)
    x <- as.matrix(x)
    if (!identical(all.equal(x, round(x)), TRUE))
        stop("function accepts only integers (counts)")
    if (missing(col))
        col <- par("col")
    if (missing(lty))
        lty <- par("lty")
    tot <- rowSums(x) # calculates library sizes
    S <- specnumber(x) # calculates n species for each sample
    if (any(S <= 0)) {
        message("empty rows removed")
        x <- x[S > 0, , drop = FALSE]
        tot <- tot[S > 0]
        S <- S[S > 0]
    } # removes any empty rows
    nr <- nrow(x) # number of samples
    col <- rep(col, length.out = nr)
    lty <- rep(lty, length.out = nr)
    # parallel mclapply
    # set number of cores
    mc <- getOption("mc.cores", ifelse(max.cores, detectCores(), nCores))
    message(paste("Using ", mc, " cores"))
    out <- mclapply(seq_len(nr), mc.cores = mc, function(i) {
        n <- seq(1, tot[i], by = step)
        if (n[length(n)] != tot[i])
            n <- c(n, tot[i])
        drop(rarefy(x[i, ], n))
    })
    Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
    Smax <- sapply(out, max)
     plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab,
       type = "n", ...)
    if (!missing(sample)) {
      abline(v = sample)
      rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"),
         y = z, xout = sample, rule = 1)$y)
      abline(h = rare, lwd = 0.5)
      }
    for (ln in seq_along(out)) {
      N <- attr(out[[ln]], "Subsample")
      lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
      }
    if (label) {
      ordilabel(cbind(tot, S), labels = rownames(x), ...)
      }
    invisible(out)
}

# try running sub-samples in parallel rather than samples
# return dataframe for alternative plots
quickerRareCurve <- function (x, step = 1, sample, xlab = "Sample Size",
  ylab = "Species", label = TRUE, col, lty, max.cores = T, nCores = 1, ...)
{
    require(parallel)
    x <- as.matrix(x)
    if (!identical(all.equal(x, round(x)), TRUE))
        stop("function accepts only integers (counts)")
    if (missing(col))
        col <- par("col")
    if (missing(lty))
        lty <- par("lty")
    tot <- rowSums(x) # calculates library sizes
    S <- specnumber(x) # calculates n species for each sample
    if (any(S <= 0)) {
        message("empty rows removed")
        x <- x[S > 0, , drop = FALSE]
        tot <- tot[S > 0]
        S <- S[S > 0]
    } # removes any empty rows
    nr <- nrow(x) # number of samples
    col <- rep(col, length.out = nr)
    lty <- rep(lty, length.out = nr)
    # parallel mclapply
    # set number of cores
    mc <- getOption("mc.cores", ifelse(max.cores, detectCores(), nCores))
    # write function for rarefying one sample
    parallelRarefy <- function(x, libSize){
      rareIters <- unlist(mclapply(seq_len(libSize), mc.cores = mc, function(l){ rarefy(x, sample = l)}))
    }

    # lapply parallelised rarefy function over rows of matrix
    out <- lapply(seq_len(nr), function(i){
      parallelRarefy(x[i, ], tot[i])})

    # plotting function
    maxLib <- max(tot)
    maxSp <- max(sapply(out, max))
    plot(c(0, maxLib), c(0, maxSp), type = "n",
      xlab = xlab, ylab = ylab)
    for(ln in seq_len(length(out))) {
        N <- seq_len(length(out[[ln]]))
        lines(N, out[[ln]], lty = 1, lwd=0.8)

    }

  return(out)
}
