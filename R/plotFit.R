
locatePt <- function(ptx, pty, datx, daty) {
  absdif <- abs(ptx - datx)
  indmin <- which(absdif == min(absdif))       # multiple
  indmin[which.min(abs(pty - daty[indmin]))]   # single
}

#' Plot the fit for standards and the samples
#'
#' Produces a plot that includes points for standards, proposed fit, removed outliers, bounds for "flat" portions of the curve, and values for samples and for the background.
#'
#' @details to be added
#'
#' @param fitpar  values of function parameters
#' @param FUNmod  model function
#' @param iout    indices of removed standard points.
#' @param smpflag character vector, flags for each sample
#' @param trimval for final results, the values at which the samples are trimmed
#' @param trimext integer vector of length two indicating if the values are trimmed at the extremum (lower and upper)
#' @inheritParams processLum
#' @inheritParams fitStd
#'
#' @export

plotFit <- function(xvar, yvar, std, fitpar = NULL, FUNmod = NULL, iout = NULL,
                    bg = NULL, smp = NULL, smpflag = NULL, trimval = NULL,
                    trimext = NULL, ptcol = "firebrick3",
                    rugcols = c("cadetblue", "purple", "red"), ...) {
  if (!hasArg(ylim)) {
    ylim <- range(std[, yvar], smp, log(bg), na.rm = TRUE)
    plot(std[, xvar], std[, yvar], col = ptcol, xaxt = "n", ylim = ylim, ...)
  } else {
    plot(std[, xvar], std[, yvar], col = ptcol, xaxt = "n", ...)
  }
  if(!is.null(smp)) {
    rug(smp, side = 2, col = rugcols[1])
    if (!is.null(smpflag)) {
      rug(smp[grep("lower|upper", smpflag)], side = 2, col = rugcols[2])
      rug(smp[grep("min|max",     smpflag)], side = 2, col = rugcols[3])
    }
  }

  axis(side = 1, at = std[, xvar], cex.axis = 0.7, tcl = -0.1,
       labels = parse(text = paste("frac(1, ", std[, "Dil"], ")", sep = "")))
  abline(h = log(bg), lty = 3)
  if (!is.null(iout)) {
    points(std[iout, xvar], std[iout, yvar], col = 2, pch = 4, cex = 2)
  }
  if (!is.null(trimval)) {
    abline(h = trimval, col = rugcols[2:3][trimext], lty = 4)
    legend("right", inset = 0.03, box.col = "white", bg = "white", cex = 0.9,
           col = rugcols[3:2], lty = 4, lwd = 1.5,
           title = "trimmed at", legend = c("extrema", "bounds"))
  }
  if (is.null(fitpar)) {
    legend("bottom", inset = 0.03, box.col = "grey", box.lwd = 0.8,
           bg = "white", cex = 0.9, col = c(ptcol, 1), lty = c(NA, 3),
           lwd = c(1.5, 1), pch = c(1, NA),
           legend = c("standards", "background"))
  } else {
    x <- seq(min(std[, xvar]), max(std[, xvar]), length = 200)
    y <- FUNmod(x, fitpar)
    lines(x, y, lty = 2, lwd = 2)
    legend("bottom", inset = 0.05, box.col = "grey", box.lwd = 0.8,
           bg = "white", cex = 0.9, col = c(ptcol, 1, 1), lty = c(NA, 3, 2),
           lwd = c(1.5, 1, 2), pch = c(1, NA, NA),
           legend = c("standards", "background", "fit"))
  }
}
