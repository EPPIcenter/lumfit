#' Process raw Luminex files
#'
#' Process data for a single antigen: fit a standard curve, establish bounds,
#' normalize samples, and save a plot showing the fit and the samples.
#'
#' details to be added
#'
#' @param antigen character string.
#' @param fname name of the file that contains raw data.
#' @param fdir directory where the file is located (alternatively, full path can
#'   be included in \code{fname}).
#' @param plotdir directory for the plots to be saved.
#' @param xvar,yvar character strings for the variables used to fit a standard
#'   curve.
#' @param model the model to be fit.
#' @param Alow lower asymptote for the sigmoid model. If \code{NULL}, the lower
#'   asymptote will be estimated, adding an extra parameter to the model. To fix
#'   the asymptote at the level of background, specify \code{"bg"}. Numeric
#'   value of \code{Alow} will force the asymptote to be fixed at the provided
#'   level.
#' @param asym if \code{TRUE}, asymmetry in the fit is allowed, adding an extra
#'   parameter to the model.
#' @param trim.flat logical value determining how the values of \code{yvar} are
#'   trimmed. If \code{TRUE}, they will be trimmed at the bounds where the curve
#'   starts to flatten out (automatically determined as maxima of the third
#'   derivative of the function). If \code{FALSE}, \code{yvar} will be trimmed
#'   at extrema, defined as the range of standards or asymptotes of the fit
#'   (whichever are less extreme).
#' @param interactive logical value. If \code{TRUE}, the user is prompted to
#'   evaluate the standards (and/or the fit) and asked to remove outliers if
#'   needed. \code{TRUE} value takes precedence over \code{rm.before} and
#'   \code{rm.after}: if both are \code{FALSE}, \code{rm.after} is reset to
#'   \code{TRUE}.
#' @param monot.prompt if \code{TRUE}, the user is prompted to evaluate the
#'   standards and possibly remove outliers if the standards are not monotonic
#'   (increasing). \code{FALSE} value is ignored if \code{interactive} is
#'   \code{TRUE}.
#' @param rm.before logical value indicating if potential outliers should be
#'   removed before the model is fitted. Ignored if \code{interactive} is
#'   \code{FALSE}.
#' @param rm.after logical value indicating if potential outliers should be
#'   removed after the model is fitted. Ignored if \code{interactive} is
#'   \code{FALSE}.
#' @param maxrm maximum number of outliers to remove.
#' @param set.bounds if \code{TRUE}, the user is given the option to manually
#'   set the bound that is not set automatically. In that case, the prompt
#'   appears even if \code{interactive} is \code{FALSE}.
#' @param overwrite.bounds logical value indicating the option to overwrite
#'   automatically set bounds. Ignored if \code{interactive} is \code{FALSE}.
#' @param ifix sorted integer vector of length 3 with indices of standards to be
#'   used for getting starting values for optimization.
#' @param dtype character string for data type in the file.
#' @param stdstr character string indicating standards in the file's "Sample"
#'   column. Not case sensitive. If \code{""} (empty string), standards will be
#'   determined by the pattern "1/" only.
#' @param bgstr character string indicating background in the file's "Sample"
#'   column. Not case sensitive.
#' @param stddil a vector of standard dilutions. If \code{NULL}, dilutions are
#'   inferred from the file's "Sample" column. Non-null value can be used to
#'   exclude some dilutions from model fitting.
#' @param smpdil dilution used for the samples.
#' @param optmethod method to be used in optimization.
#' @param maxit maximum number of iterations in optimization.
#' @param nwells number of wells. If \code{NULL}, inferred from the file.
#' @param nsep number of lines separating different data types in the file.
#' @param ncolmax maximum number of columns in the file.
#' @param width,height optional parameters for the final saved plot.
#' @param ptcol color of the standard points on the plot.
#' @param rugcols vector of three colors for the rugplot, which indicates sample
#'   values (inside the bounds, between the bounds and extrema, and beyond
#'   extrema).
#' @param ... further graphical parameters.
#'
#' @return A list of length two. The first element is a data frame that contains
#'   the results; the second is a character string with a flag containing
#'   information about removed points, failure to fit the model, manually set
#'   bounds, and/or an optional custom note provided by the user during an
#'   interactive model-fitting procedure.
#'
#' @examples
#'
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis legend lines locator mtext plot points rug
#' @importFrom methods hasArg
#' @importFrom stats optim
#' @importFrom utils read.csv
#'
#' @export

processLum <- function(antigen, fname, fdir = NULL, plotdir = NULL,
                       xvar = "logConc", yvar = "logMFI", model = "sigmoid",
                       Alow = NULL, asym = TRUE, trim.flat = TRUE,
                       interactive = TRUE, monot.prompt = FALSE,
                       rm.before = FALSE, rm.after = interactive, maxrm = 2,
                       set.bounds = interactive, overwrite.bounds = FALSE,
                       ifix = NULL, dtype = "Median", stdstr = "std|stand",
                       bgstr  = "blank|background", stddil = NULL,
                       smpdil = 1000, optmethod = "Nelder-Mead", maxit = 5e3,
                       nwells = NULL, nsep = 2, ncolmax = 105,
                       # nv3 = 10, nv4 = 100, # then include in fitStd() as well
                       width = 6, height = 6, ptcol = "firebrick3",
                       rugcols = c("cadetblue", "purple", "red"), ...) {
  MFI <- read_data(paste(fdir, fname, sep = ""), dtype, nwells, nsep, ncolmax)
  if (length(grep("1/", MFI$Sample)) == 0) {
    stop('Please provide standard concentrations in "1/dilution" format,
         e.g "1/200"')
  }
  pdate <- read.csv(paste(fdir, fname, sep = ""), header = FALSE, skip = 2,
                    nrow = 1)[1, 2]
  pdate <- as.Date(as.character(pdate), "%m/%d/%Y")
  if (nchar(stdstr) > 0 && length(grep(stdstr, MFI$Sample)) == 0) {
    warning(paste("File ", fname, ' does not contain "', stdstr, '"; using
                  "1/" to determine standards', sep = ""))
    stdstr <- ""
  }
  ctrl <- extractStd(MFI, stdstr, bgstr, stddil, smpdil, antigen, yvar)
  std <- ctrl$std

  main <- paste(antigen, ", ", gsub(".csv", "", fname), sep = "")
  finfit <- fitStd(std, xvar, yvar, model, Alow, asym, interactive,
                   monot.prompt, rm.before, rm.after, maxrm, set.bounds,
                   overwrite.bounds, ctrl$bg, ctrl$smp, optmethod, maxit,
                   info = main, ifix = ifix, rugcols = rugcols, ptcol = ptcol,
                   main = main, ...)
  if (grepl("sig", model)) {
    FUNmod <- fsig
    FUNinv <- fsigInv
  } else {
    FUNmod <- flin
    FUNinv <- flinInv
  }
  smps <- normalizeSmp(ctrl$out, antigen, fname, pdate, yvar, FUNinv,
                       finfit$par, finfit$bounds, finfit$flag, trim.flat)

  fplot <- paste(plotdir, antigen, "_", gsub(".csv", "", fname), ".pdf",
                 sep = "")
  pdf(file = fplot, width = width, height = height)
  if (!is.null(finfit$par)) {
    trimval <- smps[1, c("trim_lo", "trim_up")]
    trimext <- (trimval == smps[1, c("min", "max")]) + 1  # at extrema
    plotFit(xvar, yvar, std, finfit$par, FUNmod, finfit$iout, ctrl$bg,
            smps[ctrl$ismp, yvar], smps[ctrl$ismp, "Flag"], trimval, trimext,
            ptcol, rugcols,
            main = main, ylab = yvar, xlab = "Concentration", ...)
  } else {
    plotFit(xvar, yvar, std, iout = finfit$out, bg = ctrl$bg, smp = ctrl$smp,
            ptcol = ptcol, rugcols = rugcols,
            main = main, ylab = yvar, xlab = "Concentration", ...)
  }
  dev.off()
  return(list(smps = smps, fitflag = finfit$flag))
}
