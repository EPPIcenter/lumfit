#' Extract Standards
#'
#' Separates specified serial dilutions for fitting, background, and samples
#'
#' details
#'
#' @param MFI data frame with variable "Sample" that contains information about
#'   standard concentrations in "1/dilution" format (e.g. "1/200").
#' @param dilut standard dilutions to use for standard curve fitting. If
#'   \code{NULL}, all the dilutions indicated in a "Sample" variable are used.
#' @inheritParams processLum
#'
#' @return A list with standards for fitting, background values, sample values,
#'   indices for samples, and a data frame - for a specified antigen.
#'
#' @export

extractStd <- function(MFI, stdstr, bgstr, dilut, smpdil, antigen, yvar) {
  mfi <- suppressWarnings(as.numeric(MFI[, grep(antigen, colnames(MFI))]))
  istd <- grep("1/", MFI$Sample)
  if (length(istd) == 0) {    # handled by processLum(); for standalone use only
    stop('Please provide standard concentrations in "1/dilution" format,
         e.g "1/200"')
  }
  if (stdstr != "") {
    istd <- istd[grep(tolower(stdstr), tolower(MFI$Sample)[istd])]
    if (length(istd) == 0) {  # handled by processLum(); for standalone use only
      # no warning, back to original istd
      istd <- grep("1/", MFI$Sample)
    }
  }
  sname <- MFI$Sample[istd]
  locdil <- regexpr("1/", sname)
  sdilut <- as.numeric(substr(sname, locdil + 2, nchar(sname)))
  if (!is.null(dilut)) {
    idil <- sdilut %in% dilut
  } else {
    idil <- seq_along(sdilut)
  }
  std <- data.frame(Dil = sdilut[idil], Conc = 1/sdilut[idil],
                    MFI = mfi[istd[idil]])
  ibg  <- grep(tolower(bgstr),  tolower(MFI$Sample))
  #  ineg <- grep(tolower(negstr), tolower(MFI$Sample))
  ismp <- (1:nrow(MFI))[-c(istd, ibg)]
  std$logConc <- log(std$Conc)
  bg <- mfi[ibg]
  std$logMFI    <- log(std$MFI)
  std <- std[order(std$Conc), ]  # sorted by increasing concentration

  dfout                <- MFI[, 1:2]
  dfout$logMFI         <- log(mfi)
  dfout$Dilution       <- smpdil
  dfout$Dilution[ibg]  <- 0
  dfout$Dilution[istd] <- sdilut
  return(list(std = std, bg = bg, out = dfout, smp = dfout[ismp, yvar],
              ismp = ismp))
}

#' Read Raw File
#'
#' Extracts the specified data type from the file
#'
#' details
#'
#' @inherit processLum params
#'
#' @return A data frame containing data for all the antigens in the file.
#'
#' @export

read_data <- function(fname, dtype = "Median", nwells = NULL, nsep = 2,
                      ncolmax = 105) {
  all <- read.csv(fname, header = FALSE, blank.lines.skip = FALSE,
                  as.is = TRUE, col.names = 1:ncolmax)
  itype <- grep("Data", all[[1]])
  types <- all[itype, 2]
  if (is.null(nwells)) {
    nwells <- itype[2] - itype[1] - nsep - 1
  }
  nskip <- itype[tolower(types) %in% tolower(dtype)]
  if (length(nskip) == 0) {                       # what if > 1 (can't happen?)
    stop("Specified data type not found")
  }
  read.csv(fname, as.is = TRUE, skip = nskip, nrows = nwells)
}
