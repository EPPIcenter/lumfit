% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitStd.R
\name{fitStd}
\alias{fitStd}
\title{Fit Standards}
\usage{
fitStd(std, xvar, yvar, model = "sigmoid", Alow = NULL, asym = TRUE,
  interactive = TRUE, monot.prompt = FALSE, rm.before = FALSE,
  rm.after = interactive, maxrm = 2, set.bounds = FALSE,
  overwrite.bounds = FALSE, bg = NULL, smp = NULL,
  optmethod = "Nelder-Mead", maxit = 5000, info = "", ifix = NULL,
  rugcols, ptcol, ...)
}
\arguments{
\item{std}{data frame of standard dilutions}

\item{xvar}{character strings for the variables used to fit a standard
curve.}

\item{yvar}{character strings for the variables used to fit a standard
curve.}

\item{model}{the model to be fit.}

\item{Alow}{lower asymptote for the sigmoid model. If \code{NULL}, the lower
asymptote will be estimated, adding an extra parameter to the model. To fix
the asymptote at the level of background, specify \code{"bg"}. Numeric
value of \code{Alow} will force the asymptote to be fixed at the provided
level.}

\item{asym}{if \code{TRUE}, asymmetry in the fit is allowed, adding an extra
parameter to the model.}

\item{interactive}{logical value. If \code{TRUE}, the user is prompted to
evaluate the standards (and/or the fit) and asked to remove outliers if
needed. \code{TRUE} value takes precedence over \code{rm.before} and
\code{rm.after}: if both are \code{FALSE}, \code{rm.after} is reset to
\code{TRUE}.}

\item{monot.prompt}{if \code{TRUE}, the user is prompted to evaluate the
standards and possibly remove outliers if the standards are not monotonic
(increasing). \code{FALSE} value is ignored if \code{interactive} is
\code{TRUE}.}

\item{rm.before}{logical value indicating if potential outliers should be
removed before the model is fitted. Ignored if \code{interactive} is
\code{FALSE}.}

\item{rm.after}{logical value indicating if potential outliers should be
removed after the model is fitted. Ignored if \code{interactive} is
\code{FALSE}.}

\item{maxrm}{maximum number of outliers to remove.}

\item{set.bounds}{if \code{TRUE}, the user is given the option to manually
set the bound that is not set automatically. In that case, the prompt
appears even if \code{interactive} is \code{FALSE}.}

\item{overwrite.bounds}{logical value indicating the option to overwrite
automatically set bounds. Ignored if \code{interactive} is \code{FALSE}.}

\item{bg}{background values}

\item{smp}{sample values}

\item{optmethod}{method to be used in optimization.}

\item{maxit}{maximum number of iterations in optimization.}

\item{info}{information about a particular run for warning messages.}

\item{ifix}{sorted integer vector of length 3 with indices of standards to be
used for getting starting values for optimization.}

\item{rugcols}{vector of three colors for the rugplot, which indicates sample
values (inside the bounds, between the bounds and extrema, and beyond
extrema).}

\item{ptcol}{color of the standard points on the plot.}

\item{...}{further graphical parameters.}
}
\value{
A list containing parameters of the fit and bounds of the fit (named
  vectors), as well as indices of removed points (if any) and flags r.
}
\description{
Fits a specified function to MFI summary of serial dilutions. Optionally an
interactive procedure that allows to remove outliers, evaluate resulting
fits, perform revisions, and record a message regarding the fit.
}
\details{
to be added
}
\examples{

}
