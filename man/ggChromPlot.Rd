% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaboVis.R
\name{ggChromPlot}
\alias{ggChromPlot}
\title{Chromatogram plot}
\usage{
ggChromPlot(
  object,
  filenum = NA,
  mz = NA,
  ppm = NA,
  rtint = NA,
  pheno_var = 1,
  chromtype = "max",
  logscale = TRUE,
  interactive = TRUE
)
}
\arguments{
\item{object}{A OnDiskXCMSnExp or XCMSnExp object}

\item{filenum}{An optional numeric vector indicating the files to plot}

\item{mz}{A numeric indicating the mz value to plot (requires the ppm param) or a numeric two with the mz limits to plot}

\item{ppm}{An optional numeric indicating the ppm error for the mz param}

\item{rtint}{An optional length two numeric vector indicating the retention time interval to plot}

\item{pheno_var}{A numeric or string indicating the variable from the phenodata
to use as file groups.}

\item{chromtype}{The type of chromatogram to create. "max" for the maximum value and "sum" for the sum of intensities.}

\item{logscale}{A boolean to apply a log10 transformation of the intensities}

\item{interactive}{A boolean indicating if the plot will be converted to a plotly}
}
\value{
A ggplot2 or plotly line plot of the chromatogram.
}
\description{
Creates a ggplot2 chromatogram or EIC from an XCMS object
}
\author{
Jordi Rofes Herrera
}
