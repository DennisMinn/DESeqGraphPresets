% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplots.R
\name{mk_volcanoPlot}
\alias{mk_volcanoPlot}
\title{Make Volcano Plot}
\usage{
mk_volcanoPlot(
  format,
  res,
  pCutoff,
  lfcCutoff,
  ntop = 20,
  geneList,
  unit,
  title,
  subtitle,
  xlim,
  ylim
)
}
\arguments{
\item{format}{Specifies how to subset the \code{res}. \code{"list"} will use desired genes symbols while
\code{"padj"} will use the lowest \code{'padj'} values .}

\item{res}{Contrast data produced by \code{DESeq2::res()} or \code{DESeq2::lfcShrink()}}

\item{lfcCutoff}{Log2(foldchange) threshhold}

\item{ntop}{Limits the number of genes included in \code{format = "padj"}.}

\item{geneList}{Input, desired list of genes}

\item{unit}{String to denote the comparison used for \code{res}}

\item{title}{String for title for graph}

\item{subtitle}{Input subtitle for graph}

\item{xlim}{Defines the x-axis boundaries}

\item{ylim}{Defines the y-axis boundaries}

\item{pCutOff}{P-adj threshhold.}
}
\description{
Makes a volcano plot with emphasis on genes with the lowest p-adj values or specific genes. Values below either
the pCutoff or lfcCutoff will be denoted in grey, otherwise red.
}
\examples{
## Example 1: Using "p-adj"

# Volcano Plot with lowest p-adj genes, p-adj cutoff of 0.05, and log2(foldchange) of 1 (default 20 genes)
> mk_volcanoPlot(format = "padj", res = res, pCutoff = 0.05, lfcCutoff = 1) 

# Specifying lowest 20 p-adj genes
> mk_volcanoPlot(format = "padj", res = res, pCutoff = 0.05, lfcCutoff = 1, ntop = 10) 

# Specifying bounds for x-axis [-5,5] and y-axis [-9,10]
> mk_volcanoPlot(format = "padj", res = res, pCutoff = 0.05, lfcCutoff = 1, xlim = c(-5,5), ylim = c(-9, 10)) 

## Example 2: Using "list"
> selectLab <- c("TSPAN6", "TNMD", "DPM1", "SCYL3", "C1orf112", "FGR")
> mk_volcanoPlot(format = "list", res = res, pCutoff = 0.05, lfcCutoff = 1, geneList = selectLab) 

## Example 3: Exporting graph to pdf
> pdf("volcanoPlot.pdf")
> mk_volcanoPlot(format = "padj", res = res, pCutoff = 0.05, lfcCutoff = 1) 
> dev.off()
}
