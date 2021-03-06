% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap.R
\name{mk_heatmap}
\alias{mk_heatmap}
\title{Make Heatmap}
\usage{
mk_heatmap(mat, df, sampleNames, title)
}
\arguments{
\item{mat}{matrix containing the expressions of desired genes, look at \code{mk_mat} for more detail}

\item{df}{data frame formatted from your design table, look at \code{mk_df} for more detail}

\item{sampleNames}{column vector containing the names of your samples}

\item{title}{optional argument containing the title of your graph}
}
\value{
Returns a heatmap either in your plots window, alternatively can be directed into a file with \code{pdf(file = "file_name")} and \code{dev.off()}
}
\description{
Creates a heatmap using \code{pheatmap} package with preset values
}
\examples{
## Example 1

> sampleNames <- paste("sample", 1:18, sep = ".")
> mk_heatmap(mat = mat, df = df, sampleNames = sampleNames, title = "heatmap")

## Example 2

# Directs the output into a pdf
> pdf(file = "Heatmap.pdf")
> sampleNames <- paste("sample", 1:18, sep = ".")
> mk_heatmap(mat = mat, df = df, sampleNames, "heatmap")
> dev.off()

}
