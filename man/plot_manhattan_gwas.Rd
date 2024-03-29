% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_manhattan_gwas}
\alias{plot_manhattan_gwas}
\title{Plot a Manhattan plot for GWAS results}
\usage{
plot_manhattan_gwas(
  df,
  chr = "CHR",
  minuslog10_p = "P",
  text = "SNP",
  order_idx = "ORDER_IDX",
  geom_point_args = list(alpha = 1, size = 2),
  labs_args = list(x = "Chromosome", y = "-log10(p value)"),
  col_manual = c("blue4", "orange3")
)
}
\arguments{
\item{df}{A data frame.}

\item{chr}{Name of chromosome column. Should be type character or factor.}

\item{minuslog10_p}{Name of -log10(p value) column. Should be on -log10 scale.}

\item{text}{Name of column with text annotation.}

\item{order_idx}{Name of column with order to arrange SNPs. Should be type
integer.}

\item{geom_point_args}{A named list of arguments to be passed to
\code{\link[ggplot2:geom_point]{ggplot2::geom_point()}}.}

\item{labs_args}{List of arguments to be passed to \code{\link[ggplot2:labs]{ggplot2::labs()}}.}

\item{col_manual}{Character vector of colours}
}
\value{
A data frame
}
\description{
Plot a Manhattan plot for GWAS results
}
