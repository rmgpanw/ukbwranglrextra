% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ukb_synthetic_dataset.R
\name{get_ukb_synthetic_data}
\alias{get_ukb_synthetic_data}
\title{Download the UKB synthetic dataset}
\usage{
get_ukb_synthetic_data(download_dir, dataset = "tabular.md5", timeout = 600)
}
\arguments{
\item{download_dir}{Directory to download files to.}

\item{dataset}{Character. The name of one of the \code{.md5} files on the \href{https://biobank.ndph.ox.ac.uk/ukb/exinfo.cgi?src=UKB_Synthetic_Dataset.html}{UKB synthetic data webpage}.
All}

\item{timeout}{Integer. Passed to \code{options()} to determine how long
\code{download.file()} attempts to download for.}
}
\value{
Invisibly returns a list with downloaded file paths and errors for
failed downloads.
}
\description{
Download the UKB synthetic dataset
}
\examples{
\dontrun{
 get_ukb_synthetic_data(download_dir = tempdir())
}
}
