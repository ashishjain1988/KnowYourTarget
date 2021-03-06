% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGeneEffect.R
\name{getGeneEffectsInDepMap}
\alias{getGeneEffectsInDepMap}
\title{Interactive plot to show the effects of genes on cancer proliferation}
\usage{
getGeneEffectsInDepMap(
  genes = NULL,
  depMapVersion = "",
  groupBySubtypes = FALSE
)
}
\arguments{
\item{genes}{A list of genes (Hugo symbol)}

\item{depMapVersion}{DepMap data version to be used}

\item{groupBySubtypes}{Group the cancer cell lines from cancer subtypes, default=FALSE}
}
\value{
A list containing plots for individual genes. For each gene, we have three plots showing
the DepMap CRISPR-Cas9 data in the form of histogram, boxplot, and combined plot.
}
\description{
This function creates the interactive plot to show the effects of
genes on cancer proliferation using the DepMap CRISPR-Cas9 data.
different figures and plots explaining the MAF dataset.
}
\author{
Ashish Jain
}
