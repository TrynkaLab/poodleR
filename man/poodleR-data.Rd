\name{poodleR-data}
\alias{A}
\alias{b}
\alias{b_estimate_dt}
\alias{bam_readcount}
\alias{vcf}
\title{Example data objects included with poodleR}
\description{
Small example objects used in documentation, tests, and tutorials.
}
\details{
\itemize{
\item \code{A}: genotype dosage table with donor columns and an \code{rn} SNP identifier.
\item \code{b}: pooled allele-frequency table with a \code{b_estimate} column.
\item \code{b_estimate_dt}: filtered allele-frequency table aligned to \code{A}.
\item \code{bam_readcount}: the first ten columns of a bam-readcount import.
\item \code{vcf}: donor genotype data stored as a \code{vcfR} object.
}
}
\format{
These objects are stored in \code{data/} and can be loaded with
\code{data(<name>, package = "poodleR")}.
}
