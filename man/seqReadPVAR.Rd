\name{seqReadPVAR}
\alias{seqReadPVAR}
\title{Read PLINK2 pvar file}
\description{
    Read PLINK2 pvar file for variants
}

\usage{
seqReadPVAR(pvar, sel=NULL)
}

\arguments{
    \item{pvar}{a file name of a pvar file (from \link[pgenlibr]{NewPvar}),
        or a pvar object, which can be queried for variant IDs and allele
        codes}
    \item{sel}{\code{NULL}, a logical vector or a numeric vector for
        specifying the variants; \code{NULL} for including all variants}
}

\value{
    Return a data frame with the columns \code{chrom}, \code{pos},
\code{allele} and \code{rsid}.
}

\references{
    \url{https://www.cog-genomics.org/plink/2.0/}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqPGEN2GDS}}, \code{\link{seqGDS2PGEN}}
}

\examples{
pvar_fn <- system.file("extdata", "plink2_gen.pvar", package="pgen2gds")

head(seqReadPVAR(pvar_fn))
}

\keyword{GDS}
\keyword{PLINK2}
\keyword{utilities}
