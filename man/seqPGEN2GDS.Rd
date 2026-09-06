\name{seqPGEN2GDS}
\alias{seqPGEN2GDS}
\title{Reformat PLINK2 PGEN files}
\description{
    Reformats PLINK2 pgen files to GDS format.
}

\usage{
seqPGEN2GDS(pgen.fn, pvar.fn=NULL, psam.fn=NULL, out.gdsfn,
    compress.geno="LZMA_RA", compress.annot="LZMA_RA", variant.sel=NULL,
    sample.sel=NULL, start=1L, count=NA_integer_,
    ignore.chr.prefix=c("chr", "0"), reference=NULL, optimize=TRUE,
    digest=TRUE, parallel=FALSE, balancing=TRUE, verbose=TRUE)
}

\arguments{
    \item{pgen.fn}{a file name for the pgen file}
    \item{pvar.fn}{a file name for the pvar file, or \code{NULL} to use
        the default}
    \item{psam.fn}{a file name for the psam file, or \code{NULL} to use
        the default}
    \item{out.gdsfn}{the file name of output GDS file}
    \item{compress.geno}{the compression method for "genotype"; optional
        values are defined in the function \code{add.gdsn}}
    \item{compress.annot}{the compression method for the GDS variables,
        except "genotype"; optional values are defined in the function
        \code{add.gdsn}}
    \item{variant.sel}{\code{NULL} for no variant selection, a logical vector
        or a numeric vector to specify the variant selection}
    \item{sample.sel}{\code{NULL} for no sample selection, a logical vector
        or a numeric vector to specify the sample selection}
    \item{start}{the starting variant if importing part of the pgen file}
    \item{count}{the maximum count of variant if importing part of the pgen
        file, \code{NA_integer_} or any non-positive value indicates importing
        to the end}
    \item{ignore.chr.prefix}{a vector of character, indicating the prefix of
        chromosome which should be ignored, e.g., \code{"chr"}; it is not
        case-sensitive}
    \item{reference}{genome reference, like "GRCh37", "GRCh38"; it is not
        specified if \code{reference=NULL}}
    \item{optimize}{if \code{TRUE}, optimize the access efficiency by calling
        \code{\link[gdsfmt]{cleanup.gds}}}
    \item{digest}{a logical value (TRUE/FALSE) or a character (e.g., "md5");
        add hash codes to the GDS file if TRUE or a digest algorithm is
        specified}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (parallel
        processing), a numeric value indicating the number of cores, or a
        cluster object for parallel processing; \code{parallel} is passed to
        the argument \code{cl} in \code{\link[SeqArray]{seqParallel}}, see
        \code{\link[SeqArray]{seqParallel}} for more details}
    \item{balancing}{whether to perform workload balancing or not, only
        applicable when multiple cores are used; if \code{NA}, use \code{TRUE}
        as a default until \code{getOption("seqarray.balancing")} is set and
        not \code{TRUE}}
    \item{verbose}{if \code{TRUE}, show information}
}

\details{
    The sample IDs are taken from the column \code{IID} of the psam file. The
psam columns \code{FID}, \code{PAT}, \code{MAT} and \code{SEX} are stored in
the GDS variables \code{sample.annotation/family},
\code{sample.annotation/father}, \code{sample.annotation/mother} and
\code{sample.annotation/sex} (coded as "M", "F" or "") if they exist, and the
first phenotype column is stored in \code{sample.annotation/phenotype} (the
same variables as created by \code{\link[SeqArray]{seqBED2GDS}}). Dosages and
phase information in the pgen file are not imported.
}

\value{
    Return the file name of SeqArray GDS file with an absolute path.
}

\references{
    \url{https://www.cog-genomics.org/plink/2.0/}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqGDS2PGEN}}, \code{\link{seqReadPVAR}}
}

\examples{
pgen_fn <- system.file("extdata", "plink2_gen.pgen", package="pgen2gds")

seqPGEN2GDS(pgen_fn, out.gdsfn="test.gds")

# delete the temporary file
unlink("test.gds", force=TRUE)
}

\keyword{GDS}
\keyword{PLINK2}
\keyword{utilities}
