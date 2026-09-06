# ===========================================================================
#
# gds2pgen.r: Format Conversion from GDS to PLINK2 PGEN
#
# Copyright (C) 2026    Xiuwen Zheng (zhengx@u.washington.edu)
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License Version 3 as
# published by the Free Software Foundation.
#
# pgen2gds is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with pgen2gds.
# If not, see <http://www.gnu.org/licenses/>.


#############################################################
# Internal functions
#

# write a data.frame as a tab-delimited PLINK2 text file with a '#' header
.write_plink_table <- function(df, fn)
{
    con <- file(fn, "wt")
    on.exit(close(con))
    writeLines(paste0("#", paste(names(df), collapse="\t")), con)
    write.table(df, con, quote=FALSE, sep="\t", row.names=FALSE,
        col.names=FALSE, na="NA")
    invisible()
}

# PLINK2 sex codes: 1 = male, 2 = female, NA = unknown
.psam_sex <- function(x)
{
    s <- rep(NA_character_, length(x))
    if (is.numeric(x))
    {
        s[x %in% 1L] <- "1"
        s[x %in% 2L] <- "2"
    } else {
        x <- toupper(trimws(as.character(x)))
        s[x %in% c("M", "MALE", "1")] <- "1"
        s[x %in% c("F", "FEMALE", "2")] <- "2"
    }
    s
}

# a character vector with missing values replaced by 'na'
.psam_str <- function(x, na)
{
    x <- as.character(x)
    x[is.na(x) | x==""] <- na
    x
}

# build the data.frame for a psam file from sample.id and sample.annotation
.make_psam <- function(gdsfile, sample.id)
{
    psam <- data.frame(IID=sample.id, stringsAsFactors=FALSE)
    nm <- "sample.annotation/family"
    if (exist.gdsn(gdsfile, nm))
    {
        psam <- cbind(FID=.psam_str(seqGetData(gdsfile, nm), "0"), psam,
            stringsAsFactors=FALSE)
    }
    nm <- "sample.annotation/father"
    if (exist.gdsn(gdsfile, nm))
        psam$PAT <- .psam_str(seqGetData(gdsfile, nm), "0")
    nm <- "sample.annotation/mother"
    if (exist.gdsn(gdsfile, nm))
        psam$MAT <- .psam_str(seqGetData(gdsfile, nm), "0")
    nm <- "sample.annotation/sex"
    if (exist.gdsn(gdsfile, nm))
        psam$SEX <- .psam_sex(seqGetData(gdsfile, nm))
    nm <- "sample.annotation/phenotype"
    if (exist.gdsn(gdsfile, nm))
        psam$PHENO1 <- seqGetData(gdsfile, nm)
    psam
}

# variant IDs for the pvar file
.pvar_id <- function(gdsfile, write.rsid)
{
    if (write.rsid == "chr_pos_ref_alt")
    {
        rsid <- gsub(":|,", "_", seqGetData(gdsfile, "$chrom_pos_allele"))
    } else {
        rsid <- seqGetData(gdsfile, "annotation/id")
        rsid[is.na(rsid)] <- ""
        if (write.rsid == "auto")
        {
            x <- rsid %in% c("", ".")
            if (any(x))
            {
                rsid[x] <- gsub(":|,", "_",
                    seqGetData(gdsfile, "$chrom_pos_allele")[x])
            }
        }
    }
    rsid[is.na(rsid) | rsid==""] <- "."
    rsid
}


#############################################################
# Format conversion from GDS to PGEN
#
seqGDS2PGEN <- function(gdsfile, out.fn,
    write.rsid=c("auto", "annot_id", "chr_pos_ref_alt"),
    half.call=c("missing", "haploid"), use.phase=TRUE, chr_prefix="",
    verbose=TRUE)
{
    # check
    stopifnot(is.character(gdsfile) | inherits(gdsfile, "SeqVarGDSClass"))
    stopifnot(is.character(out.fn), length(out.fn)==1L, !is.na(out.fn))
    write.rsid <- match.arg(write.rsid)
    half.call <- match.arg(half.call)
    stopifnot(is.logical(use.phase), length(use.phase)==1L, !is.na(use.phase))
    stopifnot(is.character(chr_prefix), length(chr_prefix)==1L,
        !is.na(chr_prefix))
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (verbose)
    {
        .cat("##< ", .tm())
        .cat("SeqArray GDS to PLINK2 PGEN:")
    }

    # open the GDS file
    if (is.character(gdsfile))
    {
        stopifnot(length(gdsfile)==1L)
        if (verbose) .cat("    open ", sQuote(basename(gdsfile)))
        gdsfile <- seqOpen(gdsfile, allow.duplicate=TRUE)
        on.exit(seqClose(gdsfile))
    }

    # dimensions (respecting the current sample and variant filters)
    dm <- seqSummary(gdsfile, "genotype", check="none", verbose=FALSE)$seldim
    ploidy <- dm[1L]; nsamp <- dm[2L]; nvar <- dm[3L]
    if (is.na(ploidy))
        stop("'ploidy' is not known.")
    else if (ploidy != 2L)
        stop("'ploidy' should be 2 for diploidy.")
    if (nsamp <= 0L) stop("There is no selected sample.")
    if (nvar <= 0L) stop("There is no selected variant.")
    if (verbose)
    {
        .cat("    # of samples: ", .pretty(nsamp))
        .cat("    # of variants: ", .pretty(nvar))
        .cat("    [Output]")
    }

    # psam file
    sample.id <- seqGetData(gdsfile, "sample.id")
    psamfn <- paste0(out.fn, ".psam")
    if (verbose) .cat("    PSAM: ", psamfn)
    .write_plink_table(.make_psam(gdsfile, sample.id), psamfn)

    # pvar file
    alt <- seqGetData(gdsfile, "$alt")
    alt[is.na(alt) | alt==""] <- "."
    pvar <- data.frame(
        CHROM = paste0(chr_prefix, seqGetData(gdsfile, "chromosome")),
        POS   = seqGetData(gdsfile, "position"),
        ID    = .pvar_id(gdsfile, write.rsid),
        REF   = seqGetData(gdsfile, "$ref"),
        ALT   = alt,
        stringsAsFactors=FALSE)
    remove(alt)
    pvarfn <- paste0(out.fn, ".pvar")
    if (verbose) .cat("    PVAR: ", pvarfn)
    .write_plink_table(pvar, pvarfn)
    remove(pvar)

    # pgen file
    pgenfn <- paste0(out.fn, ".pgen")
    if (verbose) .cat("    PGEN: ", pgenfn)
    if (use.phase && !exist.gdsn(gdsfile, "phase/data"))
        use.phase <- FALSE
    num_allele <- seqGetData(gdsfile, "$num_allele")
    max_allele <- max(num_allele, na.rm=TRUE)
    remove(num_allele)
    if (max_allele > 255L)
        stop("PLINK2 supports at most 255 alleles per variant.")
    ptr <- .Call(SEQ_PGEN_Writer_Open, pgenfn, nsamp, nvar, use.phase,
        as.integer(max_allele))
    ok <- FALSE
    on.exit({
        if (!ok) .Call(SEQ_PGEN_Writer_Close, ptr, FALSE)
    }, add=TRUE)

    # block size: about 32MB of raw genotypes and phase per block
    bsize <- as.integer(max(16L, min(4096L, 33554432 %/% (3L * nsamp))))
    hc <- if (half.call == "haploid") 1L else 0L
    if (use.phase)
    {
        FUN <- function(x)
        {
            .Call(SEQ_PGEN_Writer_Block, ptr, x[[1L]], x[[2L]], x[[3L]], hc)
            NULL
        }
        vars <- c("genotype", "phase", "$num_allele")
    } else {
        FUN <- function(x)
        {
            .Call(SEQ_PGEN_Writer_Block, ptr, x[[1L]], NULL, x[[2L]], hc)
            NULL
        }
        vars <- c("genotype", "$num_allele")
    }
    seqBlockApply(gdsfile, vars, FUN, margin="by.variant", as.is="none",
        .useraw=TRUE, bsize=bsize, .progress=verbose)

    # finalize the header and index
    cnt <- .Call(SEQ_PGEN_Writer_Close, ptr, TRUE)
    ok <- TRUE
    if (verbose)
    {
        .cat("        size: ", .pretty_size(file.size(pgenfn)))
        .cat("        variant records: ", .pretty(cnt[2L]), " dense, ",
            .pretty(cnt[3L]), " 1-bit, ", .pretty(cnt[4L]), " difflist")
        if (cnt[5L] > 0L)
            .cat("        multiallelic patches: ", .pretty(cnt[5L]))
        if (cnt[6L] > 0L)
            .cat("        phased variants: ", .pretty(cnt[6L]))
        .cat("Done.\n##> ", .tm())
    }

    # output
    invisible(normalizePath(c(pgenfn, pvarfn, psamfn)))
}
