# ===========================================================================
#
# pgen2gds.r: Format Conversion from PLINK2 PGEN to GDS
#
# Copyright (C) 2025-2026    Xiuwen Zheng (zhengx@u.washington.edu)
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

.cat <- function(...) message(paste0(...))

.tm <- function() strftime(Sys.time(), "%Y-%m-%d %H:%M:%S")

.pretty <- function(x) prettyNum(x, big.mark=",", scientific=FALSE)

.pretty_size <- function(x)
{
    stopifnot(is.numeric(x), length(x)==1L)
    if (is.na(x) || !is.finite(x))
        "NA"
    else if (x >= 1024^4L)
        sprintf("%.1fT", x/1024^4L)
    else if (x >= 1024^3L)
        sprintf("%.1fG", x/1024^3L)
    else if (x >= 1024^2L)
        sprintf("%.1fM", x/1024^2L)
    else if (x >= 1024)
        sprintf("%.1fK", x/1024)
    else
        sprintf("%gB", x)
}

# append to a gds node using the data return from a user-defined function
.append_fc_gds <- function(nd, start, count, type, fc, fc2=identity)
{
    last <- start + count
    i <- start
    while (i < last)
    {
        cnt <- min(i+10000L, last) - i
        s <- vapply(seq.int(i, length.out=cnt), fc, type)
        append.gdsn(nd, fc2(s))
        i <- i + cnt
    }
    readmode.gdsn(nd)
    invisible()
}

# read a psam file
.read_psam <- function(psam.fn)
{
    s <- readLines(psam.fn)
    # remove ## line(s)
    x <- grepl("^##", s)
    if (any(x)) s <- s[!x]
    # have a header or not
    use_header <- any(grepl("^#", s))
    # \t or whitespace
    if (any(grepl("\t", s)))
    {
        fam <- read.table(text=s, header=use_header, sep="\t", comment.char="",
            stringsAsFactors=FALSE)
    } else {
        fam <- read.table(text=s, header=use_header, sep=" ", comment.char="",
            stringsAsFactors=FALSE)
    }
    # check the header
    s <- names(fam)
    if (use_header)
    {
        s[s=="X.IID"] <- "IID"; s[s=="X.FID"] <- "FID"
        names(fam) <- s
    } else {
        nm <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")
        i <- seq_len(min(ncol(fam), length(nm)))
        s[i] <- nm[i]
        names(fam) <- s
    }
    fam
}


#############################################################
# Read PLINK2 pvar file and return a data.frame
#
seqReadPVAR <- function(pvar, sel=NULL)
{
    if (is.character(pvar))
    {
        pvar <- NewPvar(pvar)
        on.exit(ClosePvar(pvar))
    } else if (inherits(pvar, "pvar") && is(pvar$pvar, "externalptr"))
    {
        # valid pvar object, use as-is
    } else {
        stop("'pvar' should be a file name or an object ",
            "returned from pgenlibr::NewPvar().")
    }
    # selection
    if (is.null(sel))
        sel <- seq_len(GetVariantCt(pvar))
    else if (is.logical(sel))
        sel <- which(sel)
    else if (!is.numeric(sel))
        stop("'sel' should be NULL, a logical vector or a numeric vector.")

    # read
    data.frame(
        chrom  = vapply(sel, function(i) GetVariantChrom(pvar, i), ""),
        pos    = vapply(sel, function(i) GetVariantPos(pvar, i), 0L),
        allele = vapply(sel, function(i)
            paste(vapply(seq_len(GetAlleleCt(pvar, i)), function(j)
                GetAlleleCode(pvar, i, j), ""), collapse=","), ""),
        rsid   = vapply(sel, function(i) GetVariantId(pvar, i), "")
    )
}


#############################################################
# Internal: validate and process variant selection
#
.validate_variant_sel <- function(variant.sel, nvar)
{
    if (is.logical(variant.sel))
    {
        if (length(variant.sel) != nvar)
            stop("'length(variant.sel)' should be ", nvar, ".")
        variant.sel <- which(variant.sel)
        if (length(variant.sel) <= 0L)
            stop("There is no selected variant in 'variant.sel'.")
    } else if (is.numeric(variant.sel))
    {
        if (!is.integer(variant.sel))
            variant.sel <- as.integer(variant.sel)
        if (anyNA(variant.sel))
            stop("'variant.sel' should not have NA.")
        if (length(variant.sel) <= 0L)
            stop("There is no selected variant in 'variant.sel'.")
        if (anyDuplicated(variant.sel))
            stop("'variant.sel' should not have any duplicates.")
        variant.sel <- sort(variant.sel)
        if (variant.sel[1L] < 1L || variant.sel[length(variant.sel)] > nvar)
            stop("'variant.sel' should be between 1 and ", nvar, ".")
    } else {
        stop("'variant.sel' should be NULL or a logical/numeric vector.")
    }
    variant.sel
}

#############################################################
# Internal: validate and process sample selection
#
.validate_sample_sel <- function(sample.sel, nsamp)
{
    if (is.logical(sample.sel))
    {
        if (length(sample.sel) != nsamp)
            stop("'length(sample.sel)' should be ", nsamp, ".")
        sample.sel <- which(sample.sel)
    } else if (is.numeric(sample.sel))
    {
        if (!is.integer(sample.sel))
            sample.sel <- as.integer(sample.sel)
        if (anyNA(sample.sel))
            stop("'sample.sel' should not have NA.")
        if (anyDuplicated(sample.sel))
            stop("'sample.sel' should not have any duplicates.")
        sample.sel <- sort(sample.sel)
        if (length(sample.sel) > 0L)
        {
            if (sample.sel[1L] < 1L ||
                sample.sel[length(sample.sel)] > nsamp)
            {
                stop("'sample.sel' should be between 1 and ", nsamp, ".")
            }
        }
    } else {
        stop("'sample.sel' should be NULL or a logical/numeric vector.")
    }
    if (length(sample.sel) <= 0L)
        stop("There is no selected sample in 'sample.sel'.")
    sample.sel
}

#############################################################
# Internal: write variant metadata (chromosome, position, allele, id) to GDS
#
.write_variant_annot <- function(dstfile, pvar, variant.sel, start, count,
    compress.annot, ignore.chr.prefix, digest, verbose)
{
    # add chromosome
    if (verbose) message("    chromosome  ", appendLF=FALSE)
    strip_chr <- length(ignore.chr.prefix) > 0L
    if (strip_chr)
    {
        escaped <- gsub("([.\\|^$*+?(){}\\[\\]])", "\\\\\\1",
            ignore.chr.prefix)
        ignore.chr.prefix <- paste0("^(",
            paste(escaped, collapse="|"), ")")
    }
    n <- add.gdsn(dstfile, "chromosome", storage="string",
        compress=compress.annot)
    fc2 <- if (strip_chr)
        function(s) {
            r <- sub(ignore.chr.prefix, "", s, ignore.case=TRUE)
            ifelse(nchar(r) == 0L, s, r)
        } else identity
    if (is.null(variant.sel))
    {
        .append_fc_gds(n, start, count, "",
            function(i) GetVariantChrom(pvar, i), fc2)
    } else {
        .append_fc_gds(n, start, count, "",
            function(i) GetVariantChrom(pvar, variant.sel[i]), fc2)
    }
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)
    # RLE-coded chromosome
    SeqArray:::.optim_chrom(dstfile)

    # add position
    if (verbose) message("    position  ", appendLF=FALSE)
    n <- add.gdsn(dstfile, "position", storage="int32",
        compress=compress.annot)
    if (is.null(variant.sel))
    {
        .append_fc_gds(n, start, count, 0L, function(i) GetVariantPos(pvar, i))
    } else {
        .append_fc_gds(n, start, count, 0L,
            function(i) GetVariantPos(pvar, variant.sel[i]))
    }
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)

    # add allele
    if (verbose) message("    allele  ", appendLF=FALSE)
    n <- add.gdsn(dstfile, "allele", storage="string",
        compress=compress.annot)
    if (is.null(variant.sel))
    {
        .append_fc_gds(n, start, count, "", function(i) {
            paste(vapply(seq_len(GetAlleleCt(pvar, i)), function(j)
                GetAlleleCode(pvar, i, j), ""), collapse=",")
        })
    } else {
        .append_fc_gds(n, start, count, "", function(i) {
            i <- variant.sel[i]
            paste(vapply(seq_len(GetAlleleCt(pvar, i)), function(j)
                GetAlleleCode(pvar, i, j), ""), collapse=",")
        })
    }
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)

    # add annotation/id (rsid)
    if (verbose) message("    annotation/id  ", appendLF=FALSE)
    nann <- addfolder.gdsn(dstfile, "annotation")
    n <- add.gdsn(nann, "id", storage="string", compress=compress.annot)
    if (is.null(variant.sel))
    {
        .append_fc_gds(n, start, count, "", function(i) GetVariantId(pvar, i))
    } else {
        .append_fc_gds(n, start, count, "",
            function(i) GetVariantId(pvar, variant.sel[i]))
    }
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)
    sync.gds(dstfile)
    invisible()
}


#############################################################
# Format conversion from PGEN to GDS
#
seqPGEN2GDS <- function(pgen.fn, pvar.fn=NULL, psam.fn=NULL, out.gdsfn,
    compress.geno="LZMA_RA", compress.annot="LZMA_RA",
    variant.sel=NULL, sample.sel=NULL, start=1L, count=NA_integer_,
    ignore.chr.prefix=c("chr", "0"), reference=NULL, optimize=TRUE,
    digest=TRUE, parallel=FALSE, balancing=TRUE, verbose=TRUE)
{
    # check
    stopifnot(is.character(pgen.fn), length(pgen.fn)==1L)
    if (is.null(pvar.fn) && is.null(psam.fn))
    {
        fn <- gsub("\\.pgen$", "", pgen.fn, ignore.case=TRUE)
        # set pvar file
        pvar.fn <- paste0(fn, ".pvar.zst")
        if (!file.exists(pvar.fn))
        {
            pvar.fn <- paste0(fn, ".pvar.gz")
            if (!file.exists(pvar.fn))
                pvar.fn <- paste0(fn, ".pvar")
        }
        # set psam file
        psam.fn <- paste0(fn, ".psam.zst")
        if (!file.exists(psam.fn))
        {
            psam.fn <- paste0(fn, ".psam.gz")
            if (!file.exists(psam.fn))
                psam.fn <- paste0(fn, ".psam")
        }
        # check pgen
        if (!grepl("\\.pgen$", pgen.fn, ignore.case=TRUE))
            pgen.fn <- paste0(fn, ".pgen")
    }
    stopifnot(is.character(pvar.fn), length(pvar.fn)==1L)
    stopifnot(is.character(psam.fn), length(psam.fn)==1L)
    stopifnot(is.character(out.gdsfn), length(out.gdsfn)==1L)
    stopifnot(is.character(compress.geno), length(compress.geno)==1L)
    stopifnot(is.character(compress.annot), length(compress.annot)==1L)
    stopifnot(is.null(reference) | is.character(reference))
    stopifnot(is.null(variant.sel) | is.logical(variant.sel) |
        is.numeric(variant.sel))
    stopifnot(is.null(sample.sel) | is.logical(sample.sel) |
        is.numeric(sample.sel))
    stopifnot(is.numeric(start), length(start)==1L)
    stopifnot(is.numeric(count), length(count)==1L)
    stopifnot(is.character(ignore.chr.prefix))
    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(digest), length(digest)==1L)
    stopifnot(is.logical(balancing), length(balancing)==1L)
    if (is.na(balancing))
        balancing <- isTRUE(getOption("seqarray.balancing", TRUE))
    stopifnot(is.logical(verbose), length(verbose)==1L)
    use.bit1 <- FALSE

    # check files
    if (!file.exists(pgen.fn)) stop("No file ", sQuote(pgen.fn), ".")
    if (!file.exists(pvar.fn)) stop("No file ", sQuote(pvar.fn), ".")
    if (!file.exists(psam.fn)) stop("No file ", sQuote(psam.fn), ".")

    # open pgen file
    if (verbose)
    {
        .cat("##< ", .tm())
        .cat("PLINK2 PGEN to SeqArray GDS:")
        .cat("    pgen file (", .pretty_size(file.size(pgen.fn)), "):")
        .cat("        ", pgen.fn)
        .cat("    pvar file (", .pretty_size(file.size(pvar.fn)), "):")
        .cat("        ", pvar.fn)
        .cat("    genome reference: ", if (length(reference))
            paste(reference, collapse=", ") else "<unset>")
    }
    pgen.fn <- normalizePath(pgen.fn, mustWork=FALSE)
    pvar.fn <- normalizePath(pvar.fn, mustWork=FALSE)
    pvar <- pgen <- NULL
    pvar <- pgenlibr::NewPvar(pvar.fn)
    pgen <- pgenlibr::NewPgen(pgen.fn, pvar=pvar)
    on.exit({
        if (!is.null(pgen)) ClosePgen(pgen)
        if (!is.null(pvar)) ClosePvar(pvar)
    })
    nvar <- nvar_tot <- GetVariantCt(pgen)
    nsamp <- GetRawSampleCt(pgen)

    if (!is.null(variant.sel))
    {
        variant.sel <- .validate_variant_sel(variant.sel, nvar_tot)
        nvar <- length(variant.sel)
    }

    if (is.na(start) || start<1L) start <- 1L
    if (!is.integer(start)) start <- as.integer(start)
    count <- as.integer(count)
    if (is.na(count) || count<1L)
    {
        count <- as.integer(nvar - start + 1L)
    } else {
        count <- min(start + count, nvar + 1L) - start
    }

    # read psam file
    if (verbose)
    {
        .cat("    psam file (", .pretty_size(file.size(psam.fn)), "):")
        .cat("        ", psam.fn)
        .cat("        reading ...")
    }
    fam <- .read_psam(psam.fn)
    if (anyDuplicated(fam$IID) == 0L)
    {
        sample.id <- fam$IID
    } else {
        stop("Sample IDs in PLINK2 PGEN are not unique (IID column)!")
    }
    if (nsamp != nrow(fam))
    {
        stop(sprintf(
            "%d samples are in pgen, but only read %d rows in the psam file.",
            nsamp, nrow(fam)))
    }

    # sample selection
    nsamp_tot <- nsamp
    if (!is.null(sample.sel))
    {
        sample.sel <- .validate_sample_sel(sample.sel, nsamp_tot)
        sample.id <- sample.id[sample.sel]
        nsamp <- length(sample.sel)
        # reopen pgen with sample subset
        ClosePgen(pgen)
        pgen <- pgenlibr::NewPgen(pgen.fn, pvar=pvar,
            sample_subset=sample.sel)
    }

    if (verbose)
    {
        if (nsamp == nsamp_tot)
        {
            .cat("    # of samples: ", nsamp)
        } else {
            .cat("    # of selected samples: ", nsamp, " (from ",
                nsamp_tot, " samples)")
        }
        if (nvar == nvar_tot)
        {
            .cat("    # of variants: ", nvar)
        } else {
            .cat("    # of selected variants: ", nvar, " (from ",
                nvar_tot, " variants)")
        }
        message("    Output:\n        ", out.gdsfn)
        if (start!=1L || count!=nvar)
            .cat("        (starting from ", start, ", count: ", count, ")")
    }

    # the number of parallel tasks
    pnum <- SeqArray:::.NumParallel(parallel)
    ptmpfn <- NULL
    if (pnum > 1L)
    {
        if (count >= pnum)
        {
            if (isTRUE(balancing))
            {
                m <- SeqArray:::.get_bl_multiple()
                pnum <- as.integer(pnum * m)
                if (verbose)
                    message("    (workload balancing)")
            }

            fn <- sub("^([^.]*).*", "\\1", basename(out.gdsfn))
            psplit <- SeqArray:::.file_split(count, pnum, start)
            # need unique temporary file names
            ptmpfn <- character()
            while (length(ptmpfn) < pnum)
            {
                s <- tempfile(pattern=sprintf("%s_tmp%02d_",
                    fn, length(ptmpfn)+1L), tmpdir=dirname(out.gdsfn))
                file.create(s)
                if (!(s %in% ptmpfn)) ptmpfn <- c(ptmpfn, s)
            }
            if (verbose)
            {
                .cat("    output to path: ", file.path(dirname(out.gdsfn), ""))
                message(sprintf("    writing to %d files [%s]:", pnum, .tm()))
                message(paste0(sprintf("        %s [%s..%s]", basename(ptmpfn),
                    .pretty(psplit[[1L]]),
                    .pretty(psplit[[1L]] + psplit[[2L]] - 1L)),
                    collapse="\n"))
            }

            # reset memory
            ClosePgen(pgen); pgen <- NULL
            ClosePvar(pvar); pvar <- NULL
            sel <- NULL
            if (is.numeric(variant.sel))
            {
                sel <- rep(as.raw(0L), nvar_tot)
                sel[variant.sel] <- as.raw(1L)
                sel <- memCompress(sel, "gzip")
                remove(variant.sel)
            }
            gc(FALSE, reset=TRUE, full=TRUE)

            # show information
            update_info <- function(u)
            {
                .cat("        |> ", u[1L], " (", u[2L], " variants) [",
                    .tm(), " done]")
                flush.console()
                NULL
            }
            if (!isTRUE(verbose)) update_info <- "none"

            # conversion in parallel
            seqParallel(parallel, pnum,
                FUN = function(i, pnum, pgen.fn, pvar.fn, psam.fn,
                    compress.geno, compress.annot, ptmpfn, psplit, sel,
                    sample.sel)
                {
                    # set job status
                    .Call(SEQ_SetJobStatus, i, pnum)
                    on.exit(.Call(SEQ_SetJobStatus, 0L, 0L))
                    cnt <- psplit[[2L]][i]
                    if (!is.null(sel))
                        sel <- which(as.logical(memDecompress(sel, "gzip")))
                    tryCatch(
                    {
                        pgen2gds::seqPGEN2GDS(pgen.fn, pvar.fn, psam.fn,
                            ptmpfn[i], compress.geno=compress.geno,
                            compress.annot=compress.annot,
                            variant.sel=sel, sample.sel=sample.sel,
                            start=psplit[[1L]][i], count=cnt,
                            optimize=FALSE, digest=FALSE, parallel=FALSE,
                            verbose=FALSE)
                        # return the process index & # of variants
                        c(i, cnt)  
                    }, error = function(e) {
                        # capture full traceback
                        trace <- capture.output({
                            cat("Error: ", e$message, "\n", sep="")
                            traceback()
                        })
                        con <- file(paste0(ptmpfn[i], ".progress.txt"), "at")
                        writeLines(trace, con)
                        close(con)
                        stop(e$message)
                    })
                }, split = "none", .combine = update_info,
                pnum=pnum, pgen.fn=pgen.fn, pvar.fn=pvar.fn, psam.fn=psam.fn,
                compress.geno=compress.geno,
                compress.annot=compress.annot,
                ptmpfn=ptmpfn, psplit=psplit, sel=sel,
                sample.sel=sample.sel
            )
            if (verbose)
            {
                message("    done splitting (", .tm(), ")")
                message("    --------")
            }

            # reopen the file
            pvar <- pgenlibr::NewPvar(pvar.fn)
            if (!is.null(sel))
            {
                # restore variant.sel
                variant.sel <- which(as.logical(memDecompress(sel, "gzip")))
            }

        } else {
            pnum <- 1L
            message("No use of parallel environment!")
        }
    }

    # create GDS file
    dstfile <- createfn.gds(out.gdsfn)
    on.exit({
        if (!is.null(dstfile)) closefn.gds(dstfile)
    }, add=TRUE)
    put.attr.gdsn(dstfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(dstfile$root, "FileVersion", "v1.0")

    n <- addfolder.gdsn(dstfile, "description")
    put.attr.gdsn(n, "source.format", "PLINK2 PGEN Format")
    if (length(reference))
        add.gdsn(n, "reference", reference, visible=FALSE)

    # add sample.id
    if (verbose) message("    sample.id  ", appendLF=FALSE)
    n <- add.gdsn(dstfile, "sample.id", sample.id, compress=compress.annot,
        closezip=TRUE)
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)

    # add variant.id
    if (verbose) message("    variant.id  ", appendLF=FALSE)
    n <- add.gdsn(dstfile, "variant.id", seq.int(start, length.out=count),
        compress=compress.annot, closezip=TRUE)
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)

    # write variant annotations (chromosome, position, allele, id)
    .write_variant_annot(dstfile, pvar, variant.sel, start, count,
        compress.annot, ignore.chr.prefix, digest, verbose)

    # add a folder for genotypes & phase
    ngen <- addfolder.gdsn(dstfile, "genotype")
    put.attr.gdsn(ngen, "VariableName", "GT")
    put.attr.gdsn(ngen, "Description", "Genotype")
    npha <- addfolder.gdsn(dstfile, "phase")

    # add annotation folder (already created by .write_variant_annot)
    nann <- index.gdsn(dstfile, "annotation")

    # add nodes for genotypes
    if (verbose)
        .cat("    genotype [", .tm(), "] ...")
    n_g <- add.gdsn(ngen, "data", storage=ifelse(use.bit1, "bit1", "bit2"),
        valdim=c(2L, nsamp, 0L), compress=compress.geno)
    n_i <- add.gdsn(ngen, "@data", storage="uint8", compress=compress.annot,
        visible=FALSE)
    n_p <- add.gdsn(npha, "data", storage="bit1", valdim=c(nsamp, 0L),
        compress=compress.geno)

    # progress information
    progfilename <- paste0(out.gdsfn, ".progress.txt")
    progfile <- file(progfilename, "wt")
    progfile_to_rm <- FALSE
    writeLines(paste("#", .tm()), progfile)
    writeLines(paste("Input:", basename(pgen.fn)), progfile)
    if (start!=1L || count!=nvar)
        writeLines(paste0("    start: ", start, ", count: ", count), progfile)
    writeLines(paste("Output:", basename(out.gdsfn)), progfile)
    flush(progfile)
    on.exit({
        close(progfile)
        if (progfile_to_rm) unlink(progfilename, force=TRUE)
    }, add=TRUE)

    if (pnum <= 1L)
    {
        # process genotypes
        buf <- IntBuf(pgen)  # an integer vector
        ii <- integer(1L)    # variant index
        ia <- integer(1L)    # allele index
        read_gt_fc <- quote(ReadHardcalls(pgen, buf, ii, ia))
        allele_num_fc <- quote(GetAlleleCt(pvar, ii))
        # call C
        .Call(SEQ_PGEN_Geno_Import, read_gt_fc, allele_num_fc, buf, ii, ia,
            variant.sel, new.env(), dstfile$root, start, count, progfile,
            verbose)
        # remove temporary variables
        remove(read_gt_fc, allele_num_fc, buf, ii, ia)
        # close the nodes
        readmode.gdsn(n_g)
        if (verbose) message("      ", appendLF=FALSE)
        SeqArray:::.DigestCode(n_g, digest, verbose)
        readmode.gdsn(n_i)
        SeqArray:::.DigestCode(n_i, digest, FALSE)
        readmode.gdsn(n_p)
        SeqArray:::.DigestCode(n_p, digest, FALSE)
        sync.gds(dstfile)

    } else {
        varnm <- c("genotype/data", "genotype/@data", "phase/data")
        # "phase/data"
        # open all temporary files
        for (fn in ptmpfn)
        {
            if (verbose)
                message("        adding ", sQuote(basename(fn)), appendLF=FALSE)
            writeLines(paste("Adding", basename(fn)), progfile)
            flush(progfile)
            # open the gds file
            tmpgds <- openfn.gds(fn)
            n <- objdesp.gdsn(index.gdsn(tmpgds, "variant.id"))$dim
            # merge variables
            for (nm in varnm)
                append.gdsn(index.gdsn(dstfile, nm), index.gdsn(tmpgds, nm))
            # close the file
            closefn.gds(tmpgds)
            if (verbose) .cat(" [", .tm(), " done]")
            writeLines(paste0("    [", .tm(), ", ",
                prettyNum(n, big.mark=",", scientific=FALSE),
                " variants added]"), progfile)
            flush(progfile)
        }
        sync.gds(dstfile)

        # remove temporary files
        writeLines(paste("Done. #", .tm()), progfile)
    }

    # additional nodes for genotype
    n <- add.gdsn(ngen, "extra.index", storage="int32", valdim=c(3L,0L),
        compress=compress.geno, closezip=TRUE)
    put.attr.gdsn(n, "R.colnames",
        c("sample.index", "variant.index", "length"))
    add.gdsn(ngen, "extra", storage="int16", compress=compress.geno,
        closezip=TRUE)
    # additional nodes for phase data
    n <- add.gdsn(npha, "extra.index", storage="int32", valdim=c(3L,0L),
        compress=compress.geno, closezip=TRUE)
    put.attr.gdsn(n, "R.colnames",
        c("sample.index", "variant.index", "length"))
    add.gdsn(npha, "extra", storage="bit1", compress=compress.geno,
        closezip=TRUE)

    # add annotation/qual
    n <- add.gdsn(nann, "qual", storage="float", compress=compress.annot)
    SeqArray:::.append_rep_gds(n, NaN, nvar)
    readmode.gdsn(n)
    if (verbose) message("    annotation/qual", appendLF=FALSE)
    SeqArray:::.DigestCode(n, digest, verbose)

    # add filter
    n <- add.gdsn(nann, "filter", storage="int32", compress=compress.annot)
    SeqArray:::.append_rep_gds(n, as.raw(1L), nvar)
    readmode.gdsn(n)
    put.attr.gdsn(n, "R.class", "factor")
    put.attr.gdsn(n, "R.levels", c("PASS"))
    put.attr.gdsn(n, "Description", c("All filters passed"))
    if (verbose) message("    annotation/filter", appendLF=FALSE)
    SeqArray:::.DigestCode(n, digest, verbose)

    # add the INFO field
    addfolder.gdsn(nann, "info")
    # add the FORMAT field
    addfolder.gdsn(nann, "format")

    # close the file
    closefn.gds(dstfile)
    dstfile <- NULL
    progfile_to_rm <- TRUE
    # delete temporary files
    if (is.character(ptmpfn)) unlink(ptmpfn, force=TRUE)

    # optimize access efficiency
    if (verbose)
        if (optimize) .cat("Done.  # ", .tm()) else .cat("Done.")
    if (optimize)
    {
        if (verbose)
            .cat("Optimize the access efficiency ...")
        cleanup.gds(out.gdsfn, verbose=verbose)
    }
    if (verbose) .cat("##> ", .tm())

    # output
    invisible(normalizePath(out.gdsfn))
}

