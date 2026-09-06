# test-gds2pgen.R — tests for seqGDS2PGEN()

pgen_fn <- system.file("extdata", "plink2_gen.pgen", package = "pgen2gds")
pvar_fn <- system.file("extdata", "plink2_gen.pvar", package = "pgen2gds")
psam_fn <- system.file("extdata", "plink2_gen.psam", package = "pgen2gds")

# parse the header and index of a PGEN file (storage mode 0x10)
.pgen_index <- function(fn)
{
    r <- readBin(fn, "raw", file.size(fn))
    M <- readBin(r[4:7], "integer", 1L, 4L)
    N <- readBin(r[8:11], "integer", 1L, 4L)
    ctrl <- as.integer(r[12L])
    lb <- bitwAnd(ctrl, 3L) + 1L
    eight <- bitwAnd(ctrl, 4L) != 0L
    nblk <- ceiling(M / 65536)
    off <- 12L + 8L * nblk
    fpos <- readBin(r[13:20], "integer", 1L, 8L)
    if (eight)
    {
        vt <- as.integer(r[(off+1L):(off+M)])
        off <- off + M
    } else {
        nvt <- ceiling(M / 2)
        v <- as.integer(r[(off+1L):(off+nvt)])
        vt <- as.vector(rbind(bitwAnd(v, 15L), bitwShiftR(v, 4L)))[seq_len(M)]
        off <- off + nvt
    }
    lens <- vapply(seq_len(M), function(i)
        sum(as.integer(r[(off+(i-1L)*lb+1L):(off+i*lb)]) * 256^(0:(lb-1L))), 0)
    recs <- vector("list", M)
    p <- fpos
    for (i in seq_len(M))
    {
        recs[[i]] <- if (lens[i] > 0) r[(p+1L):(p+lens[i])] else raw(0)
        p <- p + lens[i]
    }
    list(magic=r[1:3], M=M, N=N, ctrl=ctrl, fpos=fpos, vrtype=vt, lens=lens,
        records=recs, size=file.size(fn))
}

# unordered genotype pairs with NA coded as -1
.sorted_geno <- function(g)
{
    g[is.na(g)] <- -1L
    apply(g, 2:3, identity) -> g
    for (k in seq_len(dim(g)[3L])) g[,,k] <- apply(g[,,k], 2L, sort)
    g
}

# a small VCF with multiallelic, phased, partially phased and rare variants
.write_test_vcf <- function(fn)
{
    writeLines(c(
        "##fileformat=VCFv4.2",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##contig=<ID=1>",
        paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", sprintf("S%02d", 1:10)), collapse="\t"),
        paste(c("1", "100", "v1_unphased", "A", "G", ".", "PASS", ".", "GT",
            "0/0", "0/1", "1/1", "./.", "0/0", "0/1", "0/0", "1/1", "0/0", "0/1"),
            collapse="\t"),
        paste(c("1", "200", "v2_allphased", "A", "G", ".", "PASS", ".", "GT",
            "0|0", "0|1", "1|0", "1|1", "0|0", "1|0", "0|1", "0|0", "1|1", "0|1"),
            collapse="\t"),
        paste(c("1", "300", "v3_partphased", "A", "G", ".", "PASS", ".", "GT",
            "0/0", "0|1", "1/0", "1|1", "0/1", "1|0", "0/1", "./.", "0|0", "1|0"),
            collapse="\t"),
        paste(c("1", "400", "v4_tri", "A", "G,T", ".", "PASS", ".", "GT",
            "0/0", "0/1", "0/2", "1/1", "1/2", "2/2", "0|2", "2|0", "1|2", "2|1"),
            collapse="\t"),
        paste(c("1", "500", "v5_quad", "A", "G,T,C", ".", "PASS", ".", "GT",
            "0/0", "0/3", "3/3", "1/3", "2|3", "3|1", "0/1", "0|2", "2|0", "./."),
            collapse="\t"),
        paste(c("1", "600", "v6_hexa", "A", "G,T,C,GG,TT", ".", "PASS", ".", "GT",
            "0/5", "5/5", "4|5", "5|4", "1/5", "0/0", "0/4", "3/4", "0/1", "1/1"),
            collapse="\t"),
        paste(c("1", "700", "v7_allref", "A", "G", ".", "PASS", ".", "GT",
            rep("0/0", 10L)), collapse="\t"),
        paste(c("1", "800", "v8_rare", "A", "G", ".", "PASS", ".", "GT",
            "0/0", "0/0", "0/0", "0/1", rep("0/0", 6L)), collapse="\t")
    ), fn)
    invisible(fn)
}

.rm_prefix <- function(prefix)
    unlink(paste0(prefix, c(".pgen", ".pvar", ".psam", ".log")), force=TRUE)


# ---- round trip: pgen -> gds -> pgen ----------------------------------------

test_that("seqGDS2PGEN reproduces the original PGEN records exactly",
{
    gds_fn <- tempfile(fileext=".gds")
    out <- tempfile()
    on.exit({
        unlink(c(gds_fn, paste0(gds_fn, ".progress.txt")), force=TRUE)
        .rm_prefix(out)
    }, add=TRUE)

    seqPGEN2GDS(pgen_fn, out.gdsfn=gds_fn, optimize=FALSE, verbose=FALSE)
    fns <- seqGDS2PGEN(gds_fn, out, verbose=FALSE)
    expect_length(fns, 3L)
    expect_true(all(file.exists(fns)))
    expect_equal(basename(fns), paste0(basename(out), c(".pgen", ".pvar", ".psam")))

    a <- .pgen_index(pgen_fn)
    b <- .pgen_index(paste0(out, ".pgen"))
    expect_identical(b$magic, as.raw(c(0x6c, 0x1b, 0x10)))
    expect_equal(b$M, a$M)
    expect_equal(b$N, a$N)
    # the same per-variant encodings and byte-identical variant records
    expect_identical(b$vrtype, a$vrtype)
    expect_identical(b$lens, a$lens)
    expect_identical(b$records, a$records)

    # pvar
    expect_identical(seqReadPVAR(paste0(out, ".pvar")), seqReadPVAR(pvar_fn))
    # psam header and IDs
    psam <- read.table(paste0(out, ".psam"), header=TRUE, comment.char="",
        sep="\t", stringsAsFactors=FALSE)
    expect_true(all(c("X.IID", "PAT", "MAT", "SEX") %in% names(psam)))
    orig <- read.table(psam_fn, header=TRUE, comment.char="", sep="\t",
        stringsAsFactors=FALSE)
    expect_identical(psam$X.IID, orig$X.IID)
    expect_identical(psam$SEX, orig$SEX)
    expect_identical(psam$PAT, orig$PAT)
})

test_that("pgenlibr reads identical hardcalls from the converted PGEN",
{
    gds_fn <- tempfile(fileext=".gds")
    out <- tempfile()
    on.exit({
        unlink(c(gds_fn, paste0(gds_fn, ".progress.txt")), force=TRUE)
        .rm_prefix(out)
    }, add=TRUE)

    seqPGEN2GDS(pgen_fn, out.gdsfn=gds_fn, optimize=FALSE, verbose=FALSE)
    seqGDS2PGEN(gds_fn, out, verbose=FALSE)

    pv1 <- pgenlibr::NewPvar(pvar_fn)
    pg1 <- pgenlibr::NewPgen(pgen_fn, pvar=pv1)
    pv2 <- pgenlibr::NewPvar(paste0(out, ".pvar"))
    pg2 <- pgenlibr::NewPgen(paste0(out, ".pgen"), pvar=pv2)
    on.exit({
        pgenlibr::ClosePgen(pg1); pgenlibr::ClosePgen(pg2)
        pgenlibr::ClosePvar(pv1); pgenlibr::ClosePvar(pv2)
    }, add=TRUE)
    expect_equal(pgenlibr::GetRawSampleCt(pg2), pgenlibr::GetRawSampleCt(pg1))
    expect_equal(pgenlibr::GetVariantCt(pg2), pgenlibr::GetVariantCt(pg1))
    b1 <- pgenlibr::IntBuf(pg1); b2 <- pgenlibr::IntBuf(pg2)
    ndiff <- 0L
    for (i in seq_len(pgenlibr::GetVariantCt(pg1)))
    {
        pgenlibr::ReadHardcalls(pg1, b1, i)
        pgenlibr::ReadHardcalls(pg2, b2, i)
        if (!identical(b1, b2)) ndiff <- ndiff + 1L
    }
    expect_equal(ndiff, 0L)
})


# ---- multiallelic and phased genotypes ------------------------------------

test_that("seqGDS2PGEN preserves multiallelic genotypes and phase",
{
    vcf_fn <- tempfile(fileext=".vcf")
    gds_fn <- tempfile(fileext=".gds")
    gds_back <- tempfile(fileext=".gds")
    out <- tempfile()
    on.exit({
        unlink(c(vcf_fn, gds_fn, gds_back, paste0(gds_back, ".progress.txt")),
            force=TRUE)
        .rm_prefix(out)
    }, add=TRUE)

    .write_test_vcf(vcf_fn)
    SeqArray::seqVCF2GDS(vcf_fn, gds_fn, verbose=FALSE)
    seqGDS2PGEN(gds_fn, out, verbose=FALSE)

    ix <- .pgen_index(paste0(out, ".pgen"))
    expect_equal(ix$M, 8L)
    expect_equal(ix$N, 10L)
    # 8-bit vrtypes since phase is present (bit 2 of the control byte)
    expect_true(bitwAnd(ix$ctrl, 4L) != 0L)
    # multiallelic patch (0x08) and phase (0x10) flags
    expect_true(all(bitwAnd(ix$vrtype[4:6], 8L) != 0L))
    expect_true(all(bitwAnd(ix$vrtype[c(2L, 3L, 4L, 5L, 6L)], 16L) != 0L))
    expect_true(all(bitwAnd(ix$vrtype[c(1L, 7L, 8L)], 16L) == 0L))

    # the pvar keeps all alleles
    pv <- seqReadPVAR(paste0(out, ".pvar"))
    expect_identical(pv$allele[4:6], c("A,G,T", "A,G,T,C", "A,G,T,C,GG,TT"))

    # round trip through seqPGEN2GDS: the same unordered genotypes
    seqPGEN2GDS(paste0(out, ".pgen"), out.gdsfn=gds_back, optimize=FALSE,
        verbose=FALSE)
    f1 <- SeqArray::seqOpen(gds_fn)
    f2 <- SeqArray::seqOpen(gds_back)
    on.exit({ SeqArray::seqClose(f1); SeqArray::seqClose(f2) }, add=TRUE)
    G1 <- SeqArray::seqGetData(f1, "genotype")
    G2 <- SeqArray::seqGetData(f2, "genotype")
    expect_identical(SeqArray::seqGetData(f2, "allele"),
        SeqArray::seqGetData(f1, "allele"))
    expect_identical(.sorted_geno(G2), .sorted_geno(G1))

    # phase and haplotype order of the biallelic variants via pgenlibr
    P1 <- SeqArray::seqGetData(f1, "phase")
    pvar <- pgenlibr::NewPvar(paste0(out, ".pvar"))
    pgen <- pgenlibr::NewPgen(paste0(out, ".pgen"), pvar=pvar)
    on.exit({ pgenlibr::ClosePgen(pgen); pgenlibr::ClosePvar(pvar) }, add=TRUE)
    expect_true(pgenlibr::HardcallPhasePresent(pgen))
    ac <- pgenlibr::IntAlleleCodeBuf(pgen)
    pp <- pgenlibr::BoolBuf(pgen)
    for (i in c(1L, 2L, 3L, 7L, 8L))
    {
        pgenlibr::ReadAlleles(pgen, ac, i, pp)
        a <- matrix(ac, nrow=2L)
        g <- G1[,,i]
        het <- !is.na(g[1L,]) & !is.na(g[2L,]) & (g[1L,] != g[2L,])
        # phased heterozygotes keep the haplotype order
        phased <- het & (P1[,i] == 1L)
        expect_true(all(as.logical(pp)[phased]))
        expect_false(any(as.logical(pp)[het & !phased]))
        expect_identical(a[, phased, drop=FALSE], g[, phased, drop=FALSE])
        # the other genotypes are the same unordered pairs
        s1 <- apply(replace(a, is.na(a), -1L), 2L, sort)
        s2 <- apply(replace(g, is.na(g), -1L), 2L, sort)
        expect_identical(s1, s2)
    }
})

test_that("use.phase=FALSE drops the phase track",
{
    vcf_fn <- tempfile(fileext=".vcf")
    gds_fn <- tempfile(fileext=".gds")
    out <- tempfile()
    on.exit({ unlink(c(vcf_fn, gds_fn), force=TRUE); .rm_prefix(out) }, add=TRUE)
    .write_test_vcf(vcf_fn)
    SeqArray::seqVCF2GDS(vcf_fn, gds_fn, verbose=FALSE)
    seqGDS2PGEN(gds_fn, out, use.phase=FALSE, verbose=FALSE)
    ix <- .pgen_index(paste0(out, ".pgen"))
    expect_true(bitwAnd(ix$ctrl, 4L) == 0L)
    expect_true(all(bitwAnd(ix$vrtype, 16L) == 0L))
    pvar <- pgenlibr::NewPvar(paste0(out, ".pvar"))
    pgen <- pgenlibr::NewPgen(paste0(out, ".pgen"), pvar=pvar)
    on.exit({ pgenlibr::ClosePgen(pgen); pgenlibr::ClosePvar(pvar) }, add=TRUE)
    expect_false(pgenlibr::HardcallPhasePresent(pgen))
})


# ---- half calls -------------------------------------------------------------

test_that("half.call controls haploid and half-missing genotypes",
{
    vcf_fn <- tempfile(fileext=".vcf")
    gds_fn <- tempfile(fileext=".gds")
    out1 <- tempfile(); out2 <- tempfile()
    on.exit({
        unlink(c(vcf_fn, gds_fn), force=TRUE)
        .rm_prefix(out1); .rm_prefix(out2)
    }, add=TRUE)
    writeLines(c(
        "##fileformat=VCFv4.2",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##contig=<ID=1>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5",
        "1\t100\tv1\tA\tG\t.\tPASS\t.\tGT\t0\t1\t0/1\t.\t1/."
    ), vcf_fn)
    SeqArray::seqVCF2GDS(vcf_fn, gds_fn, verbose=FALSE)

    hardcalls <- function(prefix)
    {
        pv <- pgenlibr::NewPvar(paste0(prefix, ".pvar"))
        pg <- pgenlibr::NewPgen(paste0(prefix, ".pgen"), pvar=pv)
        on.exit({ pgenlibr::ClosePgen(pg); pgenlibr::ClosePvar(pv) })
        b <- pgenlibr::IntBuf(pg)
        pgenlibr::ReadHardcalls(pg, b, 1L)
        as.integer(b)
    }
    seqGDS2PGEN(gds_fn, out1, half.call="missing", verbose=FALSE)
    expect_identical(hardcalls(out1), c(NA, NA, 1L, NA, NA))
    seqGDS2PGEN(gds_fn, out2, half.call="haploid", verbose=FALSE)
    expect_identical(hardcalls(out2), c(0L, 2L, 1L, NA, 2L))
})


# ---- filters, IDs, chromosome prefix ---------------------------------------

test_that("seqGDS2PGEN respects seqSetFilter()",
{
    gds_fn <- tempfile(fileext=".gds")
    out <- tempfile()
    on.exit({
        unlink(c(gds_fn, paste0(gds_fn, ".progress.txt")), force=TRUE)
        .rm_prefix(out)
    }, add=TRUE)
    seqPGEN2GDS(pgen_fn, out.gdsfn=gds_fn, optimize=FALSE, verbose=FALSE)
    f <- SeqArray::seqOpen(gds_fn)
    on.exit(SeqArray::seqClose(f), add=TRUE)
    samp <- SeqArray::seqGetData(f, "sample.id")
    SeqArray::seqSetFilter(f, sample.id=samp[c(2L, 5L, 10L, 100L)],
        variant.id=11:30, verbose=FALSE)
    seqGDS2PGEN(f, out, verbose=FALSE)

    ix <- .pgen_index(paste0(out, ".pgen"))
    expect_equal(ix$N, 4L)
    expect_equal(ix$M, 20L)
    pv <- seqReadPVAR(paste0(out, ".pvar"))
    expect_identical(pv$pos, SeqArray::seqGetData(f, "position"))
    psam <- read.table(paste0(out, ".psam"), header=TRUE, comment.char="",
        sep="\t", stringsAsFactors=FALSE)
    expect_identical(psam$X.IID, samp[c(2L, 5L, 10L, 100L)])

    # genotypes of the subset
    pvar <- pgenlibr::NewPvar(paste0(out, ".pvar"))
    pgen <- pgenlibr::NewPgen(paste0(out, ".pgen"), pvar=pvar)
    on.exit({ pgenlibr::ClosePgen(pgen); pgenlibr::ClosePvar(pvar) }, add=TRUE)
    G <- SeqArray::seqGetData(f, "genotype")
    b <- pgenlibr::IntBuf(pgen)
    for (i in 1:20)
    {
        pgenlibr::ReadHardcalls(pgen, b, i)
        expect_identical(as.integer(b), as.integer(colSums(G[,,i] != 0L)))
    }
})

test_that("write.rsid and chr_prefix control the pvar columns",
{
    gds_fn <- tempfile(fileext=".gds")
    out1 <- tempfile(); out2 <- tempfile()
    on.exit({
        unlink(c(gds_fn, paste0(gds_fn, ".progress.txt")), force=TRUE)
        .rm_prefix(out1); .rm_prefix(out2)
    }, add=TRUE)
    seqPGEN2GDS(pgen_fn, out.gdsfn=gds_fn, variant.sel=1:5, optimize=FALSE,
        verbose=FALSE)
    f <- SeqArray::seqOpen(gds_fn)
    on.exit(SeqArray::seqClose(f), add=TRUE)

    seqGDS2PGEN(f, out1, write.rsid="chr_pos_ref_alt", chr_prefix="chr",
        verbose=FALSE)
    pv <- read.table(paste0(out1, ".pvar"), header=TRUE, comment.char="",
        sep="\t", colClasses="character")
    expect_identical(pv$X.CHROM,
        paste0("chr", SeqArray::seqGetData(f, "chromosome")))
    expect_identical(pv$ID, gsub(":|,", "_",
        SeqArray::seqGetData(f, "$chrom_pos_allele")))

    seqGDS2PGEN(f, out2, write.rsid="annot_id", verbose=FALSE)
    pv <- read.table(paste0(out2, ".pvar"), header=TRUE, comment.char="",
        sep="\t", colClasses="character")
    expect_identical(pv$ID, SeqArray::seqGetData(f, "annotation/id"))
    expect_identical(pv$X.CHROM, SeqArray::seqGetData(f, "chromosome"))
})


# ---- plink2 validation (optional) -------------------------------------------

test_that("plink2 --validate accepts the converted PGEN",
{
    skip_if(Sys.which("plink2") == "", "plink2 is not available")
    vcf_fn <- tempfile(fileext=".vcf")
    gds_fn <- tempfile(fileext=".gds")
    out <- tempfile()
    on.exit({ unlink(c(vcf_fn, gds_fn), force=TRUE); .rm_prefix(out) }, add=TRUE)
    .write_test_vcf(vcf_fn)
    SeqArray::seqVCF2GDS(vcf_fn, gds_fn, verbose=FALSE)
    seqGDS2PGEN(gds_fn, out, verbose=FALSE)
    rv <- system2("plink2", c("--pfile", out, "--validate", "--out", out),
        stdout=TRUE, stderr=TRUE)
    expect_equal(attr(rv, "status"), NULL)
    expect_true(any(grepl("done", rv, fixed=TRUE)))
})


# ---- error handling ---------------------------------------------------------

test_that("seqGDS2PGEN rejects invalid arguments",
{
    expect_error(seqGDS2PGEN(pgen_fn, 42), "is.character")
    expect_error(seqGDS2PGEN(pgen_fn, "x", half.call="foo"), "should be one of")
    expect_error(seqGDS2PGEN(pgen_fn, "x", write.rsid="foo"), "should be one of")
    expect_error(seqGDS2PGEN(pgen_fn, "x", use.phase=NA), "is.na")
})


# ---- import side: sample.annotation ----------------------------------------

test_that("seqPGEN2GDS stores psam columns in sample.annotation",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE), add=TRUE)
    seqPGEN2GDS(pgen_fn, out.gdsfn=out, sample.sel=1:20, optimize=FALSE,
        verbose=FALSE)
    f <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(f), add=TRUE)
    expect_true(gdsfmt::exist.gdsn(f, "sample.annotation/sex"))
    expect_true(gdsfmt::exist.gdsn(f, "sample.annotation/father"))
    expect_true(gdsfmt::exist.gdsn(f, "sample.annotation/mother"))
    sex <- SeqArray::seqGetData(f, "sample.annotation/sex")
    expect_length(sex, 20L)
    expect_true(all(sex %in% c("M", "F", "")))
    psam <- read.table(psam_fn, header=TRUE, comment.char="", sep="\t",
        stringsAsFactors=FALSE)
    expect_identical(sex, c("", "M", "F")[psam$SEX[1:20] + 1L])
})

test_that("seqReadPVAR accepts a pgenlibr pvar object",
{
    pv <- pgenlibr::NewPvar(pvar_fn)
    on.exit(pgenlibr::ClosePvar(pv))
    df <- seqReadPVAR(pv, sel=1:3)
    expect_equal(nrow(df), 3L)
    expect_identical(df, seqReadPVAR(pvar_fn, sel=1:3))
})
