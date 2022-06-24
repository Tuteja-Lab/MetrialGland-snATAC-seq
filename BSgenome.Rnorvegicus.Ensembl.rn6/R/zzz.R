###
###

.pkgname <- "BSgenome.Rnorvegicus.Ensembl.rn6"

.seqnames <- c(1:20, "X", "Y")

.circ_seqs <- NULL

.mseqnames <- NULL

.onLoad <- function(libname, pkgname)
{
    if (pkgname != .pkgname)
        stop("package name (", pkgname, ") is not ",
             "the expected name (", .pkgname, ")")
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make and export BSgenome object.
    bsgenome <- BSgenome(
        organism="Rattus norvegicus",
        common_name="Rat",
        provider="Ensembl",
        provider_version="rn6",
        release_date="Sept. 2020",
        release_name="rnor6",
        source_url="ftp://ftp.ensembl.org/pub/release-98/fasta/rattus_norvegicus/",
        seqnames=.seqnames,
        circ_seqs=.circ_seqs,
        mseqnames=.mseqnames,
        seqs_pkgname=pkgname,
        seqs_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)

    old_objname <- "Rnorvegicus"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

