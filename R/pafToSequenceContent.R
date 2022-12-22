# ## Get PAF to plot ##
# paf.file <- system.file("extdata", "test_getFASTA.paf", package="SVbyEye")
# ## Read in PAF
# paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = 'cg')
# ## Bin PAF
# paf.table <- pafToBins(paf.table = paf.table, binsize = 10000)
#
# asm.fasta <- system.file("extdata", "test_getFASTA_query.fasta", package="SVbyEye")
#
# paf2sequenceContent <- function(paf.table=NULL, sequence.content='gc') {
#
#   query.gr <- GenomicRanges::GRanges(seqnames = unique(paf.table$q.name),
#                                       ranges = IRanges::IRanges(start = paf.table$q.start,
#                                                               end = paf.table$q.end)
#   )
#   regions2FASTA(gr = query.gr, asm.fasta = asm.fasta)
#
# }
