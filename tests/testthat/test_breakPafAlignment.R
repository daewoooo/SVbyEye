message("============================")
message("Check PAF breaking at indels")

## Test sample with no strand state changes
paf.file <- system.file("extdata", "test2.paf", package = "SVbyEye")
## Read in PAF alignment
paf.aln <- readPaf(paf.file = paf.file)
## Break PAF alignment at indels of 1 kbp and longer
paf.broken <- breakPafAlignment(paf.aln = paf.aln, min.deletion.size = 1000, min.insertion.size = 1000)
paf.aln <- paf.broken$M
paf.indels <- paf.broken$SVs
## Make checks
expect_that(nrow(paf.aln), equals(2))
expect_that(paf.aln$q.start, equals(c(343761, 289054)))
expect_that(paf.aln$q.end, equals(c(350771, 343760)))
expect_that(paf.aln$t.start, equals(c(1012155, 1052905)))
expect_that(paf.aln$t.end, equals(c(1019110, 1107668)))
expect_that(paf.aln$t.end, equals(c(1019110, 1107668)))
expect_that(paf.indels$q.start, equals(343761))
expect_that(paf.indels$t.start, equals(1019111))
expect_that(paf.indels$cg, equals("33794D"))
