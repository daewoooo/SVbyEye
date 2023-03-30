message("============================")
message("Check lifting ranges to PAF alignments")

## Define range(s) to lift
roi.gr <- as("chr17:46645907-46697277", "GRanges")
## Get PAF alignments to lift to
paf.file1 <- system.file("extdata", "test_lift1.paf", package = "SVbyEye")
paf.file2 <- system.file("extdata", "test_lift2.paf", package = "SVbyEye")
## Load PAF alignments
paf.table1 <- readPaf(paf.file = paf.file1)
paf.table2 <- readPaf(paf.file = paf.file2)
## Lift target range to query coordinates
lift1 <- liftRangesToAlignment(paf.table = paf.table1, gr = roi.gr, direction = "target2query")
lift2 <- liftRangesToAlignment(paf.table = paf.table2, gr = roi.gr, direction = "target2query")
## Make checks
## Expected mapping position based on minimap: 1015588-1066973
expect_that(start(lift1), equals(1015588))
expect_that(end(lift1), equals(1066972))
## Expected mapping position based on minimap: 557941-609322
expect_that(start(lift2), equals(557942))
expect_that(end(lift2), equals(609322))
