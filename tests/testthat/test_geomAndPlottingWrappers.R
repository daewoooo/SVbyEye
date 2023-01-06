message("=================================================")
message("Check geom_miropeat and plotMiro wrapper function")

## Get PAF to plot ##
paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
## Read in PAF
paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
## Make a plot ##
plt <- plotMiro(paf.table = paf.table, color.by = "direction")
## Make checks
expect_true(is(plt, "ggplot"))
expect_is(plt$layers[[1]]$geom, "GeomPolygon")
expect_is(plt$layers[[1]]$stat, "StatMiropeats")

message("=================================================")
message("Check geom_wide_arc and plotSelf wrapper function")

## Get PAF to plot
paf.file <- system.file("extdata", "test2.paf", package = "SVbyEye")
## Read in PAF
paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
## Make a plot
## Plot alignment as horizontal dotplots and color by alignment directionality
plt1 <- plotSelf(paf.table = paf.table, shape = "segment")
## Plot alignment as arcs and color by alignment directionality
plt2 <- plotSelf(paf.table = paf.table, shape = "arc")
## Make checks
expect_true(is(plt1, "ggplot"))
expect_is(plt1$layers[[1]]$geom, "GeomSegment")
expect_true(is(plt2, "ggplot"))
expect_is(plt2$layers[[1]]$geom, "GeomPolygon")
expect_is(plt2$layers[[1]]$stat, "StatWideArc")

message("==============================")
message("Check plotAVA wrapper function")

## Get PAF to plot
paf.file <- system.file("extdata", "test_ava.paf", package = "SVbyEye")
## Read in PAF
paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
## Make a plot
## Color by alignment directionality
plt <- plotAVA(paf.table = paf.table, color.by = "direction")
## Make checks
expect_true(is(plt, "ggplot"))
expect_is(plt$layers[[1]]$geom, "GeomPolygon")
expect_is(plt$layers[[1]]$stat, "StatMiropeats")
