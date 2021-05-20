message("==========================")
message("Check breakpoint detection")

## Test sample with no strand state changes
file <- list.files(pattern='nobreaks_')
fragments <- get(load(file))
breakpoints <- runBreakpointr(bamfile = fragments, ID=basename(file), chromosomes='chr8', windowsize=1000, binMethod='reads', conf=0.99)
breaks <- breakpoints$breaks
counts <- breakpoints$counts
expect_that(length(breaks), equals(0))
expect_that(counts$states, equals('wc'))

## Test sample with strand state changes
file <- list.files(pattern='withbreaks_')
fragments <- get(load(file))
breakpoints <- runBreakpointr(bamfile = fragments, ID=basename(file), chromosomes='chr8', windowsize=1000, binMethod='reads', conf=0.99)
breaks <- breakpoints$breaks
counts <- breakpoints$counts
expect_that(length(breaks), equals(2))
expect_that(breaks$deltaW, equals(c(981,964)))
expect_that(counts$states, equals(c('cc','ww','cc')))
