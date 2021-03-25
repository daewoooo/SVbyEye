## Load required libraries
library(SVbyEye)
library(ggplot2)
library(scales)
library(gggenes)

#paf.file <- '/home/porubsky/SVbyEye/chr12_17521907_18108016_misorient.paf'
paf.file <- '/home/porubsky/SVbyEye_devel/CHM13_putativeNovelINV_5MbExtraFlanks.paf'

#coords.old <- paf2coords(paf.file = paf.file, chromosome = 'chr7')
coords <- paf2coords(paf.file = paf.file, min.mapq = 20, min.align.len = 10000)
coords <- coords[coords$align.id == 'chr1:137260403-147334366__chr1',]

## Get y-axis labels
q.range <- range(coords$seq.pos[coords$seq.id == 'query'])
t.range <- range(coords$seq.pos[coords$seq.id == 'target'])
q.labels <- pretty(q.range)
t.labels <- pretty(t.range)
q.breaks <- SVbyEye::q2t(x = q.labels, q.range = q.range, t.range = t.range)
t.breaks <- t.labels

## Get x-axis labels
seq.labels <- c(unique(coords$seq.name[coords$seq.id == 'query']), 
                unique(coords$seq.name[coords$seq.id == 'target']))

n <- 100
strength <- 0.5
coords.plt <- coords[order(coords$group, coords$y), ]
coords1 <- coords.plt[c(TRUE, FALSE, TRUE, FALSE), ]
coords2 <- coords.plt[c(FALSE, TRUE, FALSE, TRUE), ]
coords1 <- add.control.points(coords1, strength)
coords1 <- ggforce::StatBezier$compute_panel(coords1, scales, n)
coords2 <- add.control.points(coords2[rev(seq_len(nrow(coords2))), ], strength)
coords2 <- ggforce::StatBezier$compute_panel(coords2, scales, n)
coords.plt <- rbind(coords1, coords2)

plt <- ggplot(data = coords.plt) + 
  geom_polygon(aes(x=x, y=y, group=group), color=alpha('gray', 0.25), fill=alpha('gray', 0.5)) +
  scale_y_continuous(breaks = c(1, 2), labels = seq.labels) +
  scale_x_continuous(breaks = q.breaks, labels = comma(q.labels),
                     sec.axis = sec_axis(trans = y ~ ., breaks = t.breaks, labels = comma(t.labels)))

## Add arrows
start <- coords$x[c(T, T, F, F)]
end <- coords$x[c(F, F, T, T)]
y <- coords$y[c(T, T, F, F)]
#direction <- coords$direction[c(T, T, F, F)]
group <- coords$group[c(T, T, F, F)]
plt.df <- data.frame(start=start, end=end, y=y, group=group)
plt.df$direction <- ifelse(plt.df$start < plt.df$end, '+', '-')

plt <- plt + geom_gene_arrow(data=plt.df, aes(xmin = start, xmax = end, y = y, color= direction, fill = direction), arrowhead_height = unit(3, 'mm')) +
  theme_bw()

## Add SD arrowheads
SDs.df <- read.table("/media/porubsky/DPBACKUP/HGSVC_WorkingBackup//GRCh38_annot/segDup_GRCh38.gz")
SDs.df <- SDs.df[,c(2,3,4,27)]
colnames(SDs.df) <- c('seqnames', 'start', 'end', 'fracMatch')
SDs.gr <- makeGRangesFromDataFrame(SDs.df, keep.extra.columns = TRUE)
## Keep only SDs >=98% identical
#SDs.gr <- SDs.gr[SDs.gr$fracMatch >= 0.98]

## Convert genomic ranges to polygon coords
#sub.gr <- GRanges(seqnames='chr7', ranges=IRanges(start=53256586, end=62876688))

SDs.gr <- SDs.gr[start(SDs.gr) > min(plt.df$start) & end(SDs.gr) < max(plt.df$end) & seqnames(SDs.gr) == 'chr1']
strand(SDs.gr) <- sample(c('+', '-'), length(SDs.gr), 1)
plt.df <- as.data.frame(SDs.gr)

library(ggnewscale)
plt + new_scale_fill() + new_scale_color() +
  geom_arrowhead(data=plt.df, aes(xmin=start, xmax=end, y=2.05, color=fracMatch, fill=fracMatch))
  #geom_polygon(data=plt.df, aes(x=x, y=3, group=group, color=fracMatch, fill=fracMatch))



library(SVbyEye)
library(ggplot2)
paf.file <- '/home/porubsky/SVbyEye_devel/chr12_17521907_18108016_misorient.paf'

coords <- paf2coords(paf.file = paf.file, min.mapq = 20, min.align.n = 1, min.align.len = 10000)

## Get y-axis labels
q.range <- range(coords$seq.pos[coords$seq.id == 'query'])
t.range <- range(coords$seq.pos[coords$seq.id == 'target'])
q.labels <- pretty(q.range)
t.labels <- pretty(t.range)
q.breaks <- SVbyEye::q2t(x = q.labels, q.range = q.range, t.range = t.range)
t.breaks <- t.labels

## Get x-axis labels
seq.labels <- c(unique(coords$seq.name[coords$seq.id == 'query']), 
                unique(coords$seq.name[coords$seq.id == 'target']))

ggplot2::ggplot(coords) +
  geom_miropeats(aes(x, y, group = group), fill='gray', alpha=0.5, color='black') +
  scale_y_continuous(breaks = c(1, 2), labels = seq.labels) +
  scale_x_continuous(breaks = q.breaks, labels = comma(q.labels),
                     sec.axis = sec_axis(trans = y ~ ., breaks = t.breaks, labels = comma(t.labels)))


########################################################################################################
# n <- 100
# strength <- 0.5
# 
# infile <- '/home/porubsky/SVbyEye/CHM13_putativeNovelINV_5MbExtraFlanks.paf'
# paf <- read.table(infile, stringsAsFactors = FALSE)
# ## Keep only first 12 columns
# paf <- paf[,c(1:12)]
# ## Add header
# header <- c('q.name', 'q.len', 'q.start', 'q.end', 'strand', 't.name', 't.len', 't.start', 't.end', 'n.match', 'aln.len', 'mapq') 
# colnames(paf) <- header
# ## Flip start-end if strand == '-'
# paf[paf$strand == '-', c('t.start','t.end')] <- rev(paf[paf$strand == '-', c('t.start','t.end')])
# 
# ## Keep only chromosome 7
# #paf <- paf[grep(paf$q.name, pattern = 'chr7'),]
# paf <- paf[grep(paf$t.name, pattern = 'chr7'),]
# ## Keep min alignment length of 1kb
# paf <- paf[paf$aln.len >= 1000,]
# 
# ## Sync scales
# # dist <- max(diff(sort(c(paf$q.start, paf$q.end, paf$t.start, paf$t.end))))
# # paf$q.start <- paf$q.start + dist
# # paf$q.end <- paf$q.end + dist
# q.range <- range(c(paf$q.start, paf$q.end))
# t.range <- range(c(paf$t.start, paf$t.end))
# paf$q.start.trans <- q2t(x = paf$q.start, q.range = q.range, t.range = t.range)
# paf$q.end.trans <- q2t(x = paf$q.end, q.range = q.range, t.range = t.range)
# 
# ## Get y-axis labels
# q.labels <- pretty(q.range)
# q.breaks <- q2t(x = q.labels, q.range = range(q.labels), t.range = t.range)
# t.labels <- pretty(t.range)
# t.breaks <- t.labels
# 
# ## Get x-axis labels
# seq.labels <- c(unique(paf$q.name), unique(paf$t.name))
# 
# ## Vectorize data transformation ##
# #y1 <- c(rbind(paf$q.start, paf$t.start))
# #y2 <- c(rbind(paf$q.end, paf$t.end))
# y1 <- c(rbind(paf$q.start.trans, paf$t.start))
# y2 <- c(rbind(paf$q.end.trans, paf$t.end))
# 
# x1 <- rep(c(1,2), times=nrow(paf))
# x2 <- rep(c(1,2), times=nrow(paf))
# 
# group <- rep(1:nrow(paf), each=2)
# 
# coords1 <- data.frame(x=x1, y=y1, group=group)
# coords1 <- coords1[order(coords1$group, coords1$x, coords1$y), ]
# coords2 <- data.frame(x=x2, y=y2, group=group)
# coords2 <- coords2[order(coords2$group, coords2$x, coords2$y), ]
# 
# ## Add fill scale
# coords1$fill.scale <- rep(paf$n.match / paf$aln.len, each=2)
# coords2$fill.scale <- rep(paf$n.match / paf$aln.len, each=2)
# #coords1$x <- rep(c(1,2), times=nrow(paf))
# #coords2$x <- rep(c(2,1), times=nrow(paf))
# 
# coords1 <- add_controls(coords1, strength)
# coords1 <- StatBezier$compute_panel(coords1, scales, n)
# #coords2 <- add_controls(coords2, strength)
# coords2 <- add_controls(coords2[rev(seq_len(nrow(coords2))), ], strength)
# coords2 <- StatBezier$compute_panel(coords2, scales, n)
# 
# 
# coords <- rbind(coords1, coords2)
# ggplot(data = coords) + 
#   #geom_polygon(aes(x=x, y=y, group=group, fill=as.factor(group))) +
#   geom_polygon(aes(x=x, y=y, group=group, fill=fill.scale)) +
#   #geom_polygon(aes(x=x, y=y, group=group), alpha=0.5) +
#   #scale_fill_manual(values = brewer.pal(n = 3, name = 'Set1')) +
#   scale_fill_gradient(low = 'white', high='red') +
#   scale_x_continuous(breaks = c(1, 2), labels = seq.labels) +
#   scale_y_continuous(breaks = q.breaks, labels = comma(q.labels),
#                      sec.axis = sec_axis(trans = y ~ ., breaks = t.breaks, labels = comma(t.labels))) +
#   coord_flip()
# 
# 
# ## Test arrowhead
# example_genes$strand <- ifelse(example_genes$direction == 1, '+', '-')
# ggplot2::ggplot(example_genes, ggplot2::aes(xmin = start, xmax = end, y = molecule, fill = gene, strand=strand)) +
#   geom_arrowhead() + ggplot2::facet_wrap(~ molecule, scales = "free")
# 
# 
# test.df <- SDs.df[sample(1:nrow(SDs.df), size = 100),]
# test.df <- test.df[test.df$seqnames %in% c('chr1', 'chr2', 'chr3'),]
# test.df$strand <- '+'
# 
# ggplot2::ggplot(test.df, ggplot2::aes(xmin = start, xmax = end, y = seqnames, strand=strand, fill = fracMatch)) +
#   geom_arrowhead() + ggplot2::facet_wrap(~ seqnames, scales = "free")
