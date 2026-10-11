#!/usr/bin/env Rscript
# Reproducible, wholly synthetic gallery example. No biological measurements.
args <- commandArgs(FALSE)
self <- sub('^--file=', '', grep('^--file=', args, value = TRUE)[1])
root <- dirname(dirname(normalizePath(self)))
dest <- file.path(root, 'examples', 'generated', 'hierarchy')
input <- file.path(dest, 'input')
dir.create(input, recursive = TRUE, showWarnings = FALSE)
set.seed(59341)

# Irregular lobes illustrate a dense embedding; these are NOT computed UMAPs.
centers <- rbind(
  c(412,534),c(418,444),c(454,393),c(464,274),c(391,364),
  c(295,161),c(363,120),c(404,95),c(261,491),c(338,530),
  c(298,406),c(306,280),c(372,287),c(338,369),c(219,399),
  c(150,351),c(224,310),c(164,278),c(180,203),c(123,250),
  c(266,330),c(523,421),c(513,501),c(565,531),c(611,418),
  c(586,357),c(548,470),c(575,445))
sx <- c(21,29,30,32,32,28,25,18,32,33,28,18,31,22,23,28,33,32,22,25,24,25,23,16,16,25,20,18)
sy <- c(24,30,34,34,38,52,32,23,38,27,38,30,24,28,24,31,28,28,38,30,30,24,25,16,22,20,21,20)
angle <- c(-.5,-.5,.3,.1,.2,.6,.5,.8,.3,-.2,-.2,-.2,-.1,-.4,.1,-.7,-.8,-.8,-.8,-.2,-.4,-.7,.7,0,.7,.4,.7,.3)
weights <- c(1.6,1.2,1.3,2.4,1,1.3,1.1,.4,1.7,1.3,2.2,1.2,.7,.6,1,1.5,1.7,1.6,1.1,.8,.5,.8,.9,.5,.5,.6,.5,.35)
n <- round(80000 * weights / sum(weights))
n[28] <- 80000 - sum(n[1:27])
cells <- do.call(rbind, lapply(seq_len(28), function(i) {
  a <- rnorm(n[i]); b <- rnorm(n[i])
  tail <- runif(n[i]) < .1
  a[tail] <- a[tail] * 1.25
  b <- b + .1 * (a^2 - 1)
  radius <- sqrt(a*a + b*b)
  taper <- ifelse(radius > 0, tanh(radius/2.3)/(radius/2.3), 1)
  xx <- sx[i] * a * taper; yy <- sy[i] * b * taper
  data.frame(cell_id = paste0('sim_', i, '_', seq_len(n[i])),
    cluster_id = as.character(i),
    umap_1 = (centers[i,1] + xx*cos(angle[i]) - yy*sin(angle[i]))/40,
    umap_2 = -(centers[i,2] + xx*sin(angle[i]) + yy*cos(angle[i]))/40)
}))
cells <- cells[sample(nrow(cells)), ]

# Five abstract parent populations with related colors within each family.
colors <- read.csv(file.path(root, 'assets/style/reference_28_cluster_colors.csv'))$color
group <- rep(LETTERS[1:5], c(5,3,6,7,7))
within <- sequence(c(5,3,6,7,7))
clusters <- data.frame(cluster_id = as.character(1:28),
  label = sprintf('%02d  %s.%d', 1:28, group, within), color = colors,
  group_id = group, group_label = paste('Group', group))
write_tsv <- function(x, file) write.table(x, file.path(input, file), sep='\t', quote=FALSE, row.names=FALSE)
write_tsv(cells, 'cells.tsv')
write_tsv(clusters, 'clusters.tsv')
write_tsv(data.frame(cluster_id=as.character(1:28), x=centers[,1]/40, y=-centers[,2]/40), 'label_positions.tsv')

# Controlled synthetic blocks make the encodings visible and reproducible.
parent <- match(group, LETTERS[1:5])
marker <- expand.grid(cluster_id=as.character(1:28), feature=paste0('M',1:6), stringsAsFactors=FALSE)
i <- as.integer(marker$cluster_id); j <- match(marker$feature, paste0('M',1:6))
on <- j == parent[i] | (j == 6 & i %% 4 == 0)
marker$color_value <- ifelse(on, runif(nrow(marker),1,2), runif(nrow(marker),-1.2,.25))
marker$size_value <- ifelse(on, runif(nrow(marker),65,98), runif(nrow(marker),4,45))
write_tsv(marker, 'markers.tsv')

programs <- expand.grid(cluster_id=as.character(1:28), feature=paste0('P',1:6), stringsAsFactors=FALSE)
i <- as.integer(programs$cluster_id); j <- match(programs$feature, paste0('P',1:6))
base <- rbind(c(1.35,-.4,-.2,.3,-.6,.1), c(-.6,1.65,.4,-.3,-.1,-.4),
              c(.15,-.5,1.3,.6,-.35,.15), c(-.7,-.8,.1,1.55,.45,-.25),
              c(-.35,.1,-.7,.05,1.25,.95))
programs$value <- base[cbind(parent[i], j)] + .35*sin(within[i]*1.2+j) + rnorm(nrow(programs),0,.15)
write_tsv(programs, 'programs.tsv')

signals <- expand.grid(cluster_id=as.character(1:28), feature=paste0('S',1:8), stringsAsFactors=FALSE)
i <- as.integer(signals$cluster_id); j <- match(signals$feature,paste0('S',1:8))
base <- rbind(c(.3,-.25,1.2,.15,-.4,-.65,.1,.35), c(-.4,1.5,.1,-.5,.35,-.1,.55,-.6),
              c(.1,-.55,.4,1.1,.75,-.3,-.6,.2), c(-.2,.2,-.65,.1,1.6,.55,-.45,-.2),
              c(-.65,.4,-.2,-.6,.25,-.4,1.25,1.05))
signals$value <- base[cbind(parent[i], j)] + .4*sin(within[i]*.9+j*.6) + rnorm(nrow(signals),0,.18)
write_tsv(signals, 'signals.tsv')

config <- list(
  title='Synthetic atlas · 80,000 cells · 28 subtypes',
  prefix='hierarchy_atlas', input_dir='input', output_dir='output',
  layouts='portrait', mini_mode='group', mini_placement='sidebar',
  main_mm=160, particle_factor=.65, preview_dpi=180,
  row_height_mm=4.7, group_gap_mm=4.5, sidebar_min_mm=28, legend_space_mm=32,
  modules=list(
    list(type='bubble', file='markers.tsv', title='Marker features', features=paste0('M',1:6),
      limits=c(-2,2), ticks=c('−2','0','2'), legend='Simulated z-score', ramp='pathway',
      size_max=100, size_breaks=c(25,50,100), size_legend='Simulated detection (%)',
      semantics='Artificial feature z-scores and detection percentages; no expression matrix analysed.'),
    list(type='heatmap', file='programs.tsv', title='Programs', features=paste0('P',1:6),
      limits=c(-2,2), ticks=c('−2','0','2'), legend='Simulated z-score', ramp='pathway',
      semantics='Artificial program patterns solely for a visual example.'),
    list(type='heatmap', file='signals.tsv', title='Signatures', features=paste0('S',1:8),
      limits=c(-2,2), ticks=c('−2','0','2'), legend='Simulated z-score', ramp='pathway',
      semantics='Artificial signature patterns; no inference of biological activity.')))
writeLines(c('ATLAS <-', capture.output(dput(config))), file.path(dest,'config.R'))
writeLines(c('SYNTHETIC DATA — rendering demonstration only.',
  'Seed: 59341; 80,000 generated points; 28 abstract subtypes; 5 parent groups.',
  'Coordinates are not a computed UMAP. All features, values, and memberships are artificial.'),
  file.path(dest,'DATA_NOTICE.txt'))
cat('Wrote synthetic hierarchy:', nrow(cells), 'points,', nrow(clusters), 'subtypes.\n')
