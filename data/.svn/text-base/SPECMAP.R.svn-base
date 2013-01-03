#### SPECMAP.R : process SPECMAP data, originally in
#### SPECMAP_018_reversed.R.  NOTE: these data have already been
#### normalised.

local({
  fpath <- file.path('..', 'orig', 'SPECMAP.txt')
  SPECMAP <<- read.table(fpath, header = TRUE)
})

SPECMAP$CE <- -1e3 * SPECMAP$Age
SPECMAP <- SPECMAP[order(SPECMAP$CE), c('CE', 'delta18O')]
