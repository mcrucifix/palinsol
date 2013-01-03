#### load Berger 1990 table from INSOL file

fpath <- file.path('..', 'orig', 'BER90.IN')

Table4_90 <- read.table(fpath, skip=1, nrow=80)
Table1_90 <- read.table(fpath, skip=161, nrow=1000)
Table5_90 <- read.table(fpath, skip=6481, nrow=1000)

rm(list = c('fpath'))
