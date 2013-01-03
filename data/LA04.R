#### load La04 table

fpathn <- file.path('..', 'orig', 'INSOLN.LA2004.BTL.ASC')
fpathp <- file.path('..', 'orig', 'INSOLP.LA2004.BTL.ASC')

la04past     <<- read.table(fpathn, col.names=c('time','ecc','eps','varpi'))
la04future   <<- read.table(fpathp, col.names=c('time','ecc','eps','varpi'))


la04past['varpi'] <- (la04past['varpi'] - pi ) %% (2*pi)
la04future['varpi'] <- (la04future['varpi'] - pi ) %% (2*pi)


rm(list = c('fpathn','fpathp'))
