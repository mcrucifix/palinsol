#### load La_ioepl.data
fpath <- file.path("inst","extdata", "La88_IOEPL.dat")

epi <- utils::read.table(fpath, skip=2, nrow=80, fill=TRUE)
io <- utils::read.table(fpath, skip=84, nrow=80, fill=TRUE)

sec2pi <- pi/3600/180

epi <- data.frame( 
                   Freq  = epi$V2*sec2pi,
                   Amp   = epi$V3/1e8, 
                   Phases = epi$V4*pi/180)

i2o <- data.frame( 
                   Freq  = io$V2*sec2pi,
                   Amp   = io$V3/1e8, 
                   Phases = io$V4*pi/180)

class(epi) <- "discreteSpectrum"
attr(epi,"nfreq") <- 80

class(io) <- "discreteSpectrum"
attr(io,"nfreq") <- 80

La88 <- list(epi=epi, i2o=i2o)

usethis::use_data(La88, overwrite=TRUE)

