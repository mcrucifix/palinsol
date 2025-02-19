#### load La_ioepl.data
fpath <- file.path("inst","extdata", "bre78.dat")

epi <- utils::read.table(fpath, skip=0, nrow=19, fill=TRUE)
io <- utils::read.table(fpath, skip=19, nrow=15, fill=TRUE)

sec2pi <- pi/3600/180

epi <- data.frame( 
                   Freq  = epi$V3*sec2pi,
                   Amp   = epi$V2, 
                   Phases = epi$V4*pi/180)

io <- data.frame( 
                   Freq  = io$V3*sec2pi,
                   Amp   = io$V2, 
                   Phases = io$V4*pi/180)

class(epi) <- "discreteSpectrum"
attr(epi,"nfreq") <- 19

class(io) <- "discreteSpectrum"
attr(io,"nfreq") <- 15

Bre73 <- list(epi=epi, io=io)

usethis::use_data(Bre73, overwrite=TRUE)

