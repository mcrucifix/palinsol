# This code is not ready for use. This is a experimentation with the automatic
# to an attribution routine available in the sister package "gtseries". 


epi <- read.table('../extdata/La88_IOEPL.dat', col.names=c('id','Freq','Amp','Phases','name'), fill=TRUE, skip=2, nrow=80)
iom <- read.table('../extdata/La88_IOEPL.dat', col.names=c('id','Freq','Amp','Phases','name'), fill=TRUE, skip=84, nrow=80)
full = rbind(epi, iom)
withnames <- which(full$name != "")
gref <- full$Freq[withnames]
names(gref) <- full$name[withnames]
require(gtseries)
attributeTones(infreq=full$Freq, omegas=gref, frac=FALSE, tol2=1.e-2, tol1=1.e-5, keepPositives=FALSE)
# but at this point the tone attribution might be a bit too generous
# file is here in arcsec/year
full$names <- attributeTones(infreq=abs(full$Freq), omegas=abs(gref), frac=FALSE, tol2=1.e-3, tol1=1.e-4)
epi <- full[seq(80),]
iom <- full[80+seq(80),]
epi
class(epi) <- "discreteSpectrum"
class(iom) <- "discreteSpectrum"


# note: the 'sign' of the attribution is not always correct. Check code
# the attribution is reasonably convincing on frequencies but not on phases. But the latter are supper fragile.
# check what Laskar 88 did ? 
