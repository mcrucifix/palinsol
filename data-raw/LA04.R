# Copyright (c) 2012 Michel Crucifix <michel.crucifix@uclouvain.be>

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject

# the following conditions:

# The above copyright notice and this permission notice shall be
# incluudedin all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND INFRINGEMENT
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR

# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# When using this package for actual applications, always
# cite the authors of the original insolation solutions 
# Berger, Loutre and/or Laskar, see details in man pages


# ------------------------------------------------------------------
# R Code developed for R version 2.15.2 (2012-10-26) -- "Trick or Treat"
# ------------------------------------------------------------------ 

#### load La04 table
local({
fpathn <- file.path("..","inst","extdata", "INSOLN.LA2004.BTL.ASC.gz")
fpathp <- file.path("..","inst","extdata", "INSOLP.LA2004.BTL.ASC.gz")

if (! file.exists(fpathn)) # mmm. maybe the packages is already installed
{
fpathn <- system.file("extdata", "INSOLN.LA2004.BTL.ASC.gz", package="palinsol")
fpathp <- system.file("extdata", "INSOLP.LA2004.BTL.ASC.gz", package="palinsol")
}



la04past     <- utils::read.table(gzfile(fpathn), col.names=c('time','ecc','eps','varpi'))
la04future   <- utils::read.table(gzfile(fpathp), col.names=c('time','ecc','eps','varpi'))


la04past['varpi'] <- (la04past['varpi'] - pi ) %% (2*pi)
la04future['varpi'] <- (la04future['varpi'] - pi ) %% (2*pi)

#LA04 <<- list(la04past=la04past, la04future = la04future)

la04past <<- la04past
la04future <<- la04future
})

la04past <- .GlobalEnv$la04past
la04future <- .GlobalEnv$la04future

LA04 <- list(la04past=la04past, la04future=la04future)

usethis::use_data(LA04, overwrite=True)
