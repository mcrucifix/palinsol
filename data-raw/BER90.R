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

#### load Berger 1990 table from INSOL file

local({
fpath <- file.path("..","inst","extdata", "BER90.IN.gz")
fpath_table2 <- file.path("..","inst","extdata", "GenerateTable2.R")

if (! file.exists(fpath)) # mmm. maybe the packages is already installed
{
fpath <- system.file("extdata", "BER90.IN.gz", package="palinsol")
fpath_table2 <- system.file("extdata", "GenerateTable2.R", package="palinsol")
}



Table4 <- utils::read.table(gzfile(fpath), skip=1, nrow=80)
Table1 <- utils::read.table(gzfile(fpath), skip=161, nrow=1000)
Table5 <- utils::read.table(gzfile(fpath), skip=6481, nrow=1000)

source(fpath_table2)
Table2 <- .GlobalEnv$GenerateTable2(Table1, Table4, Table5, sol="BER90")

# add period
Table4 <- cbind(Table4, data.frame(V5=360*360 / Table4$V3) )

colnames(Table1) <- c('Term','Amp','Rate','Phase','Period')
colnames(Table4) <- c('Term','Amp','Rate','Phase','Period')
colnames(Table5) <- c('Term','Amp','Rate','Phase','Period')

Table1 <<- Table1
Table2 <<- Table2
Table4 <<- Table4
Table5 <<- Table5

})

Table1 <- .GlobalEnv$Table1
Table2 <- .GlobalEnv$Table2
Table4 <- .GlobalEnv$Table4
Table5 <- .GlobalEnv$Table5
#
BER90 <- list(Table1=Table1, Table2=Table2, Table4=Table4, Table5=Table5, psibar =50.41726176, estar = 23.33340950, zeta = 1.60075265)

usethis::use_data(BER90)

#
