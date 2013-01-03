### Read the edc-ch4-2008.txt  data containing the  CH4
### measures at EPICA DOME-C by the European Drilling Project
### see this file for proper references

### The original data files are in the subdirectory 'orig', to
### preserve their meta-data (R does not like such files in data)
local ({
fpath <- file.path('..', 'orig', 'EPICA', 'edc-ch4-2008.txt')
edcCH4 <<- read.table(fpath, skip = 153, header = TRUE, nrows = 2152)

edcCH4$CE <<- -edcCH4$Gas_Age
edcCH4 <<- edcCH4[order(edcCH4$CE), ]})
