### Read the edc-co2-2008.txt  data containing the  CO2
### measures at EPICA DOME-C by the European Drilling Project
### see this file for proper references

### The original data files are in the subdirectory 'orig', to
### preserve their meta-data (R does not like such files in data)

local({
  fpath <- file.path('..', 'orig', 'EPICA', 'edc-co2-2008.txt')
  edcCO2 <<- read.table(fpath, skip = 769, header = TRUE, nrows = 1095)
})

edcCO2$CE <- -edcCO2$Age.yrBP.
edcCO2$CO2 <- edcCO2$CO2.ppmv.
edcCO2 <- edcCO2[order(edcCO2$CE), c('CE', 'CO2')]
