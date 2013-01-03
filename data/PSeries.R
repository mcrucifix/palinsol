#### load Berger 1978 table from INSOL file

fpath <- file.path('..', 'orig', 'INSOL.IN')

Table4 <- read.table(fpath, skip=6, nrow=19)
Table1 <- read.table(fpath, skip=25, nrow=47)
Table5 <- read.table(fpath, skip=72, nrow=78)

rm(list = c('fpath'))

local(
{
  P<- 50.439273  
  zeta <- 3.392506
  ## PsiBar

  g <- Table4$V3
  M <- Table4$V2
  beta <- Table4$V4
  F <- Table5$V2/60./60.*pi/180
  f <- Table5$V3
  delta <- Table5$V4

  ## division in 3 groups, as in Table 13 and Table 14 of Berger & Loutre, Ac. Roy. 1990

  Fre <- c(g+P,outer(g,f,"+")+P,outer(g,f,"-")+P)
  Amp <- c(M,outer(M,F,"*")/2.,outer(M,F,"*")/2. )
  Pha <- c(beta+zeta,outer(beta,delta,"+")+zeta,outer(beta,delta,"-")+zeta)


  ## regroup similar frequencies

  tol <- 0.0001
  Ntrun <- 200.

  Order <- order(abs(Amp),decreasing=TRUE)

  Amp <- Amp[Order]
  Fre <- Fre[Order]
  ## truncates the first 200 terms
  Fre <- Fre[1:Ntrun]
  Amp <- Amp[1:Ntrun]


  N <- length(Fre)
  for (i in 1:(N-1)) {
   for (j in (i+1):N) {
   if (abs(Fre[j] - Fre[i]) < tol) {
    Amp[i] <- Amp[i]+Amp[j]
    Amp[j] <- 0
    } } }

  Order <- order(abs(Amp),decreasing=TRUE)

  Amp <- Amp[Order]
  Fre <- Fre[Order]
  Pha <- (Pha[Order]+180)%%360.
  Per <- 360*60*60/Fre/1000.
  Table2 <<- data.frame(Index=seq(1,length(Amp)),Amp=Amp,Fre=Fre,Pha=Pha,Per=abs(Per[Order])) 

})

Table2 <- Table2
