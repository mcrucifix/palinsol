times=seq(-10e3,5)*1e3
pregas <- seq(50,62,0.25)

# generate frequencies

# Frequency_models <- sapply(seq(50,62,0.25), function (i) {
#     dummy <- compute_tables(prea=i)
#     freqs <- with(dummy, c(k = psibar, 
#                p1= Table2[1,5], 
#                p2= Table2[3,5], 
#                p3= Table1[1,5]))})
# 
DevonTables  <- compute_tables(prea=60.250000)
ModernTables <- compute_tables(prea=50.273147)


plot ( ModernTables$Table2[,3] - DevonTables$Table2[,3])

 o1 <- function(t)  table_based_solution(t, DevonTables)
 o2 <- function(t)  table_based_solution(t, ModernTables)

 DevonExample <- sapply(times,o1)
 ModernExample <- sapply(times,o2)


ModernTables <- ModernTables[seq(12700),]
DevonTables <- DevonTables[seq(12700),]

O1 <- order(ModernTables$Table2[,5])
O2 <- order(DevonTables$Table2[,5])

# ModernSpectrum <- ModernTables$Table2[O1,c(5,2)]
# DevonSpecturm <- DevonTables$Table2[O2,c(5,2)]

pdf('modern_vs_devonian_1.pdf')

plot(times, ModernExample['esinw',], type='l');
lines(times, ModernExample['ecc',], type='l', lty=2);
lines(times, DevonExample['esinw',], type='l', col='red');
lines(times, DevonExample['ecc',], type='l', col='red');

legend('topright',c('modern','Devon'), col=c('black','red'))

dev.off()

pdf('modern_vs_devonian_2.pdf')

plot(times, ModernExample['esinw',], type='l', xlim=c(-500000,0));
lines(times, ModernExample['ecc',], type='l', lty=2);
lines(times, DevonExample['esinw',], type='l', col='red');
lines(times, DevonExample['ecc',], type='l', col='red');

legend('topright',c('modern','Devon'), col=c('black','red'))

dev.off()



