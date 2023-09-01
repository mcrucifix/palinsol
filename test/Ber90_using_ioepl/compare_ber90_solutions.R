times=seq(-1e3,1e3,2)*2e3


ber90_values_from_paper <- compute_tables (bea = 23.44579, prea = 50.273147, prega = -2514.27, ala = 54.9066, apoa = 17.3919 )
ber90_values_from_optim <- compute_tables (bea=23.44579,prea=50.27251,prega= -2464.25259,ala=54.88959,apoa=18.12131)

b1 <- sapply(times, function(t) table_based_solution(t, palinsol::BER90))   # reproduces BER90 published solution exactly
b2 <- sapply(times, function(t) table_based_solution(t, ber90_values_from_paper))   # reproduces BER90 based on published parameters
b3 <- sapply(times, function(t) table_based_solution(t, ber90_values_from_optim))   # reproduces BER90 based on optmising parameters to get the same genreal precession rate and phase as publised

PRECESS <- lapply(list(b1,b2,b3), function(sol) (sol['ecc',] * sin (sol['varpi',])))

# plot ( times, PRECESS[[1]], type='l', xlim=c(-1e6, -0.9e6)) 
# lines ( times, PRECESS[[2]], type='l', col='blue') 
# lines ( times, PRECESS[[3]], type='l', col='red') 


plot(times, PRECESS[[1]]-PRECESS[[3]])

