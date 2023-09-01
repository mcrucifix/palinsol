times <- seq(-1000,1000,2)*1e3

OT <- compute_tables()
OT2 <- compute_tables(prea=58)

# OT <- palinsol::BER90

fn <- function(t) table_based_solution(t, OT,  compute_esinw_from_table2=FALSE)['esinw']
fn2 <- function(t) table_based_solution(t, OT, compute_esinw_from_table2=TRUE)['esinw']

sol1 <- sapply(times, fn)
sol2 <- sapply(times, fn2)

plot(sol1, type='l')
lines(sol2,col='red')
lines(sol1-sol2,col='blue')

s1= OT$Table2[seq(10000),3]
s2= OT2$Table2[seq(10000),3]

plot(s1-s2)

amp <- OT$Table2[seq(10000), 2]





