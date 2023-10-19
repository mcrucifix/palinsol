require(palinsol)

generate_partial_solution <- function(nobl=2, npre=2){
  BER90_Truncated <- palinsol::BER90
  BER90_Truncated$Table2 <- BER90_Truncated$Table2[seq(npre),]
  BER90_Truncated$Table4 <- BER90_Truncated$Table4[seq(npre),]
  BER90_Truncated$Table5 <- BER90_Truncated$Table5[seq(npre),]
  BER90_Truncated$Table1 <- BER90_Truncated$Table1[seq(nobl),]
  return(BER90_Truncated)
}

partial <- function(t, partial_solution){
   return(table_based_solution(t, tab_solution = partial_solution, degree = FALSE, compute_esinw_from_table2 = TRUE))
}

# example

times <- seq(500)*2000
npre=5; nobl=2;
partial_solution <- generate_partial_solution(npre=npre, nobl=nobl)
my_solution <- function(t) partial(t, partial_solution)


O = sapply(times, my_solution)

pdf('test_partial_solutions.pdf')
plot(times, O['eps',]*180/pi, main=sprintf("obliquity with %d components", nobl), type='l' ) ; 
plot(times, O['esinw',]*180/pi, main=sprintf("prcession with %d components", npre), type='l' ) ; 
dev.off()

