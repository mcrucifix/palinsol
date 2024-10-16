load_all () 
# install_github("mcrucifix/palinsol@loic_1")


require(palinsol)

times=seq(-10e3,5)*1e3

# generate frequencies

# Frequency_models <- sapply(seq(50,62,0.25), function (i) {
#     dummy <- compute_tables(prea=i)
#     freqs <- with(dummy, c(k = psibar, 
#                p1= Table2[1,5], 
#                p2= Table2[3,5], 
#                p3= Table1[1,5]))})
# 

# data from AstroGeo22 for -300 M
# http://www.astrogeo.eu/wp-content/uploads/2022/11/AstroGeo22.txta
 #  T(Ga)     a (RE)      amin        amax         lod       lodmin      lodmax       obli       oblimin     oblimax       pre           premin        premax

 # 0.370   57.673582   57.335846   57.885557   21.818198   21.547225   21.991389   22.231292   22.227687   22.236342     61.062980     60.584835     61.826758

 estar <- 22.231292 
 prma <- 61.062980    
 ala <- prma * cos(estar * pi / 180)
 

 # this is still code in development. bea=0 to tell palinsol not to
 # take into account initial conditions
DevonTables  <- compute_tables(bea=0, ala=ala, ha = estar, prma = prma)

 # you can force estar to be the value supplied

DevonTables$estar <- estar

# print the main precession frequencies

DevonTables$Table2[seq(10),]

# print the main obliquity frequencies


DevonTables$Table1[seq(10),]

# compare with the standard BER90

BER90$Table2[seq(10),]

# compute a solution

 o1 <- function(t)  table_based_solution(t, DevonTables)
 time_solution <- sapply ( times, o1)

dev.off()



