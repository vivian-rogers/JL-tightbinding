push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")

module runs
using InBi
using SLutils
using Driver
using Constants
using UsefulFunctions
using LinearAlgebra
#main(params)

nx = 30; ny = 30; nz = 10; 
# superlattice basis vectors, in basis of a_1, a_2, a_3
SL1 = [nx; 0; 0]; SL2 = [0; ny; 0]; SL3 = [0; 0; nz]


#runtype = "domainwall"
#runtype = "bulk"
runtype = "nanopillars"
#runtype = "eggcarton"
#runtype = "afmthinfilm"

p = genSL(params, nx, ny, nz, SL1, SL3, SL3, runtype) # generate SL params




A = Agen(p,runtype,10^7*4*μₑ)


main(p,A)


end
