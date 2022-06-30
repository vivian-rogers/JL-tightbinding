push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, "./periodicweyl/")
push!(LOAD_PATH, "./src/")


module runs
using InBi
using SLutils
using Driver
using Constants
using UsefulFunctions
using LinearAlgebra
#main(params)

nx = 20; ny = 1; nz = 1; 
# superlattice basis vectors, in basis of a_1, a_2, a_3
SL1 = [nx; 0; 0]; SL2 = [0; ny; 0]; SL3 = [0; 0; nz]


#runtype = "domainwall"
#runtype = "bulk"
#runtype = "nanopillars"
#runtype = "eggcarton"
#runtype = "neelwall"
#runtype = ""
#runtype = "fmthinfilm"
runtype = "blochwall"
fieldtype = "β"

p = genSL(params, nx, ny, nz, SL1, SL3, SL3, runtype, fieldtype) # generate SL params



if(fieldtype=="A")
	A = Agen(p,runtype,10^8*4*μₑ)
else
	A = βgen(p,runtype,0.5*eV,30.0)
end

main(p,A)


end
