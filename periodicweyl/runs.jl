#push!(Base.load_path(), "./")
#push!(Base.load_path(), "../src/")
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

nx = 30; ny = 1; nz = 60; 
# superlattice basis vectors, in basis of a_1, a_2, a_3
SL1 = [nx; 0; 0]; SL2 = [0; ny; 0]; SL3 = [0; 0; nz]


#runtype = "domainwall"
#runtype = "bulk"
#runtype = "nanopillars"
#runtype = "eggcarton"
runtype = "blochwall"
fieldtype = "β"

p = genSL(params, nx, ny, nz, SL1, SL3, SL3, runtype, fieldtype) # generate SL params



if(fieldtype=="A")
	A = Agen(p,runtype,10^8*4*μₑ)
else
	A = βgen(p,runtype,0.2*eV)
end

main(p,A)


end
