push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")

module runs
using InBi
using Driver
using Constants
using UsefulFunctions
using LinearAlgebra
#main(params)

nx = 1; ny = 1; nz = 128; 
# superlattice basis vectors, in basis of a_1, a_2, a_3
SL1 = [nx; 0; 0]; SL2 = [0; ny; 0]; SL3 = [0; 0; nz]


#runtype = "domainwall"
runtype = "nanopillars"

p = genSL(params, nx, ny, nz, SL1, SL3, SL3, runtype) # generate SL params



function Agen(runtype="",M₀=0) # modify vector potential as desired
	if(runtype=="nanopillars")
		function npA(R::Vector{Float64})
			Aval = zeros(3)
			n = 4
			for i = -4:4
				for j = -4:4
					R₀ = i*p.SLa₁ .+ j*p.SLa₂ .+ p.SLa₃ .+ 5*p.a₃
					Aval += μ₀/(4*π) * (M₀*[0;0;1])×(R-R₀)/((R-R₀)⋅(R-R₀))
				end
			end
			return [0;0;0]
		end
		return npA
	elseif(runtype=="domainwall")
		# egg-carton type potential
		function dwA(R::Vector{Float64})
			B₀ = 0;
			A₃ = sin(2*π*(R⋅p.A[:,1])/(p.A[:,1]⋅p.A[:,1]))
			return B₀*[0; 0; A₃]
		end
		return dwA
	elseif(runtype=="eggcarton")
		# egg-carton type potential
		function ecA(R::Vector{Float64})
			B₀ = 1;
			A₂ = 1;
			for i = 1:2
				A₂ *= sin(2*π*(R⋅p.A[:,i])/(p.A[:,i]⋅p.A[:,i]))
			end
			return B₀*[0; A₂; 0]
		end
		return ecA
	elseif(runtype=="bulk")
		function bA(R::Vector{Float64})
			return [0;0;0]
		end
		return bA
	else
		function noA(R::Vector{Float64})
			return [0;0;0]
		end
		return noA
	end
end

A = Agen(runtype,10^4*4*μₑ)


main(p,A)


end
