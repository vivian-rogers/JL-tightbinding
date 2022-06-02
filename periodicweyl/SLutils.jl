push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")

module SLutils
using Constants
using LinearAlgebra

export pruneHoppingType, Agen

function pruneHoppingType(runtype::String="")
	if(runtype == "nanopillars" || "eggcarton")
		#return []
		return ["z"]
	end
	if(runtype == "domainwall")
		return ["z","y"]
	end
	if(runtype == "device")
		return ["x","y","z"]
	end
	return []
end

function Agen(p,runtype::String="",M₀::Float64=0) # modify vector potential as desired
	if(runtype=="nanopillars")
		function npA(R::Vector{Float64})
			Aval = zeros(3)
			n = 20
			for i = -n:(n+1)
				for j = -n:(n+1)
					R₀ = i*p.SLa₁ .+ j*p.SLa₂ .+ p.SLa₃ .+ 30*p.a₃
					#println("R₀ = R₀")
					#println("R₀ = $(round.(R₀,sigdigits=3)), R = $(round.(R,sigdigits=3))")
					Aval += μ₀/(4*π) * (M₀*[0;0;1])×(R-R₀)/((R-R₀)⋅(R-R₀))^(3/2)
				end
			end
			#println("A = $Aval")
			#return [0;0;0]
			return Aval
		end
		return npA
	elseif(runtype=="afmthinfilm")
		# egg-carton type potential
		function afA(R::Vector{Float64})
			λ = 2*nm; B₀ = 200 # tesla
			#By = cos(R[1]*π/p.a₁[1])*cos(R[2]*π/p.a₂[2])*(1/2)^(-(i-R[3] + p.SL₃[3])*λ)
			A₃ = -(p.a₁[1]/π)*sin(R[1]*π/p.a₁[1])*cos(R[2]*π/p.a₂[2])*(1/2)^(R[3]/λ)
			#A₃ = sin(2*π*(R⋅p.A[:,1])/(p.A[:,1]⋅p.A[:,1]))
			println("R₀ = $R")
			return B₀*[0; 0; A₃]
		end
		return afA
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
			B₀ = 8*100;
			A₂ = norm(p.SLa₁)/(2*π);
			for i = 1:2
				A₂ *= sin(2*π*(R⋅p.A[:,i])/(p.A[:,i]⋅p.A[:,i]))
			end
			return B₀*[0; 0; A₂]
		end
		return ecA
	elseif(runtype=="bulk")
		function bA(R::Vector{Float64})
			B₀ = 1200 # Tesla
			return [0;0;B₀*R[2]] # B₀ in +x
		end
		return bA
	else
		function noA(R::Vector{Float64})
			return [0;0;0]
		end
		return noA
	end
end

#main(params)

end
