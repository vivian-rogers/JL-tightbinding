push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")

module SLutils
using Constants
using LinearAlgebra
using Arpack
using Optim

export pruneHoppingType, Agen, βgen, findBandgap

function pruneHoppingType(runtype::String="")
	println("runtype = $runtype")
        if(runtype ∈ ["nanopillars", "eggcarton", "afmthinfilm", "fmthinfilm", "fmdotsP", "fmdotsAP", "neelwall", "blochwall" ])
		return []
		#return ["z"]
	end
	if(runtype == "domainwallthin")
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
			return Float64.(Aval)
		end
		return npA
	elseif(runtype=="afmthinfilm")
		# egg-carton type potential
		function afA(R::Vector{Float64})
			λ = 2*nm; B₀ = 200 # tesla
			#By = cos(R[1]*π/p.a₁[1])*cos(R[2]*π/p.a₂[2])*(1/2)^(-(i-R[3] + p.SL₃[3])*λ)
			A₃ = -(p.a₁[1]/π)*sin(R[1]*π/p.a₁[1])*cos(R[2]*π/p.a₂[2])*(1/2)^(-(R[3]-p.SLa₃[3])/λ)
			#A₃ = sin(2*π*(R⋅p.A[:,1])/(p.A[:,1]⋅p.A[:,1]))
			#println("R₀ = $R")
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
			B₀ = 4000 # Tesla
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


function βgen(p,runtype::String,β₀::Float64=0.2*eV, θ::Float64=360, startDWs::Float64=-5*nm)
	C = ħ/m₀
        β₀ = C^-1 * β₀
	if(runtype=="fmdotsP")
		function fmdotPβ(R::Vector{Float64})
			rad = 0.25; λ = 2*nm
			SLa₁ = p.A[:,1]; SLa₂ = p.A[:,2]
                        coeff = C^-1*β₀*exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)*[1;0;0]
                        #println("Coeff = $(round.(C*coeff,sigdigits=3)), R = $R")
                        for i = 0:1
				for j = 0:1
					R₀ = i*SLa₁ + j*SLa₂
					if( (R[1]-R₀[1])^2 + (R[2]-R₀[2])^2 < (0.25*norm(SLa₂))^2)
						return 1*coeff
					end
				end
			end
			if( (R[1]-0.5*SLa₁[1])^2 + (R[2]-0.5*SLa₂[2])^2 < (0.25*norm(SLa₂))^2)
				return 1*coeff
			else
				return 0*coeff
			end
		end
                return fmdotPβ
	elseif(runtype=="multiblochdws")
            function mbdwβ(R::Vector{Float64})
                θ = 0
                error = 1
                # see: O'Handley on domain walls
                #a = 5*nm; l = 1*nm
                #a = 30*nm; l = 5*nm
                i = 0
                decay= exp(-((p.nz-1)*p.a₃[3] - R[3])/p.λ)
                while error > 10^-5
                    θ₋ = deepcopy(θ)
                    θ += 2*atan(π*exp((p.startDWs - R[1] - p.DWspacing*i - p.DWwidth)/p.DWwidth))
                    i += 1
                    error = abs(θ-θ₋)
                end
                return β₀*[0;cos(θ);sin(θ)]*decay
            end
            return mbdwβ
	elseif(runtype=="fmdotsAP")
            function fmdotAPβ(R::Vector{Float64})
			rad = 0.25; λ = 2*nm
			SLa₁ = p.A[:,1]; SLa₂ = p.A[:,2]
                        coeff = C^-1*β₀*exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)*[1;0;0]
                        #println("Coeff = $(round.(coeff,sigdigits=3)), R = $R")
                        for i = 0:1
				for j = 0:1
					R₀ = i*SLa₁ + j*SLa₂
					if( (R[1]-R₀[1])^2 + (R[2]-R₀[2])^2 < (0.25*norm(SLa₂))^2)
						return 1*coeff
					end
				end
			end
			if( (R[1]-0.5*SLa₁[1])^2 + (R[2]-0.5*SLa₂[2])^2 < (0.25*norm(SLa₂))^2)
				return -1*coeff
			else
				return 0*coeff
			end
		end
		return fmdotAPβ
	elseif(runtype=="neelwall")
		function dwnβ(R::Vector{Float64})
			α = 51 # arbitrary constant to smooth square wave
			λ = 2*nm # penetration depth of ferromagnetism into slab
                        xmag = cos(2*π*(θ/360)*(R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))^α
                        ymag = (1-xmag^2)*sign(sin(2*π*(θ/360)*R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
			decay= C^-1*exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)
			return β₀*[xmag;ymag;0]*decay
		end
		return dwnβ
	elseif(runtype=="simplejunction")
		function sjβ(R::Vector{Float64})
			α = 51 # arbitrary constant to smooth square wave
                        half = p.A[1,1]/2
			λ = 2*nm # penetration depth of ferromagnetism into slab
			decay= exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)
                        if(R[1] > half*(1-0.00001))
                            return β₀*[0;cosd(θ);sind(θ)]*decay
                        else
                            return β₀*[0;1;0]*decay
                        end
		end
		return sjβ
	elseif(runtype=="blochdw")
		function dwbβ(R::Vector{Float64})
			α = 51 # arbitrary constant to smooth square wave
			λ = 2*nm # penetration depth of ferromagnetism into slab
                        ymag = cos(2*π*(θ/360)*(R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
                        zmag = (1-ymag^2)*sign(sin(2*π*(θ/360)*R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
			decay=1
                        #decay= C^-1*exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)
			return β₀*[0;ymag;zmag]*decay
		end
		return dwbβ
	elseif(runtype=="neellattice")
		function latnβ(R::Vector{Float64})
			α = 51 # arbitrary constant to smooth square wave
			λ = 2*nm # penetration depth of ferromagnetism into slab
                        ymag = cos(2*π*(p.θ/360)*(R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
                        xmag = sin(2*π*(p.θ/360)*(R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
                        decay= exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)
			return β₀*[xmag;ymag;0]*decay
		end
		return latnβ
	elseif(runtype=="blochlattice")
		function latbβ(R::Vector{Float64})
			α = 51 # arbitrary constant to smooth square wave
			λ = 2*nm # penetration depth of ferromagnetism into slab
                        ymag = cos(2*π*(p.θ/360)*(R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
                        zmag = sin(2*π*(p.θ/360)*(R⋅p.A[:,1]/(p.A[:,1]⋅p.A[:,1])))
                        decay= exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)
			return β₀*[0;ymag;zmag]*decay
		end
		return latbβ
	elseif(runtype=="fmthinfilm")
		function fmβ(R::Vector{Float64})
			λ = 5*nm
			decay= C^-1*exp(-((p.nz-1)*p.a₃[3] - R[3])/λ)
			return β₀*[0;1;0]*decay
		end
		return fmβ
	else
		function noβ(R::Vector{Float64})
			return [0;0;0]
		end
	end
end
	#main(params)


function findBandgap(p::NamedTuple)
    startk = [0.1;0.01].*(rand(2) .-0.5) # yes this is a bit jank,,,
    function bandgap(H::Function)
            # assumes that the gap lies around E = 0
            function gap(k::Vector{Float64})
                    newk = p.B*vcat(k,[0])
                    Energies, Eigstates = eigs(H(newk), nev=1, maxiter=100000, which=:SM) 
                    Energies = real.(Energies)
                    
                    #assume particle-hole symmetry...
                    ΔE = 2*minimum(abs.(Energies))
                    #conduction = Energies[findall(E->(real(E)>=0),Energies)]
                    #valence = Energies[findall(E->(real(E)<=0),Energies)]
                    #println("C = $conduction")
                    #println("V = $valence")
                    #ΔE = abs(minimum(conduction)-maximum(valence))
                    
                    #println("ΔE = $ΔE")
                    return ΔE
            end
            lower = -0.5*[1;1]; upper=-lower
            #results = optimize(gap,startk)
            #method=NelderMead()
            f_tol = 10^-7
            method=NelderMead()
            #method=SimulatedAnnealing(temperature=(f(i)=(0.2)^i))

            results = optimize(gap,startk,method,Optim.Options(f_tol=10^-7, iterations=1000)
                              )
            #results = optimize(gap,startk,NelderMead())
            #results = optimize(gap,lower,upper,NelderMead())
            #results = optimize(gap,lower,upper,ParticleSwarm(lower=lower,upper=upper,n_particles=5))
            display(results)
            return minimum(results)
    end
end

end
