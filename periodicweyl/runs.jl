push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")

module runs
using InBi
using Driver


#main(params)

nx = 1; ny = 1; nz = 1; 
# superlattice basis vectors, in basis of a_1, a_2, a_3
SL1 = [nx; 0; 0]; SL2 = [0; ny; 0]; SL3 = [0; 0; nz]

p = genSL(params,nx, ny, nz, SL1, SL3, SL3) # generate SL params

function A(R::Vector{Float64}) # modify vector potential as desired
	B₀ = 0 # 1 Tesla
	x = R[1]; y = R[2]; z = R[3];
	# egg-carton type potential
	A₂ = 0;
	for i = 1:3
		A₂ *= sin(2*π*(R⋅p.A[:,i])/(p.A[:,i]⋅p.A[:,i]))
	end
	return B₀*[0; A₂; 0]
end




main(p,A)


end
