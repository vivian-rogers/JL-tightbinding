push!(LOAD_PATH, "./src/")


using Plots
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators

B = 0.01*[3;2;-4] # B field in hartree units
H = Hₙ⊗I(2) + I(6)⊗(B⋅σ)
E = 27.2*eigvals(H)
index = [i for i = 1:length(E)]
display(plot(index,E, seriestype = :scatter))

#println(N)

