push!(LOAD_PATH,"./")
push!(LOAD_PATH,"./src/")


using LinearAlgebra
using SparseArrays
using UsefulFunctions


n = 150
A = sparse(rand(ComplexF64,n,n))

AexactInv = exactInv(A)
Ainv = grInv(A,5)

println("Difference between A exact and A recInv = $(norm(AexactInv-Ainv)))")

display(Ainv)
