σ₁ = [0 1; 1 0]; σ₂ = [0 -im; im 0]; σ₃ = [1 0 ; 0 -1];
τ₁ = σ₁; τ₂ = σ₂; τ₃ = σ₃;
σ = Vector{Matrix}(undef,3); σ[1] = σ₁; σ[2] = σ₂; σ[3] = σ₃;
⊗(a,b) = kron(a,b)

function matdot(A::Vector{Matrix},B::Vector{Float64})
    sum = zeros(ComplexF64,size(A[1]))
    for i in eachindex(B)
        sum .+= B[i]*A[i]
    end
    return sum
end

function H(k::Vector{Float64})
    vf = 1; β = 1; m = 0.5;
    return vf*τ₁⊗matdot(σ,k) .+ m*τ₃⊗I(2) .+ β*I(2)⊗σ₃ 
end


