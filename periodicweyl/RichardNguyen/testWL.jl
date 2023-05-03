using LinearAlgebra

σx = [0 1; 1 0]
σy = [0 -im; im 0]
σz = [1 0; 0 -1]
# Define the Hamiltonian H(k)
function H(k)
    kx = k[1]; ky = k[2]; kz = k[3]
    # Define the Pauli matrices

    # Define the Hamiltonian
    H = sin(kx)*σx + sin(ky)*σy + sin(kz)*σz

    return H
end

# Define the Brillouin zone
N = 1000
#dK = 2π/N
ΔK = 0.001*π

Kx = LinRange(-ΔK, ΔK, N)
Ky = zeros(N)
Kz = LinRange(-ΔK, ΔK, N)

dK = 2*ΔK/(N)
kpath1 = [(Kx[i] , ΔK,  -ΔK) for i in 1:N]
kpath2 = [(ΔK    , ΔK, Kz[i] ) for i in 1:N]
kpath3 = [(-Kx[i], ΔK, ΔK    ) for i in 1:N]
kpath4 = [(-ΔK     , ΔK, -Kz[i]) for i in 1:N]

kpoints = hcat(kpath1,kpath2,kpath3,kpath4)

σ = [σx,σy,σz]
# Define the number of time steps

# Initialize the Wilson loop matrix

W=I(2)

eigvecs_list = []
for i in 1:(4*N)
    vecs = eigvecs(H(kpoints[i]))
    push!(eigvecs_list,vecs)
end

for i in 1:(4*N)
    k = kpoints[i]
    #show(k)
    #println()
    U_pre = eigvecs_list[mod(i-2,N)+1]
    U_0 = eigvecs_list[i]
    U_post = eigvecs_list[mod(i,N)+1]
    # Compute the gauge potential
    
    #U_k = exp(im*dK*(U_0)*(U_post-U_pre)')
    # Update the Wilson loop matrix
    global W *= U_0'*U_post
end

#=
W=I(2)

# Loop over all k-points in the path
for i in 1:N
    k = kpoints[i]

    # Compute the Hamiltonian and its eigenvectors
    Hk = H(k)
    vecs = eigvecs(Hk)

    # Loop over all eigenvectors
    for l in 1:2
        # Compute the gauge potential
        dU = zeros(2,2)
        #for m in 1:3
        U = exp(-im*(k[1]*σ[1]+k[2]*σ[2]+k[3]*σ[3])) * vecs[:,l] * vecs[:,l]'
        dU += im * U * σ[m] * U'
        #end

        # Compute the parallel transport matrix
        U_k = exp(-im*dK*Hk) * exp(dK*dU)

        # Update the Wilson loop matrix
        global W *= U_k
    end
end
=#
# Compute the trace of the Wilson loop
Tr_W = real(tr(W))

# Print the trace of the Wilson loop
println("The Wilson loop is: ")
display(W)
println("The trace of the Wilson loop is: ", Tr_W)

