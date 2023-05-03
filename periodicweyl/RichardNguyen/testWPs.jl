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

function berry_monopoles(H, kpoints)
    # Calculate the Berry connection and curvature
    Nk = size(kpoints, 2)
    Ux, Uy, Uz = zeros(Nk, Nk), zeros(Nk, Nk), zeros(Nk, Nk)
    Fx, Fy, Fz = zeros(Nk, Nk), zeros(Nk, Nk), zeros(Nk, Nk)
    for i in 1:Nk
        for j in 1:Nk
            if i == j
                Ux[i, j] = 1.0
                Uy[i, j] = 1.0
                Uz[i, j] = 1.0
            else
                k = kpoints[:, j] - kpoints[:, i]
                k /= norm(k)
                Hk = H(kpoints[:, i])
                Hk_prime = H(kpoints[:, j])
                Uk = eigvecs(Hk)
                Ukp = eigvecs(Hk_prime)
                Ux[i, j] = dot(Uk[:, 1], Ukp[:, 1])
                Uy[i, j] = dot(Uk[:, 2], Ukp[:, 2])
                Uz[i, j] = dot(Uk[:, 3], Ukp[:, 3])
                Fx[i, j] = imag(log(Ux[i, j])) / (2*pi)
                Fy[i, j] = imag(log(Uy[i, j])) / (2*pi)
                Fz[i, j] = imag(log(Uz[i, j])) / (2*pi)
            end
        end
    end
    
    # Calculate the monopoles of the Berry curvature
    Nc = 0
    for i in 1:Nk
        for j in (i+1):Nk
            for k in (j+1):Nk
                r = Fx[i, j] + Fy[j, k] + Fz[k, i]
                r += Fx[j, k] + Fy[k, i] + Fz[i, j]
                r = round(r)
                if r != 0
                    Nc += sign(r)
                end
            end
        end
    end
    
    return Nc
end
N = 20
kxs = LinRange(-π,π,N)
kys = LinRange(-π,π,N)
kzs = LinRange(-π,π,N)
kpoints = [[kxs[i],kys[j],kzs[k]] for i = 1:N, j = 1:N, k = 1:N]
kpoints = reduce(vcat,kpoints)
print(berry_monopoles(H, kpoints))
