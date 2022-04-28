





first construct the 2D grid

build atomic orbital hamiltonian in basis |a>, |b>, |c>, etc
Hₐ =    [
	ϵ1 0
	0 ϵ2
	]
build hopping term hamiltonian between orbitals between adjacent atoms

Hₕ = [
      t11 t12
      t21 t22
     ]

consider H in basis |a,↑>, |a,↓>, |b,↑>, |b,↓>

H = Hₐ⊗𝟙₂ + B₃*𝟙₂⊗S₃ #atom orbital energies and zeeman splitting in z 

# B > 0 for conical helimagnet, B = 0 for helimagnet below ~130 K

for i = 1:3 
	# r₁ is the vector from central atom to 1 right
	# R is the rotation matrix around Z
	r = R(i*π/3)*r₁
	H += exp(im*k*R(i*π/3)*r)*Hₕ
end

now consider atom1 orbital - atom2 orbital exchange interactions
Not truly interested in accurate hamiltonian, only correct ground state

# orbital interactions inside atom: 
# THIS IS WRONG
H_self = [
      # |a,↑>, |a,↓>, |b,↑>, |b,↓>
	0	Jab	0	0
	Jab	0	0	0
	0	0	0	Jab
	0	0	Jab	0
	] 

#orbital interactions between atom and next atom
# 8 x 8 matrix with orbital exchange couplings, in |atom, orbital, spin > basis
H_nn = [
	0	0		0
	0	0	0	0
	0	0	0	0
	0	0	0	0
	]

  →   →
→   →   →
  →   →  

OR build up eigendecomposition with "correct" ground state wavefunctions

	ground state wfc at atom 1:
        # |a,→>, |a,←>, |b,→>, |b,←>
	wfc_a = [0.95; 0.3; 0.8; 0.5]
	#actual values aren't correct, but should be collinear in x
	wfc at next atom = 
	wfc_b = exp(im*k⋅r)*[0.95; 0.3; 0.8; 0.5]

H_nn = E_ground*|wfc_a, wfc_b><wfc_a, wfc_b| + E_excited*|wfc_excited><wfc2_excited|
with the expectation that the higher states won't be occupied at low temp
	
then next layer is linked by 
	wfc_c = exp(im*k⋅r)* 𝟙₂⊗R(ϕ)*|wfc_a>

where R(ϕ) is the rotation matrix for the spins around z, for a turn of ~30°



To get actual wfcs, we know (say |ψ> = [a;b]. μ ≈ 9 μB for Ho atoms, but will be smaller for orbitals)

θₜ = 10° for conical helimagnet phase, 0 in helimagnet phase
ϕ = 30° for conical hemimagnet, smaller and temp-dependent for helimagnet
<ψ|σ₁|ψ> = μ*cos(ϕ)*cos(θₜ) = b'a + a'b
<ψ|σ₂|ψ> = μ*sin(ϕ)*cos(θₜ) = im*a'b - im*b'a
<ψ|σ₃|ψ> = μ*sin(θₜ) = a'a - b'b

-> solve system of equations to get spins 


	





