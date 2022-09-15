
module Constants

hartree = false
if(hartree)
	#hartree units
	ħ = 1 #hbar
	h = ħ * 2*π
	m₀ = 1 #electron mass
	q = 1 #electron charge
	ϵ₀ = 1 #permittivity
	r₀ = 1 #bohr radii
	Ry = 1/2 #rydberg
	au = 1 #hartree energy
	Å = 1.8897 * r₀	#angstrom, in hartree units
	nm = 10 * Å 		#angstrom, in hartree units
	qHo = 67;
	eV = au/27.2
	μB = 1
else
	ħ = 1.05457E-34
	h = ħ * 2*π
	m₀ = 9.10938E-31
	q = 1.60218E-19
	ϵ₀ = 8.854E-12
	metre = 1
	au = 27.2
	eV = 1.0
	cm = 1E-2
	μ₀ = 1.2566*10^-6 # H/m vacuum permeability
	nm = 1E-9
	Å = 1E-10
	r₀ = 5.29E-11
	Ry = m₀*q^4 / (8*h^2*ϵ₀^2)
	μₑ = 9.28*10^-24 # electron magnetic moment in A*m^2
	μB = 5.788838E-5 # bohr magneton in eV/T
	kB = 8.617E-5 # boltzmann constant in eV/K
end

#unit conversions
μₑ = 9.28*10^-24 # electron magnetic moment in A*m^2
#Å = 1.8897 * r0 	#angstrom, in hartree units
nm = 10 * Å 		#angstrom, in hartree units
kT = 0.02585*eV 	#DUBIOUS!!! kT @ 293K in hartree units??
#metre = 10^10 * Å
vf = 10^6
# for Si in particular
m = 0.25*m₀
ϵ = 11.9*ϵ₀




#exports all units
for n in names(@__MODULE__; all=true)
               if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
                   @eval export $n
               end
end

end
