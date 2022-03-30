module GrapheneParameters

using Constants

#a = 1
a = 2.46*Å
t = 2.8*eV
ε = 0*eV
a₁ = a*[cosd(-30);sind(-30);0]
a₂ = a*[cosd(30);sind(30);0]


r = sqrt((tand(30)*2.46/2)^2 + (2.46/2)^2)*Å
r₁ = r*[cosd(-30);sind(-30);0]
r₂ = r*[cosd(30);sind(30);0]

#exports all units
for n in names(@__MODULE__; all=true)
               if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
                   @eval export $n
               end
end

end
