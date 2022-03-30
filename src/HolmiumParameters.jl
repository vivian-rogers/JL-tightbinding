module HolmiumParameters

using Constants

a = 3.57*Å
c = 5.64*Å
t = 2*eV
ε = 1*eV
a₁ = a*[cosd(-30);sind(-30);0]
a₂ = a*[cosd(30);sind(30);0]
a₃ = c*[0;0;1]

r = sqrt((tand(30)*a/2)^2 + (a/2)^2)
r₁ = r*[cosd(30); sind(30); 0]

x₂ = r/(2*cosd(3)); y₂= √(r^2 - x₂^2);
r₂ = [-x₂;0; y₂]
r₃ = [x₂;0; -y₂]

#exports all units
for n in names(@__MODULE__; all=true)
               if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
                   @eval export $n
               end
end

end
