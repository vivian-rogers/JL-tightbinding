
using LinearAlgebra
using SparseArrays

A_gen(N) = rand(ComplexF64,N,N)
B_gen(N) = rand(ComplexF64,N,N)

function f_prepostproc(N)
   f = f_generator(N)
   println("postproc finished")
   return f
end

function f_generator(N)
   A = A_gen(N)
   B = B_gen(N)
   function f(x)
      return sparse(A .+ B*x) #
  end
  println("f is complete, returning")
  return f
end 

function dostuff(f)
	println("entering dostuff")
	example = tr(f(1))
	println("Ex: $example")
	return
end


function main(N)
   println("Running with N = $N")
   f = f_prepostproc(N)
   dostuff(f)
end



N = 100^2
main(N)
