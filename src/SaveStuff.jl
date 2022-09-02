module SaveStuff

using DelimitedFiles
using JLD2

function mktxt(name::String,contents::String)
    open(name, "w") do file
        write(file, contents)
    end
    println("Saved output to $name")
end

function mkdelim(name::String, data) # put in data as [x1 x2 ...]
    open(name, "w") do io
        writedlm(io, data,)
    end
    println("Saved output to $name")
end

# export all 
for n in names(@__MODULE__; all=true)
   if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
       @eval export $n
   end
end


end

