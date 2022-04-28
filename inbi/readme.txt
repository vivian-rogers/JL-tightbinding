Hi!


I wrote these codes to try to model the behavior of indium bismuth, but I suppose this framework can be expanded to many systems. 
This is partially based on some old codes I wrote a while ago, anid there are a bunch of useful functions I spent time tuning in ../src/
Such as plotting, band structures with projection, LDOS, and more

to run, define the routines you want to run in Driver.jl under main() and run "julia -i InBi.jl" 

To understand how the code works, the runSCF function calls the genNNs function which will generate all of the hopping terms as defined.
These are added up into "static" terms and "periodic boundaries" terms, which are connected by bloch's theorem in the H(k) definition. 
runSCF returns a function H(k); julia is at least partially a functional language so this is useful to pass the H(k) around the code.
This can be passed to the bands routine, DOS routine, etc, though it may need tuning. 

Material parameters are defined in InBi.jl, though the notation is certainly lacking for now....



############### INSTALL ######################

How do? 

to install julia on ubuntu, use "sudo apt-get install julia"
on arch based systems, use "yay -Syu julia-git" because the community repo install is buggy. Install yay if not yet installed.

from there, one can type "julia" then "]" and type:
add PyPlot (requires working matplotlib/numpy/etc etc python3 install on machine)
add Plots
add LinearAlgebra
add ColorSchemes
add Arpack
add SparseArrays
add Printf

I also suggest installing julia-vim for manipulating documents with the nice unicode, if you do like using vim. Else, I suggest finding an IDE or text editor that supports unicode. I can change the code, but the fancy symbols make it easier I think.

These things would be automated by running "bash init.sh", which you may also look at.  


Truly, you (we (I)) should be using an IDE or jupyter notebooks. I'd like to set that up at some point, if you can figure it out quickly do let me know!
