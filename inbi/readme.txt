Hi!


I wrote these codes to try to model the behavior of indium bismuth, but I suppose this framework can be expanded to many systems. 
This is partially based on some old codes I wrote a while ago, and there are a bunch of useful functions I spent time tuning in ../src/
Such as plotting, band structures with projection, LDOS, and more

to run, define the routines you want to run in Driver.jl under main() and run "julia -i InBi.jl" 

To understand how the code works, the runSCF function calls the genNNs function which will generate all of the hopping terms as defined.
These are added up into "static" terms and "periodic boundaries" terms, which are connected by bloch's theorem in the H(k) definition. 
runSCF returns a function H(k); julia is at least partially a functional language so this is useful to pass the H(k) around the code.
This can be passed to the bands routine, DOS routine, etc, though it may need tuning. 

Material parameters are defined in InBi.jl, though the notation is certainly lacking for now....
