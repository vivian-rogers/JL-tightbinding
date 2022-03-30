module GMaterialParameters

using Constants
using UsefulFunctions 
#a = 1
# for ScN



#define parameters to pass around

params =(

	#real-space lattice vectors
	
	a = 4.5*Å,
	a₁ = a*(1/2)*[0;1;1], a₂ = a*(1/2)*[1;0;1], a₃ = a*(1/2)*[1;1;0],
	A = hcat(a₁, a₂, a₃),
	B = transpose((2*π)*inv(A)),

	#reciprocal lattice
	b₁ = B[:,1], b₂ = B[:,2], b₃ = B[:,3],

	# define atoms
	atoms = [
			 Atom("Sc", [0,0,0], ["dxy","dxz","dz2","dyz","dx2-y2"], 21,1),
			 Atom("N", (a₁ + a₂ +a₃)/2, ["px","py","pz"], 7,3),
			 #Atom("Sc", [0,0,0], [-1,-1,-1,-1,0,0,0,0,0], 21,1),
			 #Atom("N", (a₁ + a₂ +a₃)/2, [-1,0,0,0,1,1,1,1,1],7,3),
			 ],
	
	orbitals = [
							 "s",
							 "px",
							 "py",
							 "pz",
							 "dxy",
							 "dxz",
							 "dz2",
							 "dyz",
							 "dx2-y2",
		 ],
							 

	# Nearest-neighbor vector collections
	NNs = [nn1_vects(a),nn2_vects(a)],
	
	# The atoms whose wfcs overlap with each other
	connections = [
					Connection("Sc","N",1,1),
					Connection("Sc","Sc",1,2),
					Connection("Sc","Sc",2,1),
					],
					

	# orbital eigenbasis
	orbitalEigenBasis = Dict(
							 "s" => 1,
							 "px" => 2,
							 "py" => 3,
							 "pz" => 4,
							 "dxy" => 5,
							 "dxz" => 6,
							 "dz2" => 7,
							 "dyz" => 8,
							 "dx2-y2" => 9
							 ),
	
	atomicEigenBasis = Dict(
							"Sc"=>1,
							"N"=>2,
							),

	kdict = Dict(
			"γ" => B*[0;     0],
			"κ" => B*[2/3; 1/3],
			"κ'" => B*[-2/3; -1/3],
			"κ'2" => B*[1/3; -2/3],
			"μ" => B*[1/2; 0],
			"μ'" => B*[-1/2; 0],
			),
	)

