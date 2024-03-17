variable up equal 1.0e-6
variable atomjiggle equal 1.0e-5
units metal
variable cfac equal 1.0e-4
variable cunits string GPa
variable etol equal 0.0 
variable ftol equal 1.0e-10
variable maxiter equal 100
variable maxeval equal 1000
variable dmax equal 1.0e-2
boundary p p p
lattice fcc 3.62 origin 0.0 0.0 0.0 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region box prism 0 6 0 6 0 6 0.0 0.0 0.0 units lattice
create_box 2 box
create_atoms 2 box
