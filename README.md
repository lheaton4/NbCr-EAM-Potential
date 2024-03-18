This repository contains the Python coding infrastructure used to train 'NbCr_Heaton.eam.alloy', the embedded-atom model (EAM) interatomic potential developed and used in 'A Study of the Thermomechanical Behavior of NbCr Solid Solutions and the Stability of NbCr2 Laves Phases Using Molecular Dynamics Simulations' by Lucas A. Heaton and Adib J. Samin. The crack propagation script ('in.crack') employed in the work is also included. The file simulates crack propagation in a pure Nb matrix, reading in 'Nb_large.lmp'.

'fitNbCr.py' is the only script that needs to be executed. The code applies Python's non-linear least squares fitting algorithm to optimize all 42 parameters based on model predictions generated by 'runLAMMPS.py'. Upon running the code, a 'CrNb.eam.alloy' file will be generated which is the setfl formatted file that all LAMMPS runs reference. It is updated with each new iteration according to the new parameters. In addition to this file, a number of output files and LAMMPS scripts are generated. These are clearly labeled and may be deleted after optimization. Eight LAMMPS inputs scripts external to 'runLAMMPS.py' are used to compute the lattice constants and cohesive energies of the four structures used for training. These can be easily written in 'runLAMMPS.py' instead, if so desired. The C11, C12, and C44 elastic coefficients are calculated using the input scripts developed by Aidan Thompson in the LAMMPS distribution. The DFT data used to perform the fits is contained within 'fitNbCr.py' and is labeled as 'DFTdata'.



