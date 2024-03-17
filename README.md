This repository contains the Python infrastructure used to train 'NbCr_Heaton.eam.alloy', an embedded-atom model (EAM) interatomic potential used in 'A Study of the Thermomechanical Behavior of NbCr 
Solid Solutions and the Stability of NbCr2 Laves Phases Using Molecular Dynamics Simulations' by Heaton and Samin. 

'fitNbCr.py' is the code used to train the potential. The code applies Python's non-linear least squares fitting algorithm to optimize all 42 parameters, 21 for each material. All DFT data is contained
in the same file and is labeled as 'DFTdata'. 'runLAMMPS.py' provides the model predictions to 'fitNbCr.py' and makes use of the elastic coefficient example scripts developed by Aidan Thompson. 
