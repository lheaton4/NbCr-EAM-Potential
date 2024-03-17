import numpy as np
import os

def remove_first_blank_line(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        
    if lines and lines[0].strip() == '':
        lines.pop(0)

    with open(filename, 'w') as file:
        file.writelines(lines)

def computeParameters():
    
    solutionVector = np.zeros(60)
    
    def calcLatticeParameterCr_BCC():
            
        os.system('lmp -in in.latticeConstantCr_BCC > latticeConstantCr_BCC_output.txt')
    
        with open('latticeConstantCr_BCC_output.txt', 'r') as f:
            output_lines = f.readlines()
            last_line = output_lines[-2].strip()
            lattice_constant = float(last_line)
            
        return lattice_constant
    
    def calcCohesiveEnergyCr_BCC():
        
        os.system('lmp -in in.cohesiveEnergyCr_BCC > cohesiveEnergyCr_BCC_output.txt')
    
        with open('cohesiveEnergyCr_BCC_output.txt','r') as f:
            output_lines = f.readlines()
            last_line = output_lines[-2].strip()
            cohesive_energy = float(last_line)
            
        return cohesive_energy
    
    def calcElasticCoeffCr_BCC():
        lammpsScript = '''        
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
lattice bcc 2.846 origin 0.0 0.0 0.0 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region box prism 0 6 0 6 0 6 0.0 0.0 0.0 units lattice
create_box 2 box
create_atoms 2 box
'''
        with open('init.mod', 'w') as f:
            f.write(lammpsScript)
            
        remove_first_blank_line('init.mod')
    
        os.system('lmp -in in.elastic.lmp > elasticCoeffCr_BCC_output.txt')
    
        with open('elasticCoeffCr_BCC_output.txt', 'r') as f:
            output_lines = f.readlines()
            c11_bcc = float(output_lines[-29].strip()[:-4][26:])
            c12_bcc = float(output_lines[-26].strip()[:-4][26:])
            c44_bcc = float(output_lines[-23].strip()[:-4][26:])
            
        return (c11_bcc,c12_bcc,c44_bcc)
    
    latticeConstantCr_bcc = calcLatticeParameterCr_BCC()
    cohesiveEnergyCr_bcc = calcCohesiveEnergyCr_BCC()
    c11_Cr_bcc, c12_Cr_bcc, c44_Cr_bcc = calcElasticCoeffCr_BCC()
    
    solutionVector[:5] = [latticeConstantCr_bcc,cohesiveEnergyCr_bcc,c11_Cr_bcc,c12_Cr_bcc,c44_Cr_bcc]
    constants = np.linspace(2.3,3.2,10)
    
    def calcEnergyVolumeCr_bcc(a):
        lammpsScript = '''
clear
units metal
dimension 3
boundary p p p
lattice bcc '''+str(a)+''' origin 0.0 0.0 0.0 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region box block 0 6 0 6 0 6 units lattice
create_box 2 box
create_atoms 2 box
pair_style eam/alloy
pair_coeff * * CrNb.eam.alloy Nb Cr
minimize 1.0e-8 1.0e-10 100000 1000000
variable potEnergy equal "pe"
variable natoms equal "count(all)"
variable PEperAtom equal "v_potEnergy/v_natoms"
print "${PEperAtom}"
'''
        
        with open('in.energyVolumeCr_BCC'+str(a), 'w') as f:
            f.write(lammpsScript)
            
        remove_first_blank_line('in.energyVolumeCr_BCC'+str(a))
            
        return
    
    count = 5
    
    for item in constants:
        
        calcEnergyVolumeCr_bcc(item)
        os.system('lmp -in in.energyVolumeCr_BCC'+str(item)+' > energyVolumeCr_BCC'+str(item)+'_output.txt')
        with open('energyVolumeCr_BCC'+str(item)+'_output.txt', 'r') as f:
            output_lines = f.readlines()
            last_line = output_lines[-2].strip()
            solutionVector[count] = float(last_line)-cohesiveEnergyCr_bcc
        count += 1
        
    def calcLatticeParameterCr_FCC():
            
        os.system('lmp -in in.latticeConstantCr_FCC > latticeConstantCr_FCC_output.txt')
    
        with open('latticeConstantCr_FCC_output.txt', 'r') as f:
            output_lines = f.readlines()
            last_line = output_lines[-2].strip()
            lattice_constant = float(last_line)
            
        return lattice_constant
    
    def calcCohesiveEnergyCr_FCC():
        
        os.system('lmp -in in.cohesiveEnergyCr_FCC > cohesiveEnergyCr_FCC_output.txt')
    
        with open('cohesiveEnergyCr_FCC_output.txt','r') as f:
            output_lines = f.readlines()
            last_line = output_lines[-2].strip()
            cohesive_energy = float(last_line)
            
        return cohesive_energy
    
    def calcElasticCoeffCr_FCC():
        lammpsScript = '''    
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
'''
        with open('init.mod', 'w') as f:
            f.write(lammpsScript)
            
        remove_first_blank_line('init.mod')

    
        os.system('lmp -in in.elastic.lmp > elasticCoeffCr_FCC_output.txt')
    
        with open('elasticCoeffCr_FCC_output.txt', 'r') as f:
            output_lines = f.readlines()
            c11_fcc = float(output_lines[-29].strip()[:-4][26:])
            c12_fcc = float(output_lines[-26].strip()[:-4][26:])
            c44_fcc = float(output_lines[-23].strip()[:-4][26:])
            
        return (c11_fcc,c12_fcc,c44_fcc)
    
    latticeConstantCr_FCC = calcLatticeParameterCr_FCC()
    cohesiveEnergyCr_FCC = calcCohesiveEnergyCr_FCC()
    c11_Cr_FCC, c12_Cr_FCC, c44_Cr_FCC, = calcElasticCoeffCr_FCC()
    
    solutionVector[15:20] = [latticeConstantCr_FCC,cohesiveEnergyCr_FCC,c11_Cr_FCC,c12_Cr_FCC,c44_Cr_FCC]
    # solutionVector[5:10] = [latticeConstantCr_FCC,cohesiveEnergyCr_FCC,c11_Cr_FCC,c12_Cr_FCC,c44_Cr_FCC]
    
    constants = np.linspace(3.1,4,10)
    
    def calcEnergyVolumeCr_FCC(a):
        lammpsScript = '''
clear
units metal
dimension 3
boundary p p p
lattice fcc '''+str(a)+''' origin 0.0 0.0 0.0 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region box block 0 6 0 6 0 6 units lattice
create_box 2 box
create_atoms 2 box
pair_style eam/alloy
pair_coeff * * CrNb.eam.alloy Nb Cr
minimize 1.0e-8 1.0e-10 100000 1000000
variable potEnergy equal "pe"
variable natoms equal "count(all)"
variable PEperAtom equal "v_potEnergy/v_natoms"
print "${PEperAtom}"
'''
        
        with open('in.energyVolumeCr_FCC'+str(a), 'w') as f:
            f.write(lammpsScript)
            
        remove_first_blank_line('in.energyVolumeCr_FCC'+str(a))
            
        return
    
    count = 20
    
    for item in constants:
        
        calcEnergyVolumeCr_FCC(item)
        os.system('lmp -in in.energyVolumeCr_FCC'+str(item)+' > energyVolumeCr_FCC'+str(item)+'_output.txt')
        with open('energyVolumeCr_FCC'+str(item)+'_output.txt', 'r') as f:
            output_lines = f.readlines()
            last_line = output_lines[-2].strip()
            solutionVector[count] = float(last_line)-cohesiveEnergyCr_FCC
        count += 1
        
    def calcLatticeParameterNb_BCC():
            
        os.system('lmp -in in.latticeConstantNb_BCC > latticeConstantNb_BCC_output.txt')
    
        with open('latticeConstantNb_BCC_output.txt', 'r') as f:
            output_lines = f.readlines()
            last_line = output_lines[-2].strip()
            lattice_constant = float(last_line)
            
        return lattice_constant
    
    def calcCohesiveEnergyNb_BCC():
        
        os.system('lmp -in in.cohesiveEnergyNb_BCC > cohesiveEnergyNb_BCC_output.txt')
    
        with open('cohesiveEnergyNb_BCC_output.txt','r') as f:
            output_lines = f.readlines()
            last_line = output_lines[-2].strip()
            cohesive_energy = float(last_line)
            
        return cohesive_energy
    
    def calcElasticCoeffNb_BCC():
        lammpsScript = '''        
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
lattice bcc 3.322 origin 0.0 0.0 0.0 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region box prism 0 6 0 6 0 6 0.0 0.0 0.0 units lattice
create_box 2 box
create_atoms 1 box
'''
        with open('init.mod', 'w') as f:
            f.write(lammpsScript)
            
        remove_first_blank_line('init.mod')
    
        os.system('lmp -in in.elastic.lmp > elasticCoeffNb_BCC_output.txt')
    
        with open('elasticCoeffNb_BCC_output.txt', 'r') as f:
            output_lines = f.readlines()
            c11_bcc = float(output_lines[-29].strip()[:-4][26:])
            c12_bcc = float(output_lines[-26].strip()[:-4][26:])
            c44_bcc = float(output_lines[-23].strip()[:-4][26:])
            
        return (c11_bcc,c12_bcc,c44_bcc)
    
    latticeConstantNb_BCC = calcLatticeParameterNb_BCC()
    cohesiveEnergyNb_BCC = calcCohesiveEnergyNb_BCC()
    c11_Nb_BCC, c12_Nb_BCC, c44_Nb_BCC = calcElasticCoeffNb_BCC()
    
    solutionVector[30:35] = [latticeConstantNb_BCC,cohesiveEnergyNb_BCC,c11_Nb_BCC,c12_Nb_BCC,c44_Nb_BCC]
    constants = np.linspace(2.75,3.65,10)
    
    def calcEnergyVolumeNb_BCC(a):
        lammpsScript = '''
clear
units metal
dimension 3
boundary p p p
lattice bcc '''+str(a)+''' origin 0.0 0.0 0.0 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region box block 0 6 0 6 0 6 units lattice
create_box 2 box
create_atoms 1 box
pair_style eam/alloy
pair_coeff * * CrNb.eam.alloy Nb Cr
minimize 1.0e-8 1.0e-10 100000 1000000
variable potEnergy equal "pe"
variable natoms equal "count(all)"
variable PEperAtom equal "v_potEnergy/v_natoms"
print "${PEperAtom}"
'''
        
        with open('in.energyVolumeNb_BCC'+str(a), 'w') as f:
            f.write(lammpsScript)
            
        remove_first_blank_line('in.energyVolumeNb_BCC'+str(a))
            
        return
    
    count = 35
    
    for item in constants:
        
        calcEnergyVolumeNb_BCC(item)
        os.system('lmp -in in.energyVolumeNb_BCC'+str(item)+' > energyVolumeNb_BCC'+str(item)+'_output.txt')
        with open('energyVolumeNb_BCC'+str(item)+'_output.txt', 'r') as f:
            output_lines = f.readlines()
            last_line = output_lines[-2].strip()
            solutionVector[count] = float(last_line)-cohesiveEnergyNb_BCC
        count += 1
        
    def calcLatticeParameterNb_FCC():
            
        os.system('lmp -in in.latticeConstantNb_FCC > latticeConstantNb_FCC_output.txt')
    
        with open('latticeConstantNb_FCC_output.txt', 'r') as f:
            output_lines = f.readlines()
            last_line = output_lines[-2].strip()
            lattice_constant = float(last_line)
            
        return lattice_constant
    
    def calcCohesiveEnergyNb_FCC():
        
        os.system('lmp -in in.cohesiveEnergyNb_FCC > cohesiveEnergyNb_FCC_output.txt')
    
        with open('cohesiveEnergyNb_FCC_output.txt','r') as f:
            output_lines = f.readlines()
            last_line = output_lines[-2].strip()
            cohesive_energy = float(last_line)
            
        return cohesive_energy
    
    def calcElasticCoeffNb_FCC():
        lammpsScript = '''    
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
lattice fcc 4.23 origin 0.0 0.0 0.0 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region box prism 0 6 0 6 0 6 0.0 0.0 0.0 units lattice
create_box 2 box
create_atoms 1 box
'''
        with open('init.mod', 'w') as f:
            f.write(lammpsScript)
            
        remove_first_blank_line('init.mod')

    
        os.system('lmp -in in.elastic.lmp > elasticCoeffNb_FCC_output.txt')
    
        with open('elasticCoeffNb_FCC_output.txt', 'r') as f:
            output_lines = f.readlines()
            c11_fcc = float(output_lines[-29].strip()[:-4][26:])
            c12_fcc = float(output_lines[-26].strip()[:-4][26:])
            c44_fcc = float(output_lines[-23].strip()[:-4][26:])
            
        return (c11_fcc,c12_fcc,c44_fcc)
    
    latticeConstantNb_FCC = calcLatticeParameterNb_FCC()
    cohesiveEnergyNb_FCC = calcCohesiveEnergyNb_FCC()
    c11_Nb_FCC, c12_Nb_FCC, c44_Nb_FCC, = calcElasticCoeffNb_FCC()
    
    solutionVector[45:50] = [latticeConstantNb_FCC,cohesiveEnergyNb_FCC,c11_Nb_FCC,c12_Nb_FCC,c44_Nb_FCC]
    constants = np.linspace(3.75,4.65,10)
    
    def calcEnergyVolumeNb_FCC(a):
        lammpsScript = '''
clear
units metal
dimension 3
boundary p p p
lattice fcc '''+str(a)+''' origin 0.0 0.0 0.0 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region box block 0 6 0 6 0 6 units lattice
create_box 2 box
create_atoms 1 box
pair_style eam/alloy
pair_coeff * * CrNb.eam.alloy Nb Cr
minimize 1.0e-8 1.0e-10 100000 1000000
variable potEnergy equal "pe"
variable natoms equal "count(all)"
variable PEperAtom equal "v_potEnergy/v_natoms"
print "${PEperAtom}"
'''
        
        with open('in.energyVolumeNb_FCC'+str(a), 'w') as f:
            f.write(lammpsScript)
            
        remove_first_blank_line('in.energyVolumeNb_FCC'+str(a))
            
        return
    
    count = 50
    
    for item in constants:
        
        calcEnergyVolumeNb_FCC(item)
        os.system('lmp -in in.energyVolumeNb_FCC'+str(item)+' > energyVolumeNb_FCC'+str(item)+'_output.txt')
        with open('energyVolumeNb_FCC'+str(item)+'_output.txt', 'r') as f:
            output_lines = f.readlines()
            last_line = output_lines[-2].strip()
            solutionVector[count] = float(last_line)-cohesiveEnergyNb_FCC
        count += 1
            
    return solutionVector