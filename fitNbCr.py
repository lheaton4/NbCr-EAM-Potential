from runLAMMPS import computeParameters
from numpy import *
import os
from scipy.optimize import least_squares

filePath = 'updatedParameters.txt' # Nb comes first in file, followed by Cr
try:
    os.remove(filePath)
    print(f"File '{filePath}' has been deleted.")
except OSError as e:
    print(f"Error deleting the file: {e}")

DFTdata = array([   # bcc Cr (lattice constant, cohesive energy, C11, C12, C44, eqns of state...)
                    2.846, -5.96, 504, 145, 109, 
                    6.09082, 3.5569195, 1.8844245, 0.84168, 0.262358, 
                    0.0231535, 0.028515, 0.206396, 0.50288, 0.877319,
                    # fcc Cr
                    3.62, -5.57, 48, 336, -90, 
                    2.63247525, 1.561497, 0.82586525, 0.355957, 0.09603275,
                    0.00226475, 0.0366755, 0.1686625, 0.3741625, 0.63347075,
                    # bcc Nb
                    3.322, -6.987, 205, 151, 11,
                    4.592245, 2.798216, 1.569646, 0.765903, 0.279655, 
                    0.044635, 0.005409, 0.1180895, 0.344772, 0.658587,
                    # fcc Nb
                    4.23, -6.66 ,-30, 257, -72,
                    1.59098475, 0.93389075, 0.47617675, 0.18262425, 0.03578275,
                    0.0014385, 0.0621165, 0.2021285, 0.4047355, 0.65385775,
                    ])

y2 = DFTdata**2

def atomicElectronDensity11(fE11, beta11, rE11, lamb11, r):
    y = r/rE11
    denominator = 1 + (y - lamb11) ** 20
    numerator = fE11 * exp(-beta11 * (y - 1))
    return numerator / denominator


def atomicElectronDensity22(fE22, beta22, rE22, lamb22, r):
    y = r/rE22
    denominator = 1 + (y - lamb22) ** 20
    numerator = fE22 * exp(-beta22 * (y - 1))
    return numerator / denominator

def phi11(r,A11,alpha11,rE11,kappa11,B11,beta11,lamb11):
    y = r/rE11
    
    numerator1 = A11*exp(-alpha11*(y-1))
    denominator1 = 1 + ((y-kappa11)**20) # change to Nb values
    
    numerator2 = B11*exp(-beta11*(y-1))
    denominator2 = 1 + ((y-lamb11)**20) # change to Nb values
        
    return (numerator1/denominator1) - (numerator2/denominator2)


def phi22(r,A22,alpha22,rE22,kappa22,B22,beta22,lamb22):
    y = r/rE22
    
    numerator1 = A22*exp(-alpha22*(y-1))
    denominator1 = 1 + ((y-kappa22)**20)
    
    numerator2 = B22*exp(-beta22*(y-1))
    denominator2 = 1 + ((y-lamb22)**20)
        
    return (numerator1/denominator1) - (numerator2/denominator2)


def phi12(r,A11,alpha11,rE11,kappa11,B11,beta11,lamb11,A22,alpha22,rE22,kappa22,B22,beta22,lamb22,fE11,fE22):
    term1 = phi11(r,A11,alpha11,rE11,kappa11,B11,beta11,lamb11)*((atomicElectronDensity22(fE22,beta22,rE22,lamb22,r)/
                                                                  atomicElectronDensity11(fE11,beta11,rE11,lamb11,r)))
    term2 = phi22(r,A22,alpha22,rE22,kappa22,B22,beta22,lamb22)*((atomicElectronDensity11(fE11,beta11,rE11,lamb11,r)/
                                                                  atomicElectronDensity22(fE22,beta22,rE22,lamb22,r)))
    return 0.5*(term1+term2)

def embeddingFunction1(rho,Fn0_11,Fn1_11,Fn2_11,Fn3_11,F0_11,F1_11,F2_11,F3_11,FE_11,rhoS_11,aeta_11,rhoE_11):
    if rho == 0:
        return 0
    rhoN_11 = 0.85*rhoE_11
    rho0_11 = 1.15*rhoE_11
    fni = [Fn0_11,Fn1_11,Fn2_11,Fn3_11]
    Fi = [F0_11,F1_11,F2_11,F3_11]
    if rho < rhoN_11:
        F = 0
        for i in range(4):
            F += fni[i]*((rho/rhoN_11)-1)**i
    elif rho >= rhoN_11 and rho < rho0_11:
        F = 0
        for i in range(4):
            F += Fi[i]*((rho/rhoE_11)-1)**i
    elif rho >= rho0_11:
        x = (rho/rhoS_11)**aeta_11
        F = FE_11*(1-log(x))*x
    return F

def embeddingFunction2(rho,Fn0_22,Fn1_22,Fn2_22,Fn3_22,F0_22,F1_22,F2_22,F3_22,FE_22,rhoS_22,aeta_22,rhoE_22):
    if rho == 0:
        return 0
    rhoN_22 = 0.85*rhoE_22
    rho0_22 = 1.15*rhoE_22
    fni = [Fn0_22,Fn1_22,Fn2_22,Fn3_22]
    Fi = [F0_22,F1_22,F2_22,F3_22]
    if rho < rhoN_22:
        F = 0
        for i in range(4):
            F += fni[i]*((rho/rhoN_22)-1)**i
    elif rho >= rhoN_22 and rho < rho0_22:
        F = 0
        for i in range(4):
            F += Fi[i]*((rho/rhoE_22)-1)**i
    elif rho >= rho0_22:
        x = (rho/rhoS_22)**aeta_22
        F = FE_22*(1-log(x))*x
    return F

def remove_first_blank_line(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    with open(filename, 'w') as f:
        for line in lines:
            if line.strip() or lines.index(line) > 0:
                f.write(line)
                
def writeEAM(params):
    (fE11, A11, alpha11, rE11, kappa11, B11, beta11, lamb11, Fn0_11, Fn1_11, 
      Fn2_11, Fn3_11, F0_11, F2_11, F3_11, FE_11, 
      rhoS_11, aeta_11, rhoE_11,
     
      fE22, A22, alpha22, rE22, kappa22, B22, beta22, lamb22, Fn0_22, Fn1_22, 
      Fn2_22, Fn3_22, F0_22, F2_22, F3_22, FE_22, 
      rhoS_22, aeta_22, rhoE_22) = params
    
    F1_11 = 0
    F1_22 = 0
    
    filePath = 'CrNb.eam.alloy' # Nb comes first in file, followed by Cr
    try:
        os.remove(filePath)
        print(f"File '{filePath}' has been deleted.")
    except OSError as e:
        print(f"Error deleting the file: {e}")
        
    gridResolution = 2000
    cutoff = 0.6128704555450525e1
    rhoHi = 169.091775601
    
    deltaR = cutoff/gridResolution
    deltaRho = rhoHi/gridResolution
        
    rValues = linspace(0,cutoff,gridResolution)
    rhoValues = linspace(0,rhoHi,gridResolution)
    
    eamHeader = '''
Nb-Cr potential                        
Generated December 2023        
Niobium (Z = 41)  Chromium (Z = 24)                                           
    2 Nb Cr
  2000  '''+str(deltaRho)+''' 2000  0.3065885220335430E-02  0.6128704555450525E+01
  41     0.92906380000E+02         3.81         bcc
'''
    with open(filePath, 'a') as file:
        file.write(eamHeader)

        # Initialize a counter for the number of values written per line
        values_per_line = 0

        for i, rho in enumerate(rhoValues):
            # Add a space before the first value on each line
            if values_per_line == 0:
                file.write(' ')
            file.write(f'{embeddingFunction1(rho,Fn0_11,Fn1_11,Fn2_11,Fn3_11,F0_11,F1_11,F2_11,F3_11,FE_11,rhoS_11,aeta_11,rhoE_11): <40}')

            # Check if 5 values have been written and start a new line
            if values_per_line == 4 or i == len(rhoValues) - 1:
                file.write('\n')
                values_per_line = -1  # Reset the counter for the next line
            else:
                file.write(' ')
            values_per_line += 1

        for i, dist in enumerate(rValues):
            
            if values_per_line == 0:
                file.write(' ')
            file.write(f'{atomicElectronDensity11(fE11, beta11, rE11, lamb11, dist): <40}')

            
            if values_per_line == 4 or i == len(rValues) - 1:
                file.write('\n')
                values_per_line = -1  
            else:
                file.write(' ')
            values_per_line += 1

        file.write('  24     0.5199610000E+02         3.01         bcc' + '\n')
        # 24     0.5199610000E+02         3.01         FCC

        for i, rho in enumerate(rhoValues):
            if values_per_line == 0:
                file.write(' ')
            file.write(f'{embeddingFunction2(rho,Fn0_22,Fn1_22,Fn2_22,Fn3_22,F0_22,F1_22,F2_22,F3_22,FE_22,rhoS_22,aeta_22,rhoE_22): <40}')
            
            if values_per_line == 4 or i == len(rhoValues) - 1:
                file.write('\n')
                values_per_line = -1  
            else:
                file.write(' ')
            values_per_line += 1


        for i, dist in enumerate(rValues):
            if values_per_line == 0:
                file.write(' ')
            file.write(f'{atomicElectronDensity22(fE22, beta22, rE22, lamb22, dist): <40}')

            if values_per_line == 4 or i == len(rValues) - 1:
                file.write('\n')
                values_per_line = -1  
            else:
                file.write(' ')
            values_per_line += 1

        for i, dist in enumerate(rValues):
            if values_per_line == 0:
                file.write(' ')
            file.write(f'{dist*phi11(dist,A11,alpha11,rE11,kappa11,B11,beta11,lamb11): <40}')

            if values_per_line == 4 or i == len(rValues) - 1:
                file.write('\n')
                values_per_line = -1  
            else:
                file.write(' ')
            values_per_line += 1

        for i, dist in enumerate(rValues):
            if values_per_line == 0:
                file.write(' ')
            file.write(f'{dist*phi12(dist,A11,alpha11,rE11,kappa11,B11,beta11,lamb11,A22,alpha22,rE22,kappa22,B22,beta22,lamb22,fE11,fE22): <40}')

            if values_per_line == 4 or i == len(rValues) - 1:
                file.write('\n')
                values_per_line = -1  
            else:
                file.write(' ')
            values_per_line += 1

        for i, dist in enumerate(rValues):
            
            if values_per_line == 0:
                file.write(' ')
            file.write(f'{dist*phi22(dist,A22,alpha22,rE22,kappa22,B22,beta22,lamb22): <40}')

            if values_per_line == 4 or i == len(rValues) - 1:
                file.write('\n')
                values_per_line = -1  
            else:
                file.write(' ')
            values_per_line += 1

    remove_first_blank_line(filePath)
    return

def computeError(x,y,z):
    d = x-y
    r = (d)**2/z
    return r
    
def costFunction(parameters):
        
    params = parameters
    writeEAM(params)
    modelData = computeParameters()
    print(modelData)
    with open('updatedParameters.txt', 'a') as file:
        file.write(str(modelData)+ '\n')
    diff = computeError(modelData,DFTdata,y2)

    updatedParams = ', '.join(map(str, params))
    print("Updated parameters:", updatedParams)
    with open('updatedParameters.txt', 'a') as file:
        file.write(str(params)+ '\n')
    
    return diff


guessCr = array([1.8864562457808338, 0.4041382092387275, 9.81738457853451, 
                 2.4847881119845674, 0.1667403354236807, 0.6481734144115148, 
                 5.236901535343175, 0.3561747906423359, -2.534992, -0.059605, 
                 0.193065, -2.282322, -2.5400291735714244, 0.19999188841990154, 
                 -0.14845019105373414, -2.539944606877472, 20.04146299932221, 
                 0.3917500076683635, 20.041436287823505])
    
guessNb = array([3.087916143650174, 0.6134434703020438, 8.474862420653936, 
                    2.874517833905598, 0.17489861752799038, 1.0229033328371946, 
                    4.528402412356032, 0.3573468865463416, -5.542861877308831, 
                    0.655183950530759, 1.1089372840235132, -3.5847317549277773, 
                    -5.679260668751694, 1.6390902441521729, 0.22165795211659892, 
                    -5.1415277581159975, 33.787169668397134, 0.8487943002544411, 
                    33.78632095128671])

initialGuess = concatenate((guessNb,guessCr))
least_squares(costFunction,initialGuess,bounds=[-6,40],ftol=3e-16,xtol=3e-16,gtol=3e-16) 