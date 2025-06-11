import math
import numpy as np 

#Define variables to check range of half waves 
N = 10
M = 20


def uniaxialF_calc(EModulus, nu, length, width, thickness, sigma_x):
    print("Uniaxial free edge")
    sigma_crit = (EModulus*pow(math.pi,2)) / (12*(1-pow(nu,2))) * pow(thickness/length,2)   #Calculate the critical stress 
    reserveFactor = sigma_crit/ (sigma_x*1.5)  #Calculate the reserve factor for this case, 1.5 to get ultimate loads for the reserve factor
    return sigma_crit, abs(reserveFactor)

def uniaxialSS_calc(EModulus, nu, length, width, thickness, sigma_x):
    print("Uniaxial simply supported")
    global N 
    global M
    sigma_crit_it = dict()
    alpha = length/width
    #Loop over the number of half waves in length direction
    for m in range(1,M):
        k_sigma = pow((m/alpha + alpha/m),2)
        sigma_e = (EModulus*pow(math.pi,2))/(12*(1-pow(nu,2))) * pow(thickness/width,2)
        sigma_crit_it.update({m:k_sigma * sigma_e}) 	    #Collect critcical values in a dictionary 
    finalM = min(sigma_crit_it, key = sigma_crit_it.get)    #Recover the number of half waves, where the minimum critical stress occurs 
    sigma_crit = sigma_crit_it[finalM]                      #Also recover the corresponding critical stress 
    reserveFactor = sigma_crit/ (sigma_x*1.5)                    #Calculate the reserve factor with it, 1.5 to get ultimate loads for the reserve factor 
    return finalM, sigma_crit, abs(reserveFactor)

def biaxialSS_calc(EModulus, nu, length, width, thickness, sigma_x, sigma_y):
    global N 
    global M
    #Control if panel needs to be flipped
    if sigma_y < sigma_x: 
        #Swap stresses 
        sigma_temp = sigma_x
        sigma_x = sigma_y
        sigma_y = sigma_temp
        #swap dimensions
        length_temp = length
        length = width
        width = length_temp

    k_sigma_it = dict()  #Create dictionary to store the critcial stresses for different n,m
    alpha = length/width
    beta = sigma_y/sigma_x
    sigma_e = (EModulus*pow(math.pi,2))/(12*(1-pow(nu,2))) * pow(thickness/width,2)         #refernence stress
    #Looping over n half waves in width direction and over m half waves in length direction 
    for n in range(1,N):        
        for m in range(1,M):
            k_sigma = pow((m**2 + n**2 * alpha**2), 2)/ (alpha**2 * (m**2 + beta * n**2 * alpha**2))  #buckling factor 
            if k_sigma > 0:
                k_sigma_it.update({(n,m):k_sigma})                                    #critical stress dictionary in dependence of m and n 
    finalN, finalM = min(k_sigma_it, key = k_sigma_it.get)    #Select the smallest critcial stress and recover n and m 
    k_sigma_min = k_sigma_it[(finalN,finalM)] 
    sigma_crit = k_sigma_min *  sigma_e                   #And then also recover the corresponding critical stress 
    reserveFactor = sigma_crit / (sigma_x*1.5)                          #Calculate the reserve factor based on this the critical stress, 1.5 to get ultimate loads for the reserve factor
    return finalN, finalM,k_sigma_min, sigma_crit, abs(reserveFactor)

def shearSS_calc(EModulus, nu, length, width, thickness, tau_xy):
    alpha = length/width
    if alpha < 1:
        k_tau = 4 + (5.34/pow(alpha,2))
    elif alpha >= 1:
        k_tau = 5.34 + (4/pow(alpha,2))
    tau_e = (EModulus*pow(math.pi,2))/(12*(1-pow(nu,2))) * pow(thickness/width,2)
    tau_crit = tau_e * k_tau
    reserveFactor = tau_crit / (tau_xy*1.5) # 1.5 to get ultimate loads for the reserve factor
    return k_tau, tau_crit, abs(reserveFactor)

def bendingSS_calc(EModulus, nu, length, width, thickness, sigma_x):
    print("bending calculated")

def combinedBiaxialShear(EModulus, nu, length, width, thickness, sigma_x, sigma_y, tau_xy):
    finalN, finalM,k_biaxial, sigma_crit, reserveFactorBi = biaxialSS_calc(EModulus=EModulus, nu=nu, length=length, width=width, thickness=thickness, sigma_x= sigma_x, sigma_y=sigma_y)
    k_shear, tau_crit, reserveFactorShear = shearSS_calc(EModulus=EModulus, nu=nu, length=length, width=width, thickness=thickness, tau_xy=tau_xy)
    combinedReserveFactor = 1/(1/reserveFactorBi + pow(1/reserveFactorShear,2))
    return k_shear, k_biaxial, abs(combinedReserveFactor)


# Function for calling other buckling formulas with row values 
def panelBuckApply(row, EModulus, nu):
    k_shear, k_biaxial, reserveFactor = combinedBiaxialShear(EModulus=EModulus, nu=nu, length=row['length'], width=row['width'], thickness=row['thickness'], sigma_x=row['sigmaXX'], sigma_y=row['sigmaYY'], tau_xy=row['sigmaXY'] )
    return k_shear, k_biaxial, reserveFactor 


#Running test on all functions 
if __name__ == '__main__':
    # Define test data
    E = 70000
    nu = 0.3
    length = 600
    width = 200
    thickness = 5
    sigma_x = -100
    sigma_y = -20
    tau_xy = 75
    # Expected results 
    # Biaxial calculation
    finalN_bi_expect = 1 
    finalM_bi_expect = 2
    sigma_crit_bi_expect = 128.02
    biReserveFactor_expect = 1.28

    #Shear calculation 
    tau_crit_expect = 228.73
    shearReserveFactor_expect = 3.05

    # Combined shear and biaxial calculation 
    combinedReserveFactor_expect = 1.13

    # Run tests
    #uniaxialF_calc()
    #uniaxialSS_calc()
    finalN_bi, finalM_bi, sigma_crit_bi, biReserveFactor = biaxialSS_calc(EModulus=E, nu=nu, length=length, width=width, thickness=thickness, sigma_x=sigma_x, sigma_y= sigma_y)
    tau_crit, shearReserveFactor = shearSS_calc(EModulus=E, nu=nu, length=length, width= width, thickness=thickness, tau_xy=tau_xy ) 
    combinedReserveFactor = combinedBiaxialShear(EModulus=E, nu=nu, length=length, width=width, thickness=thickness, sigma_x=sigma_x, sigma_y= sigma_y, tau_xy=tau_xy)
    # Output test results
    #Biaxial test
    print('The biaxial test resulted in the following values')
    print('Number of halfwaves in width direction: '+str(finalN_bi))
    print('Number of halfwaves in longitudinal direction: '+str(finalM_bi))
    print('A critical stress of: '+str(sigma_crit_bi)+' And a reserve factor of: '+str(biReserveFactor))
    print('The test comes out to be: '+ str(finalN_bi==finalN_bi_expect and finalM_bi==finalM_bi_expect and sigma_crit_bi==sigma_crit_bi_expect and biReserveFactor==biReserveFactor_expect))
    print()
    #Shear test
    print('The shear test resulted in the following values')
    print('A critical stress of: '+str(tau_crit)+' And a reserve factor of: '+str(shearReserveFactor))
    print('The test comes out to be: '+str(shearReserveFactor==shearReserveFactor_expect and tau_crit==tau_crit_expect))
    print()
    #Combiend biaxial and shear test
    print('The combined shear and biaxial test resulted in the following values')
    print('A reserve factor of: '+str(combinedReserveFactor))
    print('The test comes out to be: '+str(combinedReserveFactor==combinedReserveFactor_expect))