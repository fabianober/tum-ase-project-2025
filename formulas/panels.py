import math
import numpy as np 

#Define variables to check range of half waves 
N = 10
M = 20

def uniaxialF_calc(EModulus, nu, length, width, thickness, sigma_x):
    print("Uniaxial free edge")
    sigma_crit = round((EModulus*pow(math.pi,2)) / (12*(1-pow(nu,2))) * pow(thickness/length,2),2)   #Calculate the critical stress 
    reserveFactor = round(sigma_crit/sigma_x, 2)  #Calculate the reserve factor for this case 
    return sigma_crit, reserveFactor

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
        sigma_crit_it.update({m:round(k_sigma * sigma_e, 2)}) 	    #Collect critcical values in a dictionary 
    finalM = min(sigma_crit_it, key = sigma_crit_it.get)    #Recover the number of half waves, where the minimum critical stress occurs 
    sigma_crit = sigma_crit_it[finalM]                      #Also recover the corresponding critical stress 
    reserveFactor = round(sigma_crit/sigma_x, 2)                      #Calculate the reserve factor with it 
    return finalM, sigma_crit, reserveFactor

def biaxialSS_calc(EModulus, nu, length, width, thickness, sigma_x, sigma_y):
    print("Biaxial simply supported")
    global N 
    global M
    sigma_crit_it = dict()
    alpha = length/width
    beta = sigma_y/sigma_x
    #Looping over n half waves in width direction and over m half waves in length direction 
    for n in range(1,N):        
        for m in range(1,M):
            k_sigma = pow((m**2 + n**2 * alpha**2), 2)/ (alpha**2 * (m**2 + beta * n**2 * alpha**2))  #buckling factor 
            sigma_e = (EModulus*pow(math.pi,2))/(12*(1-pow(nu,2))) * pow(thickness/width,2)         #refernence stress
            if k_sigma > 0:
                sigma_crit_it.update({(n,m):round(k_sigma * sigma_e,2)})                                    #critical stress dictionary in dependence of m and n 
    finalN, finalM = min(sigma_crit_it, key = sigma_crit_it.get)    #Select the smallest critcial stress and recover n and m 
    sigma_crit = sigma_crit_it[(finalN,finalM)]                     #And then also recover the corresponding critical stress 
    reserveFactor = round(sigma_crit/sigma_x,2)                              #Calculate the reserve factor based on this the critical stress
    return finalN, finalM, sigma_crit, reserveFactor

def shearSS_calc(EModulus, nu, length, width, thickness, tau_xy):
    alpha = length/width
    if alpha < 1:
        k_tau = 4 + (5.34/pow(alpha,2))
    elif alpha >= 1:
        k_tau = 5.34 + (4/pow(alpha,2))
    tau_e = (EModulus*pow(math.pi,2))/(12*(1-pow(nu,2))) * pow(thickness/width,2)
    tau_crit = round(tau_e * k_tau,2)
    reserveFactor = round(tau_crit/tau_xy,2)
    return tau_crit, reserveFactor

def bendingSS_calc(EModulus, nu, length, width, thickness, sigma_x):
    print("bending calculated")

#Running test on all functions 
if __name__ == 'main':
    # Define test data 
    # Run tests 
    # Output test results
