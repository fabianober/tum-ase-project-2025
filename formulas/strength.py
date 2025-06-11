import math

def panelStrength_calc(row, sigma_ult):
    sigma_avg = math.sqrt(row['sigmaXX']**2 + row['sigmaYY']**2 - row['sigmaXX'] * row['sigmaYY'] + 3*row['sigmaXY']**2)
    reserveFactor = abs(sigma_ult/(1.5*sigma_avg))
    return reserveFactor

def stringerStrength_calc(row, sigma_ult):
    reserveFactor = abs(sigma_ult/(1.5*row['sigmaXX']))
    return reserveFactor