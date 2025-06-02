import math

def strength_calc(row ,sigma_ult):
    sigma_avg = math.sqrt(row.sigmaXX**2 + row.sigmaYY**2 - row.sigmaXX * row.sigmaYY + 3*row.sigmaXY**2)
    reserveFactor = abs(sigma_ult/(1.5*sigma_avg))
    return round(reserveFactor, 3)