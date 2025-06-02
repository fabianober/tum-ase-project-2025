import math

def comb_teeStringer_skin(height_str, width_str, thickness_str, thickness_skin, stringer_pitch):

    # Individual moment of inertia calculations for the skin, flange, and web of a T-stringer
    I_y_skin = round((stringer_pitch * thickness_skin**3) / 12, 2)
    I_y_flange = round((width_str * thickness_str**3) / 12, 2)
    I_y_web = round((thickness_str*(height_str-thickness_str)**3)/12, 2)

    # combined moment of inertia (will be done)

    return I_y_skin, I_y_flange, I_y_web

# Example usage
res = comb_teeStringer_skin(height_str=45, width_str=40, thickness_str=3, thickness_skin=2, stringer_pitch=200)
print("I_y_skin:", res[0])
print("I_y_flange:", res[1])
print("I_y_web:", res[2])

def EulerBuckling(EModulus, I_y, area, length, sigma_applied, c=1):
    r = math.sqrt(I_y/area) 
    lmd = (c*length)/r
    sigma_crit = round(math.pi**2 * EModulus/(lmd**2))
    reserveFactor = sigma_crit/sigma_applied
    return sigma_crit, reserveFactor

def EulerJohnson(EModulus, I_y, area, length, sigma_applied, c=1):
    return False 
def RambergOsgood():
    return False