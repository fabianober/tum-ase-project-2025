import math

def crosssectional_properties_tee_skin(height_str, width_str, thickness_web, thickness_flange, thickness_skin, stringer_pitch):

    # Individual moment of inertia calculations for the skin, flange, and web of a T-stringer
    I_y_skin = (stringer_pitch * thickness_skin**3) / 12
    I_y_flange = (width_str * thickness_flange**3) / 12
    I_y_web = (thickness_web*(height_str-thickness_web)**3)/12

    # Calculate the centroid of the T-stringer
    z_skin = -thickness_skin / 2
    z_flange = thickness_flange/2
    z_web = thickness_flange + (height_str - thickness_flange) / 2

    # Calculate the area of each component
    A_skin = stringer_pitch * thickness_skin
    A_flange = width_str * thickness_flange
    A_web = thickness_web * (height_str - thickness_flange)
    A_tot = A_skin + A_flange + A_web

    z_bar = (A_skin * z_skin + A_flange * z_flange + A_web * z_web) / A_tot


    # combined moment of inertia
    contrib_skin = I_y_skin + (z_skin-z_bar)**2 * A_skin
    contrib_flange = I_y_flange + (z_flange-z_bar)**2 * A_flange
    contrib_web = I_y_web + (z_web-z_bar)**2 * A_web
    I_y_bar = contrib_skin + contrib_flange + contrib_web

    return round(A_tot, 2), round(I_y_bar, 2)

# Example usage of crosssectional_properties_tee_skin
res = crosssectional_properties_tee_skin(height_str=45, width_str=40, thickness_web=3, thickness_flange=3, thickness_skin=2, stringer_pitch=200)
print(f"Area: {res[0]}, Moment of Inertia: {res[1]}")

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