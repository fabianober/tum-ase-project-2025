import math
import helpers as hp

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

def crosssectional_properties_hat_skin(DIM1, DIM2, DIM3, DIM4, thickness, thickness_skin, stringer_pitch):
    """
    DIM1: height of the hat section
    DIM2: bottom flange width (left)
    DIM3: top horizontal part width
    DIM4: bottom flange width (right)
    thickness: thickness of the hat section (assumed uniform)
    thickness_skin: thickness of skin
    stringer_pitch: width of skin area per stringer (used in skin calc)
    """

    # Area of each part
    A_skin = stringer_pitch * thickness_skin
    A_top = DIM3 * thickness
    A_left_web = thickness * (DIM1 - thickness)
    A_right_web = thickness * (DIM1 - thickness)
    A_bottom_left = DIM2 * thickness
    A_bottom_right = DIM4 * thickness
    A_tot = A_skin + A_top + A_left_web + A_right_web + A_bottom_left + A_bottom_right

    # z-coordinates (from bottom)
    z_skin = -thickness_skin / 2
    z_bottom = thickness / 2
    z_web = (DIM1 - thickness) / 2 + thickness
    z_top = DIM1 - thickness / 2

    # Centroid (z_bar)
    z_bar = (
        A_skin * z_skin +
        A_top * z_top +
        A_left_web * z_web +
        A_right_web * z_web +
        A_bottom_left * z_bottom +
        A_bottom_right * z_bottom
    ) / A_tot

    # Moment of inertia about y-axis for each component
    I_y_skin = (stringer_pitch * thickness_skin**3) / 12
    I_y_top = (thickness * DIM3**3) / 12
    I_y_webs = 2 * (thickness * (DIM1 - thickness)**3) / 12
    I_y_bottoms = (thickness * DIM2**3) / 12 + (thickness * DIM4**3) / 12

    # Parallel axis theorem
    contrib_skin = I_y_skin + A_skin * (z_skin - z_bar)**2
    contrib_top = I_y_top + A_top * (z_top - z_bar)**2
    contrib_webs = I_y_webs + A_left_web * (z_web - z_bar)**2 + A_right_web * (z_web - z_bar)**2
    contrib_bottoms = (
        (thickness * DIM2**3) / 12 + A_bottom_left * (z_bottom - z_bar)**2 +
        (thickness * DIM4**3) / 12 + A_bottom_right * (z_bottom - z_bar)**2
    )

    I_yy = contrib_skin + contrib_top + contrib_webs + contrib_bottoms

    return round(A_tot, 2), round(I_yy, 2)


#Column Buckling formulas 
#Euler Buckling case 
def EulerBuckling(EModulus, I_y, area, length, sigma_applied, c=1):
    lmd = hp.lmd(I_y, area, length, c)
    sigma_crit = round(math.pi**2 * EModulus/(lmd**2))
    reserveFactor = sigma_crit/sigma_applied
    return sigma_crit, reserveFactor


#Euler Johnson with Crippling
def sigma_crip(EModulus, DIM1, DIM2, DIM3, sigma_yield, r = 0):
    #We have a HAT-Stringer attached to the skin
    ki = 3.6   #Support factor for relevant parts of stringer
    #Effective width of crippling-affected parts of the HAT-stringer 
    b1 = DIM1 - DIM2/2*(1 - 0.2*(r**2/DIM2**2))
    b2 = DIM3 - DIM2*(1 - 0.2*(r**2/DIM2**2))
    #slenderness of the crippling-affected parts of the HAT-stringer
    x1 = b1/DIM2 * math.sqrt(sigma_yield/(ki*EModulus))
    x2 = b2/DIM2 * math.sqrt(sigma_yield/(ki*EModulus))
    #Compute the scaling factor alpha 
    alpha = 0
    if 0.4 <= xi <= 1.095:
        alpha = 1.4-0.628*xi
    elif 1.095 < xi <=1.633:
        alpha = 0.78/xi
    elif 1.633 < xi:
        alpha = 0.69/ pow(xi,0.75)
    sigma_crippling = alpha * sigma_yield   #Compute crippling stress
    return sigma_crippling

def EulerJohnson(EModulus, I_y, area, length, height_str, thickness_flange, thickness_web, radius, sigma_yield, sigma_applied, c=1):
    lmd = hp.lmd(I_y, area, length, c)
    sigma_crippel = sigma_crip(EModulus, height_str, thickness_flange, thickness_web,sigma_yield, radius)    #returns the crippling stress of the T-stringer
    sigma_cutoff = min(sigma_crippel, sigma_yield)  #Determine the inzterpolation stress
    sigma_crit = sigma_cutoff - 1/EModulus*(sigma_cutoff/(2*math.pi))**2 * lmd**2 # interpolate crictical stress
    reserveFactor = sigma_crit/sigma_applied
    return round(sigma_crit,2), round(reserveFactor,2) 


#Ramberg Osgood
def RambergOsgoodIt(EModulus, I_y, area, length, sigma_applied, sigma_02, sigma_u, epsilon_u, c=1, tol=0.001):
    lmd = hp.lmd(I_y, area, length, c)
    n = math.log(epsilon_u / 0.002) / math.log(sigma_u / sigma_02) # exponent

    # Start from the applied stress (initial)
    sigma_crit = sigma_applied
    step = 0.2  # Initial step size
    direction = -1  # 1 for increasing, -1 for decreasing (This is arbitrary tbh)

    while True:
        # Compute tangent modulus at current stress
        denom = 1 + 0.002 * n * (EModulus / sigma_02) * ((sigma_crit / sigma_02) ** (n - 1))
        Et = EModulus / denom

        # Compute critical stress from updated tangent modulus
        sigma_new = (math.pi**2 * Et) / (lmd**2)

        diff = sigma_new - sigma_crit

        # For overshot, reverse direction *-1 and reduce step size by 50%
        if direction * diff < 0:
            step *= 0.5
            direction *= -1

        # Update sigma_crit
        sigma_crit = sigma_crit + direction * step * abs(diff)

        # Break out, if diff between sigma_new and sigma_crit is smaller than tol=0.01
        if abs(diff) < tol:
            break

    reserveFactor = sigma_crit / sigma_applied
    return round(sigma_crit, 2), round(reserveFactor, 2) #return(critical stress, reserve factor)

#Test cases for the formula 
if __name__ == '__main__':
    # Example usage of crosssectional_properties_tee_skin
    crossecProp = crosssectional_properties_tee_skin(height_str=45, width_str=40, thickness_web=3, thickness_flange=3, thickness_skin=2, stringer_pitch=200)
    print(f"Area: {crossecProp[0]}, Moment of Inertia: {crossecProp[1]}")

    # Example usage of RambergOsgoodIt
    res = RambergOsgoodIt(EModulus=72000, I_y=crossecProp[1], area=crossecProp[0], length=600, sigma_applied=200, sigma_02=280, sigma_u=350, epsilon_u=0.1)
    # we expect arround: Et=64605, sigma_crit=218.9, reserveFactor=0.78
    print(res)

    #Example for Euler Crippling 
    sigma_crit, reserveFactor = EulerJohnson(EModulus=72000, I_y = 79820.4, area=646, length=600, height_str=45, thickness_flange=3, thickness_web=3, radius = 2, sigma_yield=280, sigma_applied=200)
    sigma_crit_expect = 131.65
    reserveFactor_expect = 0.66
    print('The resulting crictical stress is: '+str(sigma_crit)+' And the corresponding reserve Factor: '+ str(reserveFactor))
    print('The test status is thus: '+str(sigma_crit==sigma_crit_expect and reserveFactor==reserveFactor_expect))