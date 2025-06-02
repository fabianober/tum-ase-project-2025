import math
import helpers

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
crossecProp = crosssectional_properties_tee_skin(height_str=45, width_str=40, thickness_web=3, thickness_flange=3, thickness_skin=2, stringer_pitch=200)
print(f"Area: {crossecProp[0]}, Moment of Inertia: {crossecProp[1]}")

#Column Buckling formulas 
#Euler Buckling case 
def EulerBuckling(EModulus, I_y, area, length, sigma_applied, c=1):
    # r = math.sqrt(I_y/area) 
    # lmd = (c*length)/r
    lmd = lmd(I_y, area, length, c)
    sigma_crit = round(math.pi**2 * EModulus/(lmd**2))
    reserveFactor = sigma_crit/sigma_applied
    return sigma_crit, reserveFactor

#Euler Johnson with Crippling
def sigma_crip(EModulus, height_str, thickness_flange, thickness_web,sigma_yield, r = 0):
    #We assume a T-Stringer attached to the skin
    ki = 0.41   #Support factor 
    #Effective width of web 
    b12 =  height_str - thickness_flange/2 * (2-0.5*thickness_web/thickness_flange-0.2*r**2/(thickness_flange*thickness_web))
    xi = b12/thickness_web * math.sqrt(sigma_yield/ki*EModulus) #slenderness of the beam 
    #Compute alpha scaling 
    alpha = 0
    if 0.4 <= xi <= 1.095:
        alpha = 1.4-0.628*xi
    elif 1.095 < xi <=1.633:
        alpha = 0.78/xi
    elif 1.633 < xi:
        alpha = 0.69/ pow(xi,0.75)

    sigma_crippling = alpha * sigma_yield   #Compute crippling stress
    return sigma_crippling

def EulerJohnson(EModulus, I_y, area, length, sigma_yield, sigma_applied, c=1):
    r = r = math.sqrt(I_y/area) #radius of gyration 
    lmd = (c*length)/r      #lambda slenders factor
    sigma_crippel = sigma_crip()    #returns the crippling stress of the T-stringer
    sigma_cutoff = min(sigma_crippel, sigma_yield)  #Determine the inzterpolation stress
    sigma_crit = sigma_cutoff - 1/EModulus*(sigma_cutoff/(2*math.pi))**2 * lmd**2 # interpolate crictical stress
    reserveFactor = sigma_crit/sigma_applied
    return sigma_crit, reserveFactor 

def RambergOsgood():
    return False


def RambergOsgoodIt(EModulus, I_y, area, length, sigma_applied, sigma_02, sigma_u, epsilon_u, c=1, tol=0.01, max_iter=1000):
    # Radius of gyration
    r = math.sqrt(I_y / area)
    # Slenderness ratio
    lambda_ = (c * length) / r
    # Ramberg-Osgood exponent
    n = math.log(epsilon_u / 0.002) / math.log(sigma_u / sigma_02)

    # Initial guess for critical stress (Euler)
    sigma_crit_old = (math.pi**2 * EModulus) / (lambda_**2)

    for i in range(max_iter):
        # Tangent modulus based on current sigma_crit
        denominator = 1 + 0.002 * n * (EModulus / sigma_02) * ((sigma_crit_old / sigma_02) ** (n - 1))
        Et = EModulus / denominator

        # Update critical stress
        sigma_crit_new = (math.pi**2 * Et) / (lambda_**2)

        # Check for convergence
        if abs(sigma_crit_new - sigma_crit_old) < tol:
            break

        sigma_crit_old = sigma_crit_new

    reserveFactor = sigma_crit_new / sigma_applied
    return round(sigma_crit_new, 2), round(reserveFactor, 2), round(Et, 2), i + 1  # Include Et and iteration count

# Example usage of RambergOsgoodIt
res = RambergOsgoodIt(EModulus=72000, I_y=crossecProp[1], area=crossecProp[0], length=600, sigma_applied=200, sigma_02=280, sigma_u=350, epsilon_u=0.1)
# we expect arround: Et=64605, sigma_crit=218.9, reserveFactor=0.78
print(res)