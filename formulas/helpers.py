import math
import columnbuckling as colbuckl

def lmd(I_y, area, length, c=1):
    r = math.sqrt(I_y/area) 
    lmd = (c*length)/r
    return round(lmd, 2)

def crosssectional_properties_tee_skin_row(row):
    return colbuckl.crosssectional_properties_tee_skin(
        height_str=row['height_str'],
        width_str=row['width_str'],
        thickness_web=row['thickness_web'],
        thickness_flange=row['thickness_flange'],
        thickness_skin=row['thickness_skin'],
        stringer_pitch=row['stringer_pitch']
    )

def crosssectional_properties_hat_skin_row(row):
    return colbuckl.crosssectional_properties_hat_skin(
        DIM1=row['DIM1'],
        DIM2=row['DIM2'],
        DIM3=row['DIM3'],
        DIM4=row['DIM4'],
        thickness_skin=row['thickness_skin'],
        stringer_pitch=row['stringer_pitch']
    )

def personal_data_provider(name):
    if name == 'yannis':
        sigma_yield = 490
        EModulus = 65241.07
        nu = 0.34 
    elif name == 'fabian':
        sigma_yield = 490
        EModulus = 65420.46
        nu = 0.34
    elif name == 'daniel':
        sigma_yield = 490
        EModulus = 65669.47
        nu = 0.34
    elif name == 'felix':
        sigma_yield = 490
        EModulus = 65143.57
        nu = 0.34
    return sigma_yield, EModulus, nu



#Running test on all functions 
if __name__ == '__main__':
    print(lmd(I_y=79820.37, area=646, length=600, c=1), "and expected: 53.97")