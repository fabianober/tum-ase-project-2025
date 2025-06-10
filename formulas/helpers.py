import math
import columnbuckling as colbuckl

def lmd(I_y, area, length, c=1):
    r = math.sqrt(I_y/area) 
    lmd = (c*length)/r
    return round(lmd, 2)

def lambda_crit(EModulus, sigma_crip, sigma_yield):
    sigma_cutoff = min(sigma_crip, sigma_yield)
    return math.sqrt(2*(math.pi**2) * EModulus / sigma_cutoff)

def r_gyr(I_y, area):
    return math.sqrt(I_y / area)

def crosssectional_properties_tee_skin_row(row):
    return colbuckl.crosssectional_properties_tee_skin(
        height_str=row['height_str'],
        width_str=row['width_str'],
        thickness_web=row['thickness_web'],
        thickness_flange=row['thickness_flange'],
        thickness_skin=row['thickness_skin'],
        stringer_pitch=row['stringer_pitch']
    )

def crosssectional_properties_hat_skin_row(row, stringer_pitch, stringer_depth):
    return colbuckl.crosssectional_properties_hat_skin(
        DIM1=row['dim1'],
        DIM2=row['dim2'],
        DIM3=row['dim3'],
        DIM4=row['dim4'],
        thickness_skin=row['thickness'],
        stringer_pitch=stringer_pitch,
        stringer_depth=stringer_depth
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

def add_component_names_to_elements(df, component_mapping_df, element_id_column='Element ID', component_name_column='Component Name'):
    """
    Add component names to a dataframe by matching element IDs with a mapping dataframe.
    
    Args:
        df: Target dataframe to add component names to
        component_mapping_df: Mapping dataframe containing element ID to component name matches
        element_id_column: Column name for element ID in both dataframes
        component_name_column: Column name for component name in the mapping dataframe
    
    Returns:
        DataFrame with added component names
    """
    # Create a dictionary for faster lookups
    component_dict = component_mapping_df.set_index(element_id_column)[component_name_column].to_dict()
    
    # Add component names using map function
    df[component_name_column] = df[element_id_column].map(component_dict)
    
    return df

#Running test on all functions 
if __name__ == '__main__':
    print(lmd(I_y=79820.37, area=646, length=600, c=1), "and expected: 53.97")