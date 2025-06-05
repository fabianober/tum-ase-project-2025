# V1 of an automatic stript transferring .strs files into the 2 escel sheets "panel" and "stringer"
# Currently there seem to be some rounding issues regarding cases of 0.0448 becoming 0.04 and not 0.05, 
# but it seems to work fine for everything else. I just thought I'd share for the moment

import pandas as pd
import re
# run "pip install openpyxl" in CMD before execution

def scientific_format(value):
    return f"{float(value):.2E}"

def extract_loadcase_section(content, subcase_label):
    pattern = rf"\$SUBCASE\s+{subcase_label}\s+LC{subcase_label}.*?\$SUBCASE|\Z"
    matches = re.findall(pattern, content, re.DOTALL)
    return matches

def parse_stringers(section, subcase_label):
    data = []
    bar_data = re.search(r"\$ELEMENT STRESS\(BAR\) \[REAL\](.*?)\$ELEMENT STRESS", section, re.DOTALL)
    if not bar_data:
        return data
    bar_lines = bar_data.group(1).splitlines()
    for line in bar_lines:
        parts = re.findall(r"[-+]?\d*\.\d+E[+-]?\d+|[-+]?\d+", line)
        if len(parts) >= 11:
            element_id = int(parts[0])
            if 40 <= element_id <= 66:
                sigma_xx = scientific_format(parts[-1])
                stringer_group = (element_id - 40) // 3 + 1
                component_name = f"stringer{stringer_group}"
                data.append([element_id, component_name, sigma_xx, f"Subcase {subcase_label} (LC{subcase_label})"])
    return data

def extract_loadcase_sections(content):
    # Extract each $SUBCASE section fully
    return re.findall(r"(\$SUBCASE\s+\d+\s+LC\d+.*?)(?=\$SUBCASE|\Z)", content, re.DOTALL)

def parse_panels(section, subcase_label):
    data = []
    plate_data = re.search(r"\$ELEMENT STRESS\(PLATE\) \[REAL\](.*?)(?=\$SUBCASE|\Z)", section, re.DOTALL)
    if not plate_data:
        return data

    plate_lines = plate_data.group(1).splitlines()

    for line in plate_lines:
        parts = re.findall(r"[-+]?\d*\.\d+E[+-]?\d+|[-+]?\d+", line)
        if len(parts) == 8:
            element_id = int(parts[0])
            if 1 <= element_id <= 30:
                component_name = f"panel{((element_id - 1) // 3) + 1}"

                xx1 = float(parts[2])
                xx2 = float(parts[3])
                yy1 = float(parts[4])
                yy2 = float(parts[5])
                xy1 = float(parts[6])
                xy2 = float(parts[7])

                xx = scientific_format((xx1 + xx2) / 2)
                yy = scientific_format((yy1 + yy2) / 2)
                xy = scientific_format((xy1 + xy2) / 2)

                data.append([
                    element_id, component_name, xx, yy, "0.00E+00", xy, "0.00E+00", "0.00E+00",
                    f"Subcase {subcase_label} (LC{subcase_label})"
                ])
    return data



# MAIN
with open("ASE_Project2025_SuperPanel_FR.strs", "r") as file:
    content = file.read()

stringer_all = []
panel_all = []

sections = extract_loadcase_sections(content)
for section in sections:
    match = re.search(r"\$SUBCASE\s+(\d+)", section)
    if match:
        subcase = int(match.group(1))
        stringer_all.extend(parse_stringers(section, subcase))
        panel_all.extend(parse_panels(section, subcase))

# To lowercase filenames
df_stringer = pd.DataFrame(stringer_all, columns=["Element ID", "Component Name", "sigmaXX", "Load Case"])
df_panel = pd.DataFrame(panel_all, columns=[
    "Element ID", "Component Name", "sigmaXX", "sigmaYY", "sigmaZZ", "sigmaXY", "sigmaXZ", "sigmaYZ", "Load Case"
])

df_stringer.to_excel("stringer.xlsx", index=False)
df_panel.to_excel("panel.xlsx", index=False)