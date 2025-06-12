import hm
import hm.entities as ent 
import pandas as pd 

model = hm.Model()

#Creating panel property csv 
filter = hm.FilterByCollection(ent.Element, ent.Component)
panels = hm.Collection(model, ent.Component, "id=1,2,3,4,5,6,7,8,9,10")
panelElements = hm.Collection(model, filter, panels)

#Create empty lists for the values 
elementID = []
elementComponentName = []
elementThickness = []
elementMass = []

#Loop through the elements and query them 
for element in panelElements:
    elementID.append(element.id)
    elementThickness.append(element.thickness)
    elementMass.append(element.mass)
    elementComponentName.append(element.component.name)

panelPropertiesdf = pd.DataFrame({
    'Element ID':elementID,
    'Component Name':elementComponentName,
    'thickness':elementThickness,
    'mass':elementMass
})
panelPropertiesdf = panelPropertiesdf.set_index('Element ID')
panelPropertiesdf.to_csv('panel_properties.csv')
# Clear the arrays for reuse 
elementID.clear()
elementComponentName.clear()
elementThickness.clear()
elementMass.clear()


#Query all thge masses from the elements
allElements = hm.Collection(model, ent.Element)
#Fill array values 
for element in allElements:
    elementID.append(element.id)
    elementMass.append(element.mass)

#Assign them to the dataFrame 
massDf = pd.DataFrame({
    'Element ID':elementID,
    'mass': elementMass
})
massDf = massDf.set_index('Element ID')
massDf.to_csv('element_masses.csv')


#Now lastly we need to query the beam cross sections 
# So we create the stringer property csv 
beamsections = hm.Collection(model, ent.Beamsection)
#Create empty arrays
ComponentName = []
DIM1 = []
DIM2 = []
DIM3 = []
DIM4 = []

# Iterate over all beam sections
for beamsection in beamsections:
    ComponentName.append(beamsection.id)
    DIM1.append(beamsection.beamsect_dim1)
    DIM2.append(beamsection.beamsect_dim1)
    DIM3.append(beamsection.beamsect_dim1)
    DIM4.append(beamsection.beamsect_dim1)

#Create data frame 
stringerPropertiesDf = pd.DataFrame({
    'beamsects':ComponentName,
    'beamsect_dim1':DIM1,
    'beamsect_dim2':DIM2,
    'beamsect_dim3':DIM3,
    'beamsect_dim4':DIM4,
})
stringerPropertiesDf = stringerPropertiesDf.dropna()
stringerPropertiesDf = stringerPropertiesDf.set_index('beamsects')
stringerPropertiesDf.to_csv('stringer_properties.csv')

