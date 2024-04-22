#!/usr/bin/env python
# coding: utf-8

# 
# 

# # PLANETARY AND MOLECULAR CONSTRAINTS
# 

# This section contains the PLANETARY CONSTRAINTS of this project

# In[ ]:


import numpy as np
import pandas as pd 
import ipywidgets as widgets
import tkinter as tk

df = pd.read_csv("exoplanet.eu_catalog (jupiter<1).csv", nrows= 5000 )

df.columns = ["Planet Name", "Planet Mass (MJup)",  "Planet Radius (RJup)", "Star Name", 
              "Star Distance to Us (pc)", "Star Mass (SunM)", "Star Radius", "Star Type", "Has Atmosphere"]

#To add column that is in MEarth
df['Planet Mass (MEarth)']=df["Planet Mass (MJup)"]*[317.907]
#To add column that is in REarth
df['Planet Radius (REarth)']=df["Planet Radius (RJup)"]*[11.2089]
#To add column that is in LightYears
df['Star Distance to Us (LightYears)']=df["Star Distance to Us (pc)"]*[3.26156]
#print(df)



#makes table look nice
df.head(5000)


#def stellarkillK(x):
    #return x == 'K'

#result = map(stellarkillK, ['Stellar Photolysis'])

#print(list(result))  


def assign_stellar_kill(value):
    if (value == 'K0V'):  
        return 'Yes'
    else:
        return 'No'


df['assign_stellar_kill'] = df['Star Type'].apply(assign_stellar_kill)





quitWhole = True
while(quitWhole):
    print('Task avalible: starType, starMass, starDistance, starRadius, planetMassE, planetRadiusE, stop')
    task = input()
    if(task == 'starType'):
        quitTerm = True
        starTermArray = []
        i = 0
        #while-loop while true if not stop
        while(quitTerm):
            print('Enter Type of Star or type stop to quit')
            starType = input()
            if(starType == 'stop'):
                print('star Selection Stopped')
                quitTerm = False

            else:
                starTermArray.append(starType)
                print('Star Selected: ' + starType)
                i = i + 1 #keeps track of length

        #or = |
        for i in range(0, i):#instead of i use len 
            df = df[(df['Star Type'].str.startswith(starTermArray[i], na=False))]
        
    elif(task == 'starMass'):
 
        # minMass or maxMass enter --, should enter condional statement for only the other limit
        print('Enter min mass or - for no min')
        minMass = input()
        print('Enter max mass or - for no max')
        maxMass = input()
        if((minMass != "-")):
            df = df[(df['Star Mass (SunM)']>= float(minMass))]
        if((maxMass != "-")):
            df = df[(df['Star Mass (SunM)']<= float(maxMass))]
            

        
    elif(task == 'starDistance'):
        # minMass or maxMass enter --, should enter condional statement for only the other limit
        print('Enter min Distance or - for no min')
        minDist = input()
        print('Enter max Distance or - for no max')
        maxDist = input()
        if((minDist != "-")):
            df = df[(df['Star Distance to Us (LightYears)']>= float(minDist))]
        if((maxDist != "-")):
            df = df[(df['Star Distance to Us (LightYears)']<= float(maxDist))]

        
             
    elif(task == 'planetMassE'):
        # minMass or maxMass enter --, should enter condional statement for only the other limit
        print('Enter min planetMassE or - for no min')
        minplanetMassE = input()
        print('Enter max planetMassE or - for no max')
        maxplanetMassE = input()
        if((minplanetMassE != "-")):
            df = df[(df['Planet Mass (MEarth)']>= float(minplanetMassE))]
        if((maxplanetMassE != "-")):
            df = df[(df['Planet Mass (MEarth)']<= float(maxplanetMassE))] 
      
         
    elif(task == 'planetRadiusE'):
        # minMass or maxMass enter --, should enter condional statement for only the other limit
        print('Enter min planetRadiusE or - for no min')
        minplanetRadiusE = input()
        print('Enter max planetRadiusE or - for no max')
        maxplanetRadiusE = input()
        if((minplanetRadiusE != "-")):
            df = df[(df['Planet Radius (REarth)']>= float(minplanetRadiusE))]
        if((maxplanetRadiusE != "-")):
            df = df[(df['Planet Radius (REarth)']<= float(maxplanetRadiusE))] 
    
    
    elif(task == 'starRadius'):
        # minMass or maxMass enter --, should enter condional statement for only the other limit
        print('Enter min starRadius or - for no min')
        minplanetRadiusE = input()
        print('Enter max starRadius or - for no max')
        maxplanetRadiusE = input()
        if((minstarRadius != "-")):
            df = df[(df['Star Radius']>= float(minstarRadius))]
        if((maxstarRadius!= "-")):
            df = df[(df['Star Radius']<= float(maxstarRadius))] 
    
    else:
        quitWhole = False
df

#df["Star Distance to Us (pc)"](df["Star Distance to Us (pc)"].min())
#df.sort_values(by=

#['Star Distance to Us (pc)'] #Sorted dataframe by the distance to earth
#df[df['Star Distance to Us (pc)'] == df['Star Distance to Us (pc)'].min()]


# In[1]:


import pandas as pd
import ipywidgets as widgets
from IPython.display import Image

# Load the data
df = pd.read_csv("exoplanet.eu_catalog (jupiter<1).csv", nrows=5000)

df.columns = ["Planet Name", "Planet Mass (MJup)",  "Planet Radius (RJup)", "Star Name", 
              "Star Distance to Us (pc)", "Star Mass (SunM)", "Star Radius", "Star Type", "Has Atmosphere"]

# To add column that is in MEarth
df['Planet Mass (MEarth)'] = df["Planet Mass (MJup)"] * [317.907]
# To add column that is in REarth
df['Planet Radius (REarth)'] = df["Planet Radius (RJup)"] * [11.2089]
# To add column that is in LightYears
df['Star Distance to Us (LY)'] = df["Star Distance to Us (pc)"] * [3.26156]

# Dictionary to store the count of messages
message_count = {
    "promising": 0,
    "conservative": 0,
}

# Dropdown menu for selecting star type
Star_Type_widget = widgets.SelectMultiple(
    options=['K', 'G', 'M','F', 'WD'],
    value=['K', 'G', 'M','F', 'WD'],
    description='Star Type:',
    style={'description_width': 'initial', 'button_color': '#FF69B4'}  # Adjust the width of the description and change button color
)

# Dropdown menu for selecting star name
#star_name_widget = widgets.SelectMultiple(
    #options=list(df['Star Name'].unique()),
    #value=list(df['Star Name'].unique()),
    #description='Star Name:',
    #style={'description_width': 'initial', 'button_color': '#FF69B4'}  # Adjust the width of the description and change button color
#)

# Dropdown menu for selecting planet name
#planet_name_widget = widgets.SelectMultiple(
   # options=list(df['Planet Name'].unique()),
   # value=list(df['Planet Name'].unique()),
   # description='Planet Name:',
   # style={'description_width': 'initial', 'button_color': '#FF69B4'}  # Adjust the width of the description and change button color
#)

# Slider for selecting star mass range
star_mass_range_widget = widgets.FloatRangeSlider(
    value=[0.5, 5],
    min=0,
    max=10,
    step=0.1,
    description='Star Mass (SunM):',
    continuous_update=False,
    style={'description_width': 'initial', 'handle_color': '#FF69B4'}  # Adjust the width of the description and change handle color
)

# Slider for selecting planet mass range
planet_mass_range_widget = widgets.FloatRangeSlider(
    value=[0, 20],
    min=0,
    max=500,
    step=0.1,
    description='Planet Mass (MEarth):',
    continuous_update=False,
    style={'description_width': 'initial', 'handle_color': '#FF69B4'}  # Adjust the width of the description and change handle color
)

# Slider for selecting planet radius range
planet_radius_range_widget = widgets.FloatRangeSlider(
    value=[0, 20],
    min=0,
    max=500,
    step=0.1,
    description='Planet Radius (REarth):',
    continuous_update=False,
    style={'description_width': 'initial', 'handle_color': '#FF69B4'}  # Adjust the width of the description and change handle color
)

# Slider for selecting star distance range
star_distance_range_widget = widgets.FloatRangeSlider(
    value=[0, 2000],
    min=0,
    max=10000,
    step=10,
    description='Star Distance to Us (LY):',
    continuous_update=False,
    style={'description_width': 'initial', 'handle_color': '#FF69B4', 'handle_size': '20px'}  # Adjust the width of the description and change handle color and size
)

# Checkbox for selecting whether the star has atmosphere
has_atmosphere_widget = widgets.Checkbox(
    value=False,
    description='Has Atmosphere:',
    style={'description_width': 'initial', 'handle_color': '#FF69B4'}  # Adjust the width of the description and change handle color
)

def filter_dataframe(change):
    selected_star_types = Star_Type_widget.value
    #selected_planet_names = planet_name_widget.value
    #selected_star_names = star_name_widget.value
    min_star_mass, max_star_mass = star_mass_range_widget.value
    min_planet_mass, max_planet_mass = planet_mass_range_widget.value
    min_planet_radius, max_planet_radius = planet_radius_range_widget.value
    min_star_distance, max_star_distance = star_distance_range_widget.value
    has_atmosphere = has_atmosphere_widget.value
    
    

    
    filtered_df = df
 
    
    # Filter by selected star types
    if 'All' not in selected_star_types:
        filtered_df = filtered_df[filtered_df['Star Type'].isin(selected_star_types)]
    
    # Filter by selected planet names
    #if 'All' not in selected_planet_names:
       # filtered_df = filtered_df[filtered_df['Planet Name'].isin(selected_planet_names)]
    
    # Filter by selected star names
    #if 'All' not in selected_star_names:
       # filtered_df = filtered_df[filtered_df['Star Name'].isin(selected_star_names)]
    
    # Filter by star mass range
    filtered_df = filtered_df[(filtered_df['Star Mass (SunM)'] >= min_star_mass) & (filtered_df['Star Mass (SunM)'] <= max_star_mass)]
    
    # Filter by planet mass range
    filtered_df = filtered_df[(filtered_df['Planet Mass (MEarth)'] >= min_planet_mass) & (filtered_df['Planet Mass (MEarth)'] <= max_planet_mass)]
    
    # Filter by planet radius range
    filtered_df = filtered_df[(filtered_df['Planet Radius (REarth)'] >= min_planet_radius) & (filtered_df['Planet Radius (REarth)'] <= max_planet_radius)]
    
    # Filter by star distance range
    filtered_df = filtered_df[(filtered_df['Star Distance to Us (LY)'] >= min_star_distance) & (filtered_df['Star Distance to Us (LY)'] <= max_star_distance)]
    
    # Filter by whether the star has atmosphere
    if has_atmosphere:
        filtered_df = filtered_df[filtered_df['Has Atmosphere'] == True]
    
    #if "Earth" in filtered_df["Planet Name"].values:
        #print("You only found Earth! You chose the only path of certainty. In the search for life we need to allow for uncertainty and keep exploring. Expand your search! ")
    
    #if "Earth" in filtered_df["Planet Name"].values:
       # earth_df = filtered_df[filtered_df["Planet Name"] == "Earth"]
        #if len(earth_df) == 1:  # Check if only Earth appears
          #  print("You only found Earth! You chose the only path of certainty. In the search for life we need to allow for uncertainty and keep exploring. Expand your search!")
        #filtered_df = earth_df
   
   # if len(filtered_df) == 1 and filtered_df.iloc[0]["Planet Name"] == "Earth":
        #print("You only found Earth!")
    #else:
    #if len(filtered_df) == 1 and filtered_df.iloc[0]["Planet Name"] == "Earth":
        #print("You only found Earth! You chose the only path of certainty. In the search for life we need to allow for uncertainty and keep exploring. Apply more constraints and expand your search!")   
     # Check if only Earth appears in the filtered DataFrame
    if len(filtered_df) == 1 and filtered_df.iloc[0]["Planet Name"] == "Earth":
        print("You only found Earth! You chose the only path of certainty. In the search for life we need to allow for uncertainty and keep exploring. Apply more constraints and expand your search!")
    else:
        # Display the filtered DataFrame if Earth is not the only result
        display(filtered_df)    
        
    # Display the filtered DataFrame
    display(filtered_df)
    
    
    # Print message based on the number of rows in the filtered DataFrame
    num_rows = len(filtered_df)
    if num_rows > 40:
        print("You have {} potential possibilities for life! Allowing for more candidates will expand your chances at finding life! However, it might include more uncertainty. The level of certainty matters when doing science.".format(num_rows))
        message_count["promising"] += 1
    elif num_rows > 0:
        print("You have {} potential possibilities for life!".format(num_rows))
        message_count["promising"] += 1
    else:
        print("You are too conservative. You are taking a very safe route, it's not showing much! You have to allow more uncertainty in your search, run your constraints again.")
        message_count["conservative"] += 1

# Function to calculate comparative feedback
def comparative_feedback():
    total_messages = sum(message_count.values())
    optimistic_percentage = (message_count["promising"] / total_messages) * 100
    conservative_percentage = (message_count["conservative"] / total_messages) * 100
    
    if optimistic_percentage > conservative_percentage:
        print("You are more optimistic than other users!")
    elif optimistic_percentage < conservative_percentage:
        print("You are more conservative than other users!")
    else:
        print("Your optimism level matches that of other users!")

# Function to handle start of the game
def start_game():
    print("Welcome to your ranking of potential life in the galaxy!")
    display(Image(filename='cosmos.jpg'))  # Display image
    print("Set your criteria using the widgets below and see what you find!")
    print("IMPORTANT NOTES for the Users: LY symbolizes Light Years, pc symbolizes parsecs ")
    print("For selecting multiple options, hold the COMMAND button")
    print("1 LY = 9.46 trillion kilometers")
    print("1 parsec = 3.086×10^(13)kilometers or 30,860,000,000,000 kilometers")
    print("SunM = 1 solar mass")
    print("MEarth = 1 Earth mass")
    print("REarth = 1 radius of the Earth")
    print("Where would life be?")
    print("You decide the constraints and I will tell you the possibilities!")
    print("Let's explore the cosmos!")
          
    

# Observe changes in widgets and apply the filter
Star_Type_widget.observe(filter_dataframe, names='value')
#planet_name_widget.observe(filter_dataframe, names='value')
#star_name_widget.observe(filter_dataframe, names='value')
star_mass_range_widget.observe(filter_dataframe, names='value')
planet_mass_range_widget.observe(filter_dataframe, names='value')
planet_radius_range_widget.observe(filter_dataframe, names='value')
star_distance_range_widget.observe(filter_dataframe, names='value')
has_atmosphere_widget.observe(filter_dataframe, names='value')

# Display all widgets and start the game
start_game()
display(Star_Type_widget)
#display(planet_name_widget)
display(star_mass_range_widget)
display(planet_mass_range_widget)
display(planet_radius_range_widget)
display(star_distance_range_widget)
display(has_atmosphere_widget)
#display(star_name_widget)
# Call comparative_feedback function after certain interactions


# The next section includes the MOLECULAR CONSTRAINTS OF THIS PROJECT
# 

# In[2]:


import numpy as np
import pandas as pd 

df = pd.read_csv("RASCALL_Molecule_Identifiers_with_energies.csv", nrows= 8000)
df.columns = ["SMILES", "rdkit_SMILES", "Molecular Weight", "Molecular Formula", "IUPAC chemical name", "EPISUITE-compatible SMILES", "InChI Code", "InChI Key", 'Bond Dissociation Energies']
#print(df)

df.head(8000)

#def assign_dissociation_energy(value):
    #molecules_with_2_fg = [['C','CCCCCBr' ]]  # List of star types to assign 'Yes' to
   # if value in molecules_with_2_fg :  
       # return '377.4 ± 0 8'
   # else:
     #   return 'No'
#'CCCCCBr', 
# Applying the function to the DataFrame
#df['assign_dissociation_energy'] = df['SMILES'].apply(assign_dissociation_energy)


#def assign_dissociation_energy(value):
    # List of molecules with the same dissociation energy
    #molecules_with_same_energy = ['CCCCCBr', 'C' 'CC(CCC)C', 'CC(C)=C(C)C', 'OCCCC', 'C=C(CC)CC', 'CC(CBr)C', 'OC(C1CC1)C', 'CC(CCl)C', 'CC=CC', 'OC(CC=C)C', 'CC(CCCl)=O', 'SC(C)C(O)=O', 'CC(C=C)Cl', 'CC(C=C)CC', 'CC=C(C)Br', 'CCOCC', 'CC(C=C)C', 'ClC(C)(CCl)C', 'O=CCCCC', 'CCC(C)Br', 'CCNCC', 'NCC(O)(CC)', 'CCCCC#C', 'C/C(C)=C\\CO', 'CCCC(C)Br', 'CCCCl', 'CCCCC', 'CC(CCl)Br', 'O=C(O)COC', 'NCCCNC', 'C[N+]([O-])=O', 'O=CCCC', 'CC(C(Cl)Cl)Cl', 'NCC', 'COC(CBr)=O', 'CCC1CO1', 'NCCC', 'ClC(C)(C=C)C', 'C=C(C)CO', 'C(C)SSC(C)', 'OCC#CC', 'O=COC(C)C', 'CCSCC', 'SCCCCC', 'NC', 'SCC(CC)C', 'CC#CC#CC', 'COCCOC', 'C=C=C(C)C', 'CCOCCBr', 'CCC[N+]([O-])=O', 'OC(C)=O', 'BrC(C)(C)C', 'OC(CCC)C', 'CN(C)N=O', 'BrC(C)(CC)C', 'OC(C#CC)C', 'ClC(C)(CC)C', 'CC(CC)=O', 'C(C)NCC(O)', 'OC(C(C)C)=O', 'CC1=NCCS1', 'OC', 'CC(C)(C)C#N', 'CC1COCC1', 'C=C(C)CCO', 'BrC(C)(CCl)C', 'CC[N+]([O-])=O', 'COC(CCl)=O', 'C=CC=C(C)C', 'OC(C=C)C', 'NCC(CC)C', 'CC#CC', 'CC=CC=CC', 'OCC(C)C', 'OC(=O)C#C(C)', 'CCCI', 'CC#CCCC', 'CC(O)/C=C/C', 'CC(CCl)=O', 'OC(C)C', 'C(C(Br)CCBr)', 'C#CC(C)=O', 'ClCC(NC)=O', 'COCC=C', 'CC(C)([N+]#[C-])C', 'CC(NC)=O', 'O=S(C)C', 'OC(CC)CC', 'C=C(CCC)C', 'COC(CC)=O', 'CC(C)=O', 'NC(CCC)C', 'NCCCC', 'OCCC', 'CCCC', 'CN(N(C)C)C', 'CCCl', 'CC(Cl)Cl', 'O=COC', 'CC(COC)=O', 'CC(N(C)C)=O', 'CC(C)Cl', 'O=C(OC)OC', 'OCCOCC', 'CC(C)I', 'CC(C)C', 'COC(C)=O', 'NCCSCC', 'CCC=C', 'O=CN(C)C', 'CSCCCl', 'N#CCC(C)=C', 'OC(CC)=O', 'COC', 'CC(CC#C)C', 'O=CCC', 'ClCC=C(C)Cl', 'CN(N)C', 'CCCC#C', 'SC(CC)C', 'N#CN(C)C', 'OC(CCl)C', 'CN(C)C', 'CCCCI', 'C=C(C)C', 'CC(C=CC)C', 'CC(Br)=C', 'COCCBr', 'CCC=CC=C', 'COCCl', 'SCCCC', 'CC(C=C)C=C', 'CC(N)=O', 'SCCC(C)C', 'CCI', 'CCC', 'CCl', 'CC(CCCl)C', 'NCCNCC', 'NC(C)C', 'CC1CCC=C1', 'C=COC(C)=O', 'CC(NCC)=O', 'NC(C)(CC)C', 'CC(CC=C)C', 'N#CCOC', 'OCC', 'ClC(C(N)=O)C', 'OC(=O)C(Cl)(C)', 'C(C(F)F)', 'SC(C)C', 'CCNCC=C', 'OC(C)C(O)=O', 'CC(C#N)=C', 'CC(C)(CC)C', 'N#CC', 'OC(C(C)C)C', 'CCSC(C)C', 'CCCCCCl', 'CNCC=C', 'CC(C#C)C', 'CC=C(C)C', 'CCC#N', 'CCC#C', 'S=C(C)N(C)C', 'CC(CC)CCl', 'CSSC', 'CC=C', 'CSC#N', 'CC(Cl)Br', 'CNCCNC', 'CCCCC=C', 'CNCC(C)C', 'N#CCSC', 'NC(C)(C)C', 'CC1CCC1', 'CC(C1CC1)=C', 'CC(C1CC1)=O', 'ClC(C)(C)C', 'CI', 'CC', 'C/C=C\\CO', 'CC(CCl)Cl', 'CN=C=O', 'NCC(C)=C', 'CCSC=C', 'C(C(Cl)(F)F)', 'COC=C', 'CC1=CCCC1', 'C/C=C\\C#N', 'CC1CCCC1', 'CC(Cl)(Cl)Cl', 'C=CC(C)=O', 'CC1=CSC=C1', 'CC(CBr)Br', 'CC=CCCC', 'S=C=NCOC', 'CC(Cl)(CCl)Cl', 'OC[Si](C)(C)C', 'CCOCCCl', 'C=C=CC', 'COCOC', 'CC(CC)C', 'NCCCCC', 'OCC#CCC', 'OCCOC', 'CC(C)(C=C)C', 'O=CC(CC)C', 'O=COCCC', 'CC(CI)C', 'CCOC=C', 'CCOC(C)=O', 'CC#CCC', 'C=C(CC=C)C', 'CC(C)(CCl)C', 'CC(C(C)Cl)Cl', 'C(CC(Br)CBr)', 'OC(=O)C(CC)', 'CC(C)C#N', 'S=C=NCCC', 'OC(=O)C(C)=C', 'CN1C=CC=C1', 'OC(C#C)C', 'CCCCBr', 'OC(C#C)CC', 'ClCOCC', 'CC(CCC)=O', 'NC(C)(C)C(O)', 'C=C(C=CC)C', 'CC(C)(CBr)C', 'OCCC(C)C', 'NC(COC)=O', 'CC(C1CC1)C', 'O=COCC', 'CCBr', 'O=CNCC', 'COC(C#C)=O', 'NC(C(C)O)=O', 'ClC=C(C)C', 'CC(CC)CC', 'CC(C)Br', 'C(F)', 'CCCCCI', 'CCCCCC', 'O=CNC', 'SCC', 'OCCN(C)C', 'CCC#CCC', 'C/C=C/CC#N', 'COC(C=C)=O', 'CC1CO1', 'CCCC(C)I', 'CC(CCCl)Cl', 'CC=C(C=C)C', 'C=C=CCC', 'NCCCOC', 'S=CN(C)C', 'CC(C#C)CC', 'CCCC=C', 'O=CC=CC', 'CC([N+]([O-])=O)C', 'C=C(C=C)C', 'CC#C', 'CC1(CC1)C', 'OCCC#CC', 'CCC1CC1', 'CNC', 'CC(CBr)CCl', 'CCC=C(C)C', 'OC(C)(C)C', 'OCCCOC', 'CC1=CC=CO1', 'C[Si](C)(C)C', 'CC(C)(C)C', 'CC(C(C)C)C', 'CC=C(Br)Br', 'CC(C(C)C)=O', 'CCC1CCC1', 'C[Si](C)(C#C)C', 'CN(CC)C', 'C=C(CCl)C', 'O=S(C)(C)=O', 'O=CCC(C)C', 'CCSCC#C', 'COC(C)(C)C', 'CC(CCBr)C', 'CC(C(C)=C)C', 'O=CC', 'S=C=NC(C)C', 'CCC(CC)Br', 'NCC(C)O', 'O=[N+]([O-])OCC', 'CBr', 'SCCC', 'SC(C)(C)C', 'OCCCCC', 'CNCC#C', 'C/C=N/O', 'C=C(CC)C', 'CCC(CC)=O', 'SCC(C)C', 'OC(C)(CC)C', 'OC(COC)C', 'CC(C)(C#C)C']  # Add more molecules here if needed
    
    # Check if the given molecule is in the list
   # if value in molecules_with_same_energy:
   #     return '377.4 ± 0.8'  # Assign the same dissociation energy
  #  else:
     #   return 'No'

# Applying the function to the DataFrame
#df['assign_dissociation_energy'] = df['EPISUITE-compatible SMILES'].apply(assign_dissociation_energy)

#def assign_dissociation_energy(value):
    # List of molecules with corresponding energy
   # molecules_with_energy = ['COC(CC)=O']

    # Check if the value is in any of the lists
   # for molecules, energy in molecules_with_energy:
    #    if value in molecules:
    #        return '377.4 ± 0.8'  # Return the corresponding energy
   # return 'No'  # Return 'No' if the molecule is not found

# Applying the function to the DataFrame
#df['assign_dissociation_energy'] = df['SMILES'].apply(assign_dissociation_energy)

#def assign_dissociation_energy(value):
    # List of molecules with corresponding energy
    #molecules_with_energy = ['OCCCl', 'CC(CCl)C', 'CC(CCCl)=O', 'CC(C=C)Cl', 'ClC(C)(CCl)C', 'CCCCl', 'CC(CCl)Br', 'OC(=O)C(Cl)=C', 'CC(C(Cl)Cl)Cl', 'C(Cl)(Cl)=C(Cl)(F)', 'ClC(C)(C=C)C', 'ClC1(Cl)CC1', 'ClCC(Cl)(Cl)Cl', 'C=CCl', 'BrCCl', 'ClC(C)(CC)C', 'OCC(Cl)=C', 'BrC(C)(CCl)C', 'COC(CCl)=O', 'OCC(Cl)(Cl)Cl', 'ClCCOC=C', 'CC(CCl)=O', 'ClCC(NC)=O', 'ClCCCCl', 'OCC(Cl)Cl', 'C(Br)(Cl)(C(Br)Cl)', 'CCCl', 'C(Cl)(Cl)(F)(F)', 'CC(Cl)Cl', 'CC(C)Cl', 'ClCC1CC1', 'C(CCl)#C(CCl)', 'CSCCCl', 'ClCC(Cl)Cl', 'OCC(O)CCl', 'ClCC=C(C)Cl', 'C(Cl)(F)(F)(F)', 'OC(CCl)C', 'BrC(Cl)Cl', 'BrC(Cl)Br', 'COCCl', 'CCl', 'CC(CCCl)C', 'ClCCl', 'C(Cl)(CCCCl)', 'ClC(C(N)=O)C', 'OC(=O)C(Cl)(C)', 'ClCI', 'C=C(Cl)Cl', 'ClC(Cl)Cl', 'CCCCCCl', 'CC(CC)CCl', 'N#CSCCl', 'CC(Cl)Br', 'ClC(C)(C)C', 'C(Cl)(CCCBr)', 'CC(CCl)Cl', 'C(C(Cl)(F)F)', 'CC(Cl)(Cl)Cl', 'CC(Cl)(CCl)Cl', 'N#CC(Cl)(Cl)Cl', 'CCOCCCl', 'C(Cl)C(=O)C(Cl)', 'ClC=C(Cl)Cl', 'ClCC(Cl)=CCl', 'ClC(Cl)=C(Cl)Cl', 'CC(C)(CCl)C', 'CC(C(C)Cl)Cl', 'C(Cl)(C(Br)CBr)', 'ClCCC(N)=O', 'ClCOCC', 'OC(CCCl)', 'ClCC=CCl', 'ClCC(CCl)Cl', 'ClCCCBr', 'C(Cl)(Cl)(F)', 'ClC=C(C)C', 'ClCCCl', 'ClCC(Br)=C', 'ClCC(C=C)Cl', 'N#CCCCCl', 'ClCCBr', 'CC(CCCl)Cl', 'C(Cl)(CCCI)', 'BrC(Cl)(Cl)Cl', 'OC(=O)C(Cl)(Cl)', 'CC(CBr)CCl', 'C(Cl)(Cl)=C(F)(F)', 'ClCC(Cl)=C', 'C=C(CCl)C', 'C(Cl)(F)(F)', 'ClC(Cl)(Cl)Cl', 'C(Cl)(Cl)(Cl)(F)', 'ClCC=C', 'N#CC(Cl)(Cl)F', 'OC(CCl)(CCl)', 'OC(CCCCl)', 'C=C(Cl)Br']

    # Check if the value is in the list
   # if value in molecules_with_energy:
        #return '377.4 ± 0.8'  # Return the corresponding energy
   # else:
     #   return 'No'  # Return 'No' if the molecule is not found

# Applying the function to the DataFrame
#df['assign_dissociation_energy'] = df['SMILES'].apply(assign_dissociation_energy)



#df


# In[ ]:


import numpy as np
import pandas as pd 

df = pd.read_csv("RASCALL_Molecule_Identifiers_with_energies.csv", nrows= 20000)
df.columns = ["SMILES", "rdkit_SMILES", "Molecular Weight", "Molecular Formula", "IUPAC chemical name", "EPISUITE-compatible SMILES", "InChI Code", "InChI Key", 'Bond Dissociation Energies']
#print(df)

df.head(8000)

def assign_dissociation_energy_2(value):
    molecules_with_2_fg =  ['C', 'CCCCCBr', 'CC(CCC)C', 'CC(C)=C(C)C', 'OCCCC', 'C=C(CC)CC', 'CC(CBr)C', 'OC(C1CC1)C', 'CC(CCl)C', 'CC=CC', 'OC(CC=C)C', 'CC(CCCl)=O', 'SC(C)C(O)=O', 'CC(C=C)Cl', 'CC(C=C)CC', 'CC=C(C)Br', 'CCOCC', 'CC(C=C)C', 'ClC(C)(CCl)C', 'O=CCCCC', 'CCC(C)Br', 'CCNCC', 'NCC(O)(CC)', 'CCCCC#C', 'C/C(C)=C\\CO', 'CCCC(C)Br', 'CCCCl', 'CCCCC', 'CC(CCl)Br', 'O=C(O)COC', 'NCCCNC', 'C[N+]([O-])=O', 'O=CCCC', 'CC(C(Cl)Cl)Cl', 'NCC', 'COC(CBr)=O', 'CCC1CO1', 'NCCC', 'ClC(C)(C=C)C', 'C=C(C)CO', 'C(C)SSC(C)', 'OCC#CC', 'O=COC(C)C', 'CCSCC', 'SCCCCC', 'NC', 'SCC(CC)C', 'CC#CC#CC', 'COCCOC', 'C=C=C(C)C', 'CCOCCBr', 'CCC[N+]([O-])=O', 'OC(C)=O', 'BrC(C)(C)C', 'OC(CCC)C', 'CN(C)N=O', 'BrC(C)(CC)C', 'OC(C#CC)C', 'ClC(C)(CC)C', 'CC(CC)=O', 'C(C)NCC(O)', 'OC(C(C)C)=O', 'CC1=NCCS1', 'OC', 'CC(C)(C)C#N', 'CC1COCC1', 'C=C(C)CCO', 'BrC(C)(CCl)C', 'CC[N+]([O-])=O', 'COC(CCl)=O', 'C=CC=C(C)C', 'OC(C=C)C', 'NCC(CC)C', 'CC#CC', 'CC=CC=CC', 'OCC(C)C', 'OC(=O)C#C(C)', 'CCCI', 'CC#CCCC', 'CC(O)/C=C/C', 'CC(CCl)=O', 'OC(C)C', 'C(C(Br)CCBr)', 'C#CC(C)=O', 'ClCC(NC)=O', 'COCC=C', 'CC(C)([N+]#[C-])C', 'CC(NC)=O', 'O=S(C)C', 'OC(CC)CC', 'C=C(CCC)C', 'COC(CC)=O', 'CC(C)=O', 'NC(CCC)C', 'NCCCC', 'OCCC', 'CCCC', 'CN(N(C)C)C', 'CCCl', 'CC(Cl)Cl', 'O=COC', 'CC(COC)=O', 'CC(N(C)C)=O', 'CC(C)Cl', 'O=C(OC)OC', 'OCCOCC', 'CC(C)I', 'CC(C)C', 'COC(C)=O', 'NCCSCC', 'CCC=C', 'O=CN(C)C', 'CSCCCl', 'N#CCC(C)=C', 'OC(CC)=O', 'COC', 'CC(CC#C)C', 'O=CCC', 'ClCC=C(C)Cl', 'CN(N)C', 'CCCC#C', 'SC(CC)C', 'N#CN(C)C', 'OC(CCl)C', 'CN(C)C', 'CCCCI', 'C=C(C)C', 'CC(C=CC)C', 'CC(Br)=C', 'COCCBr', 'CCC=CC=C', 'COCCl', 'SCCCC', 'CC(C=C)C=C', 'CC(N)=O', 'SCCC(C)C', 'CCI', 'CCC', 'CCl', 'CC(CCCl)C', 'NCCNCC', 'NC(C)C', 'CC1CCC=C1', 'C=COC(C)=O', 'CC(NCC)=O', 'NC(C)(CC)C', 'CC(CC=C)C', 'N#CCOC', 'OCC', 'ClC(C(N)=O)C', 'OC(=O)C(Cl)(C)', 'C(C(F)F)', 'SC(C)C', 'CCNCC=C', 'OC(C)C(O)=O', 'CC(C#N)=C', 'CC(C)(CC)C', 'N#CC', 'OC(C(C)C)C', 'CCSC(C)C', 'CCCCCCl', 'CNCC=C', 'CC(C#C)C', 'CC=C(C)C', 'CCC#N', 'CCC#C', 'S=C(C)N(C)C', 'CC(CC)CCl', 'CSSC', 'CC=C', 'CSC#N', 'CC(Cl)Br', 'CNCCNC', 'CCCCC=C', 'CNCC(C)C', 'N#CCSC', 'NC(C)(C)C', 'CC1CCC1', 'CC(C1CC1)=C', 'CC(C1CC1)=O', 'ClC(C)(C)C', 'CI', 'CC', 'C/C=C\\CO', 'CC(CCl)Cl', 'CN=C=O', 'NCC(C)=C', 'CCSC=C', 'C(C(Cl)(F)F)', 'COC=C', 'CC1=CCCC1', 'C/C=C\\C#N', 'CC1CCCC1', 'CC(Cl)(Cl)Cl', 'C=CC(C)=O', 'CC1=CSC=C1', 'CC(CBr)Br', 'CC=CCCC', 'S=C=NCOC', 'CC(Cl)(CCl)Cl', 'OC[Si](C)(C)C', 'CCOCCCl', 'C=C=CC', 'COCOC', 'CC(CC)C', 'NCCCCC', 'OCC#CCC', 'OCCOC', 'CC(C)(C=C)C', 'O=CC(CC)C', 'O=COCCC', 'CC(CI)C', 'CCOC=C', 'CCOC(C)=O', 'CC#CCC', 'C=C(CC=C)C', 'CC(C)(CCl)C', 'CC(C(C)Cl)Cl', 'C(CC(Br)CBr)', 'OC(=O)C(CC)', 'CC(C)C#N', 'S=C=NCCC', 'OC(=O)C(C)=C', 'CN1C=CC=C1', 'OC(C#C)C', 'CCCCBr', 'OC(C#C)CC', 'ClCOCC', 'CC(CCC)=O', 'NC(C)(C)C(O)', 'C=C(C=CC)C', 'CC(C)(CBr)C', 'OCCC(C)C', 'NC(COC)=O', 'CC(C1CC1)C', 'O=COCC', 'CCBr', 'O=CNCC', 'COC(C#C)=O', 'NC(C(C)O)=O', 'ClC=C(C)C', 'CC(CC)CC', 'CC(C)Br', 'C(F)', 'CCCCCI', 'CCCCCC', 'O=CNC', 'SCC', 'OCCN(C)C', 'CCC#CCC', 'C/C=C/CC#N', 'COC(C=C)=O', 'CC1CO1', 'CCCC(C)I', 'CC(CCCl)Cl', 'CC=C(C=C)C', 'C=C=CCC', 'NCCCOC', 'S=CN(C)C', 'CC(C#C)CC', 'CCCC=C', 'O=CC=CC', 'CC([N+]([O-])=O)C', 'C=C(C=C)C', 'CC#C', 'CC1(CC1)C', 'OCCC#CC', 'CCC1CC1', 'CNC', 'CC(CBr)CCl', 'CCC=C(C)C', 'OC(C)(C)C', 'OCCCOC', 'CC1=CC=CO1', 'C[Si](C)(C)C', 'CC(C)(C)C', 'CC(C(C)C)C', 'CC=C(Br)Br', 'CC(C(C)C)=O', 'CCC1CCC1', 'C[Si](C)(C#C)C', 'CN(CC)C', 'C=C(CCl)C', 'O=S(C)(C)=O', 'O=CCC(C)C', 'CCSCC#C', 'COC(C)(C)C', 'CC(CCBr)C', 'CC(C(C)=C)C', 'O=CC', 'S=C=NC(C)C', 'CCC(CC)Br', 'NCC(C)O', 'O=[N+]([O-])OCC', 'CBr', 'SCCC', 'SC(C)(C)C', 'OCCCCC', 'CNCC#C', 'C/C=N/O', 'C=C(CC)C', 'CCC(CC)=O', 'SCC(C)C', 'OC(C)(CC)C', 'OC(COC)C', 'CC(C)(C#C)C']
  # List of star types to assign 'Yes' to
    if value in molecules_with_2_fg :  
        return '377.4 ± 0 8'
    else:
        return '-'
                                                            
    
    
def assign_dissociation_energy_3(value):
    molecules_with_3_fg =  [ 'N', 'OCCCl', 'CC(CCl)C', 'CC(CCCl)=O', 'CC(C=C)Cl', 'ClC(C)(CCl)C', 'CCCCl', 'CC(CCl)Br', 'OC(=O)C(Cl)=C', 'CC(C(Cl)Cl)Cl', 'C(Cl)(Cl)=C(Cl)(F)', 'ClC(C)(C=C)C', 'ClC1(Cl)CC1', 'ClCC(Cl)(Cl)Cl', 'C=CCl', 'BrCCl', 'ClC(C)(CC)C', 'OCC(Cl)=C', 'BrC(C)(CCl)C', 'COC(CCl)=O', 'OCC(Cl)(Cl)Cl', 'ClCCOC=C', 'CC(CCl)=O', 'ClCC(NC)=O', 'ClCCCCl', 'OCC(Cl)Cl', 'C(Br)(Cl)(C(Br)Cl)', 'CCCl', 'C(Cl)(Cl)(F)(F)', 'CC(Cl)Cl', 'CC(C)Cl', 'ClCC1CC1', 'C(CCl)#C(CCl)', 'CSCCCl', 'ClCC(Cl)Cl', 'OCC(O)CCl', 'ClCC=C(C)Cl', 'C(Cl)(F)(F)(F)', 'OC(CCl)C', 'BrC(Cl)Cl', 'BrC(Cl)Br', 'COCCl', 'CCl', 'CC(CCCl)C', 'ClCCl', 'C(Cl)(CCCCl)', 'ClC(C(N)=O)C', 'OC(=O)C(Cl)(C)', 'ClCI', 'C=C(Cl)Cl', 'ClC(Cl)Cl', 'CCCCCCl', 'CC(CC)CCl', 'N#CSCCl', 'CC(Cl)Br', 'ClC(C)(C)C', 'C(Cl)(CCCBr)', 'CC(CCl)Cl', 'C(C(Cl)(F)F)', 'CC(Cl)(Cl)Cl', 'CC(Cl)(CCl)Cl', 'N#CC(Cl)(Cl)Cl', 'CCOCCCl', 'C(Cl)C(=O)C(Cl)', 'ClC=C(Cl)Cl', 'ClCC(Cl)=CCl', 'ClC(Cl)=C(Cl)Cl', 'CC(C)(CCl)C', 'CC(C(C)Cl)Cl', 'C(Cl)(C(Br)CBr)', 'ClCCC(N)=O', 'ClCOCC', 'OC(CCCl)', 'ClCC=CCl', 'ClCC(CCl)Cl', 'ClCCCBr', 'C(Cl)(Cl)(F)', 'ClC=C(C)C', 'ClCCCl', 'ClCC(Br)=C', 'ClCC(C=C)Cl', 'N#CCCCCl', 'ClCCBr', 'CC(CCCl)Cl', 'C(Cl)(CCCI)', 'BrC(Cl)(Cl)Cl', 'OC(=O)C(Cl)(Cl)', 'CC(CBr)CCl', 'C(Cl)(Cl)=C(F)(F)', 'ClCC(Cl)=C', 'C=C(CCl)C', 'C(Cl)(F)(F)', 'ClC(Cl)(Cl)Cl', 'C(Cl)(Cl)(Cl)(F)', 'ClCC=C', 'N#CC(Cl)(Cl)F', 'OC(CCl)(CCl)', 'OC(CCCCl)', 'C=C(Cl)Br']
  # List of star types to assign 'Yes' to
    if value in molecules_with_3_fg :  
        return '350.2 ± 1.7'
    else:
        return '-'
    
    
#'CCCCCBr', 'C','CCCCCBr' 
#NC(C1Cl)(C1Br)
# Applying the function to the DataFrame
#df['assign_dissociation_energy_2'] = df['SMILES'].apply(assign_dissociation_energy_2)
#df['assign_dissociation_energy_3'] = df['SMILES'].apply(assign_dissociation_energy_3)
#df



# Applying the function to the DataFrame
df['assign_dissociation_energy_2'] = df['SMILES'].apply(assign_dissociation_energy_2)
df['assign_dissociation_energy_3'] = df['SMILES'].apply(assign_dissociation_energy_3)

# Merging the columns into a single column
df['merged_dissociation_energies'] = df['assign_dissociation_energy_2'] + ', ' + df['assign_dissociation_energy_3']

# Drop individual columns if needed
# df.drop(['assign_dissociation_energy_2', 'assign_dissociation_energy_3'], axis=1, inplace=True)

# Display the DataFrame with the merged column

df


# In[ ]:


import pandas as pd
import ipywidgets as widgets

# Read the dataset
df = pd.read_csv("RASCALL_Molecule_Identifiers_with_energies.csv", nrows=20000)
df.columns = ["SMILES", "rdkit_SMILES", "Molecular Weight", "Molecular Formula", "IUPAC chemical name", "EPISUITE-compatible SMILES", "InChI Code", "InChI Key", 'Bond Dissociation Energies']

# Define the FloatRangeSlider widget for bond dissociation energies
bond_energy_range_widget = widgets.FloatRangeSlider(
    value=[0, 500],  # Adjust the default range as needed
    min=0,
    max=500,  # Adjust the max value based on your dataset
    step=1,
    description='Bond Dissociation Energy (kJ/mol):',
    continuous_update=False,
    style={'description_width': 'initial'}  # Adjust the width of the description
)

# Function to filter DataFrame based on selected bond dissociation energy range
def filter_by_bond_energy_range(change):
    min_energy, max_energy = change['new']
    filtered_df = df[(df['Bond Dissociation Energies'] >= min_energy) & (df['Bond Dissociation Energies'] <= max_energy)]
    display(filtered_df)  # Display the filtered DataFrame

# Observe changes in the bond energy range widget and apply the filter
bond_energy_range_widget.observe(filter_by_bond_energy_range, names='value')

# Display the bond energy range widget
display(bond_energy_range_widget)


# In[ ]:


#group 2:  ['CCCCCBr', 'CC(CCC)C', 'CC(C)=C(C)C', 'OCCCC', 'C=C(CC)CC', 'CC(CBr)C', 'OC(C1CC1)C', 'CC(CCl)C', 'CC=CC', 'OC(CC=C)C', 'CC(CCCl)=O', 'SC(C)C(O)=O', 'CC(C=C)Cl', 'CC(C=C)CC', 'CC=C(C)Br', 'CCOCC', 'CC(C=C)C', 'ClC(C)(CCl)C', 'O=CCCCC', 'CCC(C)Br', 'CCNCC', 'NCC(O)(CC)', 'CCCCC#C', 'C/C(C)=C\\CO', 'CCCC(C)Br', 'CCCCl', 'CCCCC', 'CC(CCl)Br', 'O=C(O)COC', 'NCCCNC', 'C[N+]([O-])=O', 'O=CCCC', 'CC(C(Cl)Cl)Cl', 'NCC', 'COC(CBr)=O', 'CCC1CO1', 'NCCC', 'ClC(C)(C=C)C', 'C=C(C)CO', 'C(C)SSC(C)', 'OCC#CC', 'O=COC(C)C', 'CCSCC', 'SCCCCC', 'NC', 'SCC(CC)C', 'CC#CC#CC', 'COCCOC', 'C=C=C(C)C', 'CCOCCBr', 'CCC[N+]([O-])=O', 'OC(C)=O', 'BrC(C)(C)C', 'OC(CCC)C', 'CN(C)N=O', 'BrC(C)(CC)C', 'OC(C#CC)C', 'ClC(C)(CC)C', 'CC(CC)=O', 'C(C)NCC(O)', 'OC(C(C)C)=O', 'CC1=NCCS1', 'OC', 'CC(C)(C)C#N', 'CC1COCC1', 'C=C(C)CCO', 'BrC(C)(CCl)C', 'CC[N+]([O-])=O', 'COC(CCl)=O', 'C=CC=C(C)C', 'OC(C=C)C', 'NCC(CC)C', 'CC#CC', 'CC=CC=CC', 'OCC(C)C', 'OC(=O)C#C(C)', 'CCCI', 'CC#CCCC', 'CC(O)/C=C/C', 'CC(CCl)=O', 'OC(C)C', 'C(C(Br)CCBr)', 'C#CC(C)=O', 'ClCC(NC)=O', 'COCC=C', 'CC(C)([N+]#[C-])C', 'CC(NC)=O', 'O=S(C)C', 'OC(CC)CC', 'C=C(CCC)C', 'COC(CC)=O', 'CC(C)=O', 'NC(CCC)C', 'NCCCC', 'OCCC', 'CCCC', 'CN(N(C)C)C', 'CCCl', 'CC(Cl)Cl', 'O=COC', 'CC(COC)=O', 'CC(N(C)C)=O', 'CC(C)Cl', 'O=C(OC)OC', 'OCCOCC', 'CC(C)I', 'CC(C)C', 'COC(C)=O', 'NCCSCC', 'CCC=C', 'O=CN(C)C', 'CSCCCl', 'N#CCC(C)=C', 'OC(CC)=O', 'COC', 'CC(CC#C)C', 'O=CCC', 'ClCC=C(C)Cl', 'CN(N)C', 'CCCC#C', 'SC(CC)C', 'N#CN(C)C', 'OC(CCl)C', 'CN(C)C', 'CCCCI', 'C=C(C)C', 'CC(C=CC)C', 'CC(Br)=C', 'COCCBr', 'CCC=CC=C', 'COCCl', 'SCCCC', 'CC(C=C)C=C', 'CC(N)=O', 'SCCC(C)C', 'CCI', 'CCC', 'CCl', 'CC(CCCl)C', 'NCCNCC', 'NC(C)C', 'CC1CCC=C1', 'C=COC(C)=O', 'CC(NCC)=O', 'NC(C)(CC)C', 'CC(CC=C)C', 'N#CCOC', 'OCC', 'ClC(C(N)=O)C', 'OC(=O)C(Cl)(C)', 'C(C(F)F)', 'SC(C)C', 'CCNCC=C', 'OC(C)C(O)=O', 'CC(C#N)=C', 'CC(C)(CC)C', 'N#CC', 'OC(C(C)C)C', 'CCSC(C)C', 'CCCCCCl', 'CNCC=C', 'CC(C#C)C', 'CC=C(C)C', 'CCC#N', 'CCC#C', 'S=C(C)N(C)C', 'CC(CC)CCl', 'CSSC', 'CC=C', 'CSC#N', 'CC(Cl)Br', 'CNCCNC', 'CCCCC=C', 'CNCC(C)C', 'N#CCSC', 'NC(C)(C)C', 'CC1CCC1', 'CC(C1CC1)=C', 'CC(C1CC1)=O', 'ClC(C)(C)C', 'CI', 'CC', 'C/C=C\\CO', 'CC(CCl)Cl', 'CN=C=O', 'NCC(C)=C', 'CCSC=C', 'C(C(Cl)(F)F)', 'COC=C', 'CC1=CCCC1', 'C/C=C\\C#N', 'CC1CCCC1', 'CC(Cl)(Cl)Cl', 'C=CC(C)=O', 'CC1=CSC=C1', 'CC(CBr)Br', 'CC=CCCC', 'S=C=NCOC', 'CC(Cl)(CCl)Cl', 'OC[Si](C)(C)C', 'CCOCCCl', 'C=C=CC', 'COCOC', 'CC(CC)C', 'NCCCCC', 'OCC#CCC', 'OCCOC', 'CC(C)(C=C)C', 'O=CC(CC)C', 'O=COCCC', 'CC(CI)C', 'CCOC=C', 'CCOC(C)=O', 'CC#CCC', 'C=C(CC=C)C', 'CC(C)(CCl)C', 'CC(C(C)Cl)Cl', 'C(CC(Br)CBr)', 'OC(=O)C(CC)', 'CC(C)C#N', 'S=C=NCCC', 'OC(=O)C(C)=C', 'CN1C=CC=C1', 'OC(C#C)C', 'CCCCBr', 'OC(C#C)CC', 'ClCOCC', 'CC(CCC)=O', 'NC(C)(C)C(O)', 'C=C(C=CC)C', 'CC(C)(CBr)C', 'OCCC(C)C', 'NC(COC)=O', 'CC(C1CC1)C', 'O=COCC', 'CCBr', 'O=CNCC', 'COC(C#C)=O', 'NC(C(C)O)=O', 'ClC=C(C)C', 'CC(CC)CC', 'CC(C)Br', 'C(F)', 'CCCCCI', 'CCCCCC', 'O=CNC', 'SCC', 'OCCN(C)C', 'CCC#CCC', 'C/C=C/CC#N', 'COC(C=C)=O', 'CC1CO1', 'CCCC(C)I', 'CC(CCCl)Cl', 'CC=C(C=C)C', 'C=C=CCC', 'NCCCOC', 'S=CN(C)C', 'CC(C#C)CC', 'CCCC=C', 'O=CC=CC', 'CC([N+]([O-])=O)C', 'C=C(C=C)C', 'CC#C', 'CC1(CC1)C', 'OCCC#CC', 'CCC1CC1', 'CNC', 'CC(CBr)CCl', 'CCC=C(C)C', 'OC(C)(C)C', 'OCCCOC', 'CC1=CC=CO1', 'C[Si](C)(C)C', 'CC(C)(C)C', 'CC(C(C)C)C', 'CC=C(Br)Br', 'CC(C(C)C)=O', 'CCC1CCC1', 'C[Si](C)(C#C)C', 'CN(CC)C', 'C=C(CCl)C', 'O=S(C)(C)=O', 'O=CCC(C)C', 'CCSCC#C', 'COC(C)(C)C', 'CC(CCBr)C', 'CC(C(C)=C)C', 'O=CC', 'S=C=NC(C)C', 'CCC(CC)Br', 'NCC(C)O', 'O=[N+]([O-])OCC', 'CBr', 'SCCC', 'SC(C)(C)C', 'OCCCCC', 'CNCC#C', 'C/C=N/O', 'C=C(CC)C', 'CCC(CC)=O', 'SCC(C)C', 'OC(C)(CC)C', 'OC(COC)C', 'CC(C)(C#C)C']
    


# In[ ]:


#group 3: ['OCCCl', 'CC(CCl)C', 'CC(CCCl)=O', 'CC(C=C)Cl', 'ClC(C)(CCl)C', 'CCCCl', 'CC(CCl)Br', 'OC(=O)C(Cl)=C', 'CC(C(Cl)Cl)Cl', 'C(Cl)(Cl)=C(Cl)(F)', 'ClC(C)(C=C)C', 'ClC1(Cl)CC1', 'ClCC(Cl)(Cl)Cl', 'C=CCl', 'BrCCl', 'ClC(C)(CC)C', 'OCC(Cl)=C', 'BrC(C)(CCl)C', 'COC(CCl)=O', 'OCC(Cl)(Cl)Cl', 'ClCCOC=C', 'CC(CCl)=O', 'ClCC(NC)=O', 'ClCCCCl', 'OCC(Cl)Cl', 'C(Br)(Cl)(C(Br)Cl)', 'CCCl', 'C(Cl)(Cl)(F)(F)', 'CC(Cl)Cl', 'CC(C)Cl', 'ClCC1CC1', 'C(CCl)#C(CCl)', 'CSCCCl', 'ClCC(Cl)Cl', 'OCC(O)CCl', 'ClCC=C(C)Cl', 'C(Cl)(F)(F)(F)', 'OC(CCl)C', 'BrC(Cl)Cl', 'BrC(Cl)Br', 'COCCl', 'CCl', 'CC(CCCl)C', 'ClCCl', 'C(Cl)(CCCCl)', 'ClC(C(N)=O)C', 'OC(=O)C(Cl)(C)', 'ClCI', 'C=C(Cl)Cl', 'ClC(Cl)Cl', 'CCCCCCl', 'CC(CC)CCl', 'N#CSCCl', 'CC(Cl)Br', 'ClC(C)(C)C', 'C(Cl)(CCCBr)', 'CC(CCl)Cl', 'C(C(Cl)(F)F)', 'CC(Cl)(Cl)Cl', 'CC(Cl)(CCl)Cl', 'N#CC(Cl)(Cl)Cl', 'CCOCCCl', 'C(Cl)C(=O)C(Cl)', 'ClC=C(Cl)Cl', 'ClCC(Cl)=CCl', 'ClC(Cl)=C(Cl)Cl', 'CC(C)(CCl)C', 'CC(C(C)Cl)Cl', 'C(Cl)(C(Br)CBr)', 'ClCCC(N)=O', 'ClCOCC', 'OC(CCCl)', 'ClCC=CCl', 'ClCC(CCl)Cl', 'ClCCCBr', 'C(Cl)(Cl)(F)', 'ClC=C(C)C', 'ClCCCl', 'ClCC(Br)=C', 'ClCC(C=C)Cl', 'N#CCCCCl', 'ClCCBr', 'CC(CCCl)Cl', 'C(Cl)(CCCI)', 'BrC(Cl)(Cl)Cl', 'OC(=O)C(Cl)(Cl)', 'CC(CBr)CCl', 'C(Cl)(Cl)=C(F)(F)', 'ClCC(Cl)=C', 'C=C(CCl)C', 'C(Cl)(F)(F)', 'ClC(Cl)(Cl)Cl', 'C(Cl)(Cl)(Cl)(F)', 'ClCC=C', 'N#CC(Cl)(Cl)F', 'OC(CCl)(CCl)', 'OC(CCCCl)', 'C=C(Cl)Br']

