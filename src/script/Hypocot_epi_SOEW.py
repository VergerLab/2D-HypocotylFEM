#!/usr/bin/env python
# coding: utf-8

# # 2D Hypocotyl Longitudinal Section with Supracellular Outer Epidermal Wall
# 
# This notebook contains the Finite Element Modeling (FEM) analysis code used to investigate stress and strain distribution 
# in a 2D mesh of a hypocotyl longitudinal section with four distinct subdomains:
# 
# 1. **Supracellular Outer Epidermal Wall (SOEW)**  
# 2. **Outer Epidermal Edge Filling (OEEF)**  
# 3. **Middle Lamella (ML)**  
# 4. **Inner Walls**  
# 
# > **Note**  
# > This notebook is designed to run within a *bvpy-dev_gmsh* environment.  
# 
# Cell-cell adhesion is a fundamental aspect of tissue integrity and regulated separation during development. In this analysis, we explore how adhesion mechanisms function under mechanical tension and turgor pressure, specifically considering the roles of the SOEW and OEEF subdomains.  
# 
# ### Libraries and Dependencies  
# 
# All relevant **bvpy** functions and classes required for this analysis are loaded by sourcing `lib_fun`.  
# 
# > **Important**  
# > File paths for input and output functions (`io_function`) are relative to the directory from which this notebook is executed. In this case, the notebook should be opened from the `notebook` folder.
# 

# In[ ]:


import sys
import os

# Add the directory where function.py is located to the Python path
sys.path.append(os.path.expanduser("../io_function/"))

from lib_fun import *


# ### Generating the Mesh
# 
# Once the required libraries are loaded, the next step is to generate the mesh.
# 
# - The `generate_mesh()` function uses **GMSH** to create `.msh` files from `.geo` files, which define the physical surface attributes.  
# - The `read_and_convert_gmsh()` function reads the `.msh` files and converts the data into a format compatible with **bvpy**. This is done through the `CustomDomainGMSH` class.

# In[ ]:


mesh_file = '../../data/in/2C_SOEWdomain.geo'
initial_scale = 0.15
generate_mesh(mesh_file, initial_scale)
mesh_path = '../../data/in/2C_SOEWdomain.msh'
cd = CustomDomainGmsh.read_and_convert_gmsh(path = mesh_path)


# ### Visualizing the Mesh
# 
# To ensure the mesh has been correctly loaded and interpreted, you can visualize it using the following command:
# 

# In[ ]:


plot(cd.cdata)


# ### Locate the mesh nodes related to the boundary conditions

# In[ ]:


# Boundaries condition            
xmax = cd.mesh.coordinates()[:, 0].max()
xmin = cd.mesh.coordinates()[:, 0].min()
ymax = cd.mesh.coordinates()[:, 1].max()
ymin = cd.mesh.coordinates()[:, 1].min()

left_border = Boundary(f'near(x, {xmin})') 
right_border = Boundary(f'near(x, {xmax})') 
bottom_border = Boundary(f'near(y, {ymin})') & ~ left_border & ~ right_border
top_border = Boundary(f'near(y, {ymax})')

all_borders = Boundary('all')

# The inner cells
r = 3.05
xc1 = 3.
xc2 = 9.325
xc3 = -30
xc4 = 42.325
yc1, yc2 = 9.05, 1.05
corner_in1 = Boundary(f'near((x*x -{xc1}*2*x + {xc1}*{xc1} + y*y -{yc1}*2*y + {yc1}*{yc1}), {r}*{r}, 0.61)') & all_borders
corner_in2 = Boundary(f'near((x*x -{xc2}*2*x + {xc2}*{xc2} + y*y -{yc1}*2*y + {yc1}*{yc1}), {r}*{r}, 0.61)') & all_borders
corner_in3 = Boundary(f'near((x*x -{xc3}*2*x + {xc3}*{xc3} + y*y -{yc1}*2*y + {yc1}*{yc1}), {r}*{r}, 0.61)') & all_borders
corner_in4 = Boundary(f'near((x*x -{xc4}*2*x + {xc4}*{xc4} + y*y -{yc1}*2*y + {yc1}*{yc1}), {r}*{r}, 0.61)') & all_borders
innertop_epi = Boundary(f'near(y, 12,0.1)') & all_borders & ~ Boundary(f'near(x,6.,3)') & ~ Boundary(f'near(x,44,1)') & ~ Boundary(f'near(x,-31.5,1)') & ~ left_border & ~ right_border

# if the lower intercellular space is not meshed, then: 
innerlow_epi = Boundary(f'near(y, 0,0.01)') & all_borders & ~ left_border & ~ right_border & ~ Boundary(f'near(x,6.,1)')
xc5 = 5.
xc6 = 7.325
xc7 = -32
xc8 = 44.325
corner_in5 = Boundary(f'near((x*x -{xc5}*2*x + {xc5}*{xc5} + y*y -{yc2}*2*y + {yc2}*{yc2}), 1, 0.25)') & all_borders
corner_in6 = Boundary(f'near((x*x -{xc6}*2*x + {xc6}*{xc6} + y*y -{yc2}*2*y + {yc2}*{yc2}), 1, 0.25)') & all_borders
corner_in7 = Boundary(f'near((x*x -{xc7}*2*x + {xc7}*{xc7} + y*y -{yc2}*2*y + {yc2}*{yc2}), 1, 0.25)') & all_borders
corner_in8 = Boundary(f'near((x*x -{xc8}*2*x + {xc8}*{xc8} + y*y -{yc2}*2*y + {yc2}*{yc2}), 1, 0.25)') & all_borders

innermid = Boundary(f'near(y, 5.26,4.02)') & all_borders & ~ left_border & ~ right_border

inner_border = (corner_in1 | corner_in2 | corner_in3 | corner_in4 | innertop_epi | corner_in5 | corner_in6 | corner_in7 | corner_in8 | innerlow_epi | innermid)


# ## Running Parametric Simulations
# 
# The next block of code performs a parametric simulation over a set of predefined parameters stored in a CSV file. The simulations iteratively solve the boundary value problems (BVPs) for a 2D mesh under hyperelastic material models.
# 
# #### Key Steps:
# 
# 1. **Reading Input Parameters:**
#    - Parameters are loaded from `./www/sim_paper.csv` using `pandas`. 
#    - A range of rows (from `start_index` to `stop_index`) is iterated over to extract simulation-specific parameters, such as the Young modulus of the different subdomains.
# 
# 2. **Mesh Preparation:**
#    - `change_thickness_geo()`: Adjusts the thickness of the Supracellular Outer Epidermal Wall (SOEW) in the `.geo` file.
#    - `generate_mesh()`: Regenerates the mesh based on updated parameters.
#    - `CustomDomainGmsh.read_and_convert_gmsh()`: Loads the mesh for finite element analysis.
# 
# 3. **Boundary Conditions:**
#    - Stretching is applied as boundary conditions using `dirichlet` functions, simulating imposed tensile displacement at the left and right borders.
#    - Turgor pressure is applied as boundary condition using `NormalNeumann` function, simulating turgor pressure on the inner orders.
# 
# 4. **Material Models:**
#    - Two variational formulations are defined:
#      - **Linear Elastic:** `LinearElasticForm` for simpler material response.
#      - **Hyperelastic:** `HyperElasticForm` for more complex, nonlinear material response.
#    - Heterogeneous Young's modulus values are assigned to the mesh subdomains using `HeterogeneousParameter`.
# 
# 5. **Solving the Problems:**
#    - **Linear Elastic Model:** Solved with a fixed tolerance using `ln_inflation.solve()`.
# > **INFO**
# > Solving the linear elastic model just before the Hyperelastic Model benefit the nonlinear solver with a precondition state.
#    - **Hyperelastic Model:** Solved iteratively, reducing tolerance (`tolerance /= 10`) if convergence errors occur, until a solution is found or a minimum tolerance threshold (`1e-14`) is reached.
# 
# 6. **Saving Results:**
#    - Simulation results are saved in `.xdmf` format using the `xdmf_save()` function. The filename includes the simulation number (`simu`) and scaling factor (`scale`).
# 
# > **Warning**
# >
# > Exception handling ensures that the solver attempts to converge by gradually reducing the tolerance if an error occurs.
# 
# The output files are stored in the `./out/SOEW/` directory with descriptive filenames for post-processing.
# 

# In[ ]:


poiss = 0.3
df = pd.read_csv('./www/sim_paper.csv')
start_index = 0
stop_index = 6
# Iterate over each row in the DataFrame starting from the specified index
for index, row in df.iloc[start_index:stop_index].iterrows():
    # Iterate over each column in the row
    for column in df.columns:
        # extract variables
        exec(f"{column} = row['{column}']")
    print('--------------------')
    print(row)
    print(simu)
    print('--------------------')

    # Define the different parameters
    # Boundary conditions
    change_thickness_geo(file_path = mesh_file, OCW_thickness= SOEW_thickness)
    scale = initial_scale
    generate_mesh(mesh_file, scale)
    cd = CustomDomainGmsh.read_and_convert_gmsh(path = mesh_path)
    
    stretch = [NormalNeumann(val=-TP, boundary=inner_border),
           dirichlet([-Force,0.0,0.], boundary= left_border),
           dirichlet([Force,0.,0.], boundary= right_border)]
    # Create heterogeneous Young's modulus parameter
    young_values_by_labels = {1:Cuti, 2:SOEW,3:OEEF, 4:1, 5:ML, 6:IW,7:IW}
    # 1 = Cuticle layer
    # 2 = Surpacellular outer epidermis wall
    # 3 = Outer epidermal edge filling
    # 4 = Inner intercellular space
    # 5 = Middle lamella
    # 6 = Inner wall of cell 1
    # 7 = Inner wall of cell 2

    heterogeneous_young = HeterogeneousParameter(cd.cdata, young_values_by_labels)
    heterogeneous_Hyperelastic_response = HyperElasticForm(young=heterogeneous_young, poisson = poiss,
                                                   source=[0., 0., 0.],
                                                   plane_stress=True)
    heterogeneous_Linearelastic_response = LinearElasticForm(young=heterogeneous_young, poisson = poiss,
                                                   source=[0., 0., 0.],
                                                   plane_stress=True)
    # heterogeneous LINEAR elastic boundary value problem
    ln_inflation = BVP(domain=cd, vform=heterogeneous_Linearelastic_response, bc=stretch)
    
    # Solve the linear problem
    tolerance = 1e-10  # Initial tolerance
    ln_inflation.solve(linear_solver='superlu', absolute_tolerance=tolerance, relative_tolerance=tolerance)


    # heterogeneous Hyper elastic boundary value problem
    nl_inflation = BVP(domain=cd, vform=heterogeneous_Hyperelastic_response, bc=stretch)
    # Sole the non-linear problem
    while tolerance >= 1e-14:
        print(tolerance)
        try:
            nl_inflation.solve(linear_solver="superlu", absolute_tolerance=tolerance, relative_tolerance=tolerance)
            nl_inflation, heterogeneous_Hyperelastic_response 
            # Save the solution to an XDMF file
            nl_bulging = nl_inflation.solution
            filename = f"../../data/out/SOEW/sim_{simu}_{round(scale*100)}.xdmf"
            xdmf_save(path=filename, solution=nl_bulging, vform=heterogeneous_Hyperelastic_response)
            break  # Exit the loop if successful
        except Exception as e:
            error_str = str(e)
            tolerance /= 10  # decrease tolerance on each iteration

