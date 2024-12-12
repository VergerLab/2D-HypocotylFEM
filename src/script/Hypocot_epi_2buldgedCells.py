#!/usr/bin/env python
# coding: utf-8

# # 2D Hypocotyl longitudinal section in 2 buldged cells
# 
# 
# This notebook contains the Finite Element Modeling (FEM) analysis code used to explore stress and strain distribution in 2D mesh of a hypocotyl longitudinal section with two subdomains: the cell wall and the adhesive layer at the interface.
# 
# > **NOTE**
# > The notebook should be run within a *bvpy-dev_gmsh* environment.
# 
# Cell-cell adhesion plays a critical role in tissue integrity and controlled separation during development. Here, we investigate how adhesion mechanisms may take place under tension and turgor pressure within the "adhesive layer model".
# 
# ### Load libraries and dependencies
# 
# All bvpy functions and classes that are related to this notebook will be loaded by sourcing lib_fun.
# 
# > **NOTE**
# > The path to the "input", "output" functions (io_function) is relative to where the notebook was open.
# > In this case the notebook was open from the notebook folder.

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
# > **Note**  
# > - The `generate_mesh()` function uses **GMSH** to create `.msh` files from `.geo` files, which define the physical surface attributes.  
# > - The `read_and_convert_gmsh()` function reads the `.msh` files and converts the data into a format compatible with **bvpy**. This is done through the `CustomDomainGMSH` class.

# In[ ]:


mesh_file = '../../data/in/2C_buldged.geo'
initial_scale = 0.08
generate_mesh(mesh_file, initial_scale)
mesh_path = '../../data/in/2C_buldged.msh'
cd = CustomDomainGmsh(fname = mesh_path)


# ### Visualizing the Mesh
# 
# To ensure the mesh has been correctly loaded and interpreted, you can visualize it using the following command:
# 

# In[ ]:


plot(cd.cdata)


# ### Locate the mesh nodes related to the boundary conditions
# 
# The different tbc are the curves that coorspond to the curvature of the edge due to initial bulging of the cells  

# In[ ]:


# Locate Boundaries            
xmax = cd.mesh.coordinates()[:, 0].max()
xmin = cd.mesh.coordinates()[:, 0].min()
ymax = cd.mesh.coordinates()[:, 1].max()
ymin = cd.mesh.coordinates()[:, 1].min()

left_border = Boundary(f'near(x, {xmin})') 
right_border = Boundary(f'near(x, {xmax})') 
bottom_border = Boundary(f'near(y, {ymin}, 0.19)')
top_border = Boundary(f'near(y, {ymax}, 0.19)')
tbc1 = Boundary(f'near(x*0.17+10.4-y, 0, 0.1)')
tbc2 = Boundary(f'near(x*-0.17-0.05-y, 0, 0.1)')
tbc3 = Boundary(f'near(x*0.17-6.95-y, 0, 0.1)')
tbc4 = Boundary(f'near(x*-0.17+17.3-y, 0, 0.1)')
tbc5 = Boundary(f'near((x-20.25)*(x-20.25)/20+10.5-y, 0, 0.12)')
tbc6 = Boundary(f'near((x-20.25)*(x-20.25)/50+10.7-y, 0, 0.1)')
tbc7 = Boundary(f'near(-((x-20.25)*(x-20.25))/20-y-0.15, 0, 0.12)')
tbc8 = Boundary(f'near(-((x-20.25)*(x-20.25))/50-y-0.3, 0, 0.1)')
tbc9 = Boundary(f'near(-((x-20.25)*(x-20.25))/5-y, 0, 0.05)')

all_borders = Boundary('all')
inner_border = all_borders & ~ left_border & ~ right_border & ~ top_border & ~ bottom_border & ~ tbc1 & ~ tbc2 & ~ tbc3 & ~ tbc4 & ~tbc5 & ~tbc6 & ~tbc7 & ~tbc8 & ~tbc9


# ### Setting Boundary Conditions and Material Properties
# 
# In the following cell, we define the boundary conditions, material properties, and the governing equations for the simulation. These parameters control how the mesh responds to external forces, displacements, and internal material properties.
# 
# #### Boundary Conditions
# - **Turgor Pressure (`TP`)**: A normal force applied to the inner boundary of the mesh to simulate the internal pressure within cells.  
# - **Displacement (`ApplyDispl`)**: A prescribed displacement applied at the left and right boundaries to simulate mechanical stretching.
# 
# The boundary conditions are implemented using:
# - `NormalNeumann`: Applies the turgor pressure as a normal force.  
# - `dirichlet`: Fixes specific displacement values at defined boundaries.
# 
# #### Material Properties
# - **Poisson's Ratio (`poiss`)**: Specifies the material's lateral expansion relative to its longitudinal compression.  
# - **Young's Modulus (`young_values_by_labels`)**: Defines the stiffness of each subdomain, with values assigned based on the subdomain label.  
# 
# #### Governing Equations
# - The **HyperElasticForm** is used to model the material's nonlinear elastic behavior. This is essential for simulating realistic deformations under stress and strain.  
# - Plane stress conditions are assumed for this 2D analysis.  
# 
# #### Nonlinear Boundary Value Problem (BVP)
# The `BVP` class is used to set up and solve the system, combining the domain, variational formulation, and boundary conditions into a single framework. The `nl_inflation` object encapsulates this setup for further analysis.
# 

# In[ ]:


poiss = 0.3
TP = 0.5
ApplyDispl = 0.5

# Boundary conditions
stretch = [NormalNeumann(val=-TP, boundary=inner_border),
           dirichlet([-ApplyDispl,0.0,0.], boundary= left_border), # applied displ
           dirichlet([ApplyDispl,0.,0.], boundary= right_border)] 

# Young's Modulus
young_values_by_labels = {1:10000, 2:5000}
heterogeneous_young = HeterogeneousParameter(cd.cdata, young_values_by_labels)
elastic_potential = StVenantKirchoffPotential(young=heterogeneous_young, poisson=poiss)
heterogeneous_Hyperelastic_response = HyperElasticForm(potential_energy=elastic_potential, source=[0., 0., 0.], plane_stress=True)

# Set up the BVP
nl_inflation = BVP(domain=cd, vform=heterogeneous_Hyperelastic_response, bc=stretch)


# ### Solving the Nonlinear System
# 
# In this next step, we solve the nonlinear system representing the deformation of the mesh. The goal is to iteratively compute the equilibrium configuration of the mesh.
# 
# #### Solver Configuration:
# - **`tolerance`**: Determines the precision of the solution. Both absolute and relative tolerances are set to `1e-11` to ensure high accuracy in the convergence process.  
# - **`nl_stretch.solve()`**: This method solves the nonlinear boundary value problem using iterative algorithms.  
# - **`linear_solver='superlu'`**: Specifies the use of the SuperLU linear solver for efficient and robust performance.

# In[ ]:


tolerance = 1e-11  # Initial tolerance
nl_inflation.solve(linear_solver='superlu', absolute_tolerance=tolerance, relative_tolerance=tolerance)


# ### Saving the Simulation Results
# 
# The results of the nonlinear stretch simulation are saved to an **XDMF** file for visualization and analysis. 
# 
# **Parameters**:
#   - **`path`**: Specifies the file path for saving the results. In this case, the file is saved to `../../data/out/2BuldgedCells/two_cells_buldged.xdmf`.
#   - **`solution`**: The computed solution from the `nl_inflation` object after solving it, representing the final configuration of the mesh after deformation.
#   - **`vform`**: The variational formulation (`heterogeneous_Hyperelastic_response`) that defines the material and mechanical properties applied during the simulation.
# 
# Saving the results in this format allows for efficient post-processing and visualization in software like Paraview, enabling the analysis of stress and strain distributions across the hypocotyl epidermis model.

# In[ ]:


xdmf_save(path="../../data/out/2BuldgedCells/two_cells_buldged_01.xdmf", solution=nl_inflation.solution, vform=heterogeneous_Hyperelastic_response)

