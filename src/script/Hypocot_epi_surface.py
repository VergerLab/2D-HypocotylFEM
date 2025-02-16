#!/usr/bin/env python
# coding: utf-8

# # 2D Surface of a hypocotyl epidermis with staggered cell files
# 
# This notebook contains the Finite Element Modeling (FEM) analysis code used to investigate stress and strain distribution 
# in a 2D mesh representing the surface of a hypocotyl epidermis with staggered cell files, segmented into two subdomains: the cell and the interface between cells.
# 
# > **Note**  
# > This notebook is designed to run within a *bvpy-dev_gmsh* environment.  
# 
# Cell-cell adhesion is a fundamental aspect of tissue integrity. In this analysis, we explore how adhesion mechanisms function under mechanical tension, specifically considering the roles of the cell interface mechanical properties relatively to the cell wall domain.  
# 
# ### Libraries and Dependencies  
# 
# All relevant **bvpy** functions and classes required for this analysis are loaded by sourcing `lib_fun`.  
# 
# > **Important**  
# > File paths for input and output functions (`io_function`) are relative to the directory from which this notebook is executed. In this case, the notebook should be opened from the `notebook` folder.

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
# > - The `CustomDomainGMSH` class reads the `.msh` files and converts the data into a format compatible with **bvpy**.

# In[ ]:


mesh_file = '../../data/in/Hypocot_epi_surface.geo'
initial_scale = 0.01
generate_mesh(mesh_file, initial_scale)
mesh_path = '../../data/in/Hypocot_epi_surface.msh'
cd = CustomDomainGmsh(fname = mesh_path)


# ### Visualizing the Mesh
# 
# To ensure the mesh has been correctly loaded and interpreted, you can visualize it using the following command:

# In[ ]:


plot(cd.cdata)


# ### Locate the mesh nodes related to the boundary conditions
# 

# In[ ]:


# Locate Boundaries            
ymax = cd.mesh.coordinates()[:, 1].max()
ymin = cd.mesh.coordinates()[:, 1].min()

bottom_border = Boundary(f'near(y, {ymin})')
top_border = Boundary(f'near(y, {ymax})')


# ### Setting Boundary Conditions and Material Properties
# 
# In the following cell, we define the boundary conditions, material properties, and the governing equations for the simulation. These parameters control how the mesh responds to external forces, displacements, and internal material properties.
# 
# #### Boundary Conditions
# - **Turgor Pressure (`TP`)**: A normal force applied to the inner boundary of the mesh to simulate the internal pressure within cells.  
# - **Displacement (`ApplyDispl`)**: A prescribed displacement applied at the top and bottom boundaries to simulate mechanical stretching.
# 
# The boundary conditions are implemented using:
# - `NormalNeumann`: Applies the turgor pressure as a normal force.  
# - `dirichlet`: Fixes specific displacement values at defined boundaries.
# 
# #### Material Properties
# - **Poisson's Ratio (`poiss`)**: Specifies the material's lateral expansion relative to its longitudinal compression.  
# - **Young's Modulus (`young_values_by_labels`)**: Defines the stiffness of each subdomain, with values assigned based on the subdomain label:
#     - **1** is the cell wall domain
#     - **2** is the cell interface domain
# 
# #### Governing Equations
# - The **HyperElasticForm** is used to model the material's nonlinear elastic behavior. This is essential for simulating realistic large deformations.
# - Plane stress conditions are assumed for this 2D analysis.  
# 
# #### Nonlinear Boundary Value Problem (BVP)
# The `BVP` class is used to set up and solve the system, combining the domain, variational formulation, and boundary conditions into a single framework. The `nl_stretch` object encapsulates this setup.

# In[ ]:


poiss = 0.3
TP = 0.5
ApplyDispl = 2.5

# Boundary conditions
stretch = [dirichlet([0.,-ApplyDispl,0.], boundary= bottom_border), # applied displ
           dirichlet([0.0,ApplyDispl,0.], boundary= top_border)] 

# Young's Modulus
young_values_by_labels = {1:10000, 2:5000}
heterogeneous_young = HeterogeneousParameter(cd.cdata, young_values_by_labels)
elastic_potential = StVenantKirchoffPotential(young=heterogeneous_young, poisson=poiss)
heterogeneous_Hyperelastic_response = HyperElasticForm(potential_energy=elastic_potential, source=[0., 0., 0.],
                                                       plane_stress=True)

# Set up the BVP
nl_stretch = BVP(domain=cd, vform=heterogeneous_Hyperelastic_response, bc=stretch)


# ### Solving the Nonlinear System
# 
# In this next step, we solve the nonlinear system representing the deformation of the mesh. The goal is to iteratively compute the equilibrium configuration of the mesh.
# 
# #### Solver Configuration:
# - **`tolerance`**: Determines the precision of the solution. Both absolute and relative tolerances are set to `1e-11` to ensure high accuracy in the convergence process.  
# - **`nl_stretch.solve()`**: This method solves the nonlinear boundary value problem using iterative algorithms.  
# - **`linear_solver='superlu'`**: Specifies the use of the SuperLU linear solver for efficient and robust performance.  
# 

# In[ ]:


tolerance = 1e-11  # Initial tolerance
nl_stretch.solve(linear_solver='superlu', absolute_tolerance=tolerance, relative_tolerance=tolerance)


# ### Saving the Simulation Results
# 
# The results of the nonlinear stretch simulation are saved to an **XDMF** file for visualization and analysis. 
# 
# **Parameters**:
#   - **`path`**: Specifies the file path for saving the results. In this case, the file is saved to `../../data/out/Epi_Surface/Hypocot_surface_softML.xdmf`.
#   - **`solution`**: The computed solution from the `nl_stretch` object after solving it, representing the final configuration of the mesh after deformation.
#   - **`vform`**: The variational formulation (here, `heterogeneous_Hyperelastic_response`) that defines the material and mechanical properties applied during the simulation.
# 
# Saving the results in this format allows for efficient post-processing and visualization in software like Paraview, enabling the analysis of stress and strain distributions across the hypocotyl epidermis model.
# 

# In[ ]:


xdmf_save(path='../../data/out/Epi_Surface/Hypocot_surface_softML.xdmf', solution=nl_stretch.solution, vform=heterogeneous_Hyperelastic_response)


# ### Changing the material Properties
# - The **Young's Modulus** of each subdomain is inversed compared to the previous set-up.
#     - **1** is the cell wall domain
#     - **2** is the cell interface domain
#     
# ### Seting-up the new Boundary Value Problem (BVP)

# In[ ]:


# Young's Modulus
young_values_by_labels = {1:5000, 2:10000}
heterogeneous_young = HeterogeneousParameter(cd.cdata, young_values_by_labels)
elastic_potential = StVenantKirchoffPotential(young=heterogeneous_young, poisson=poiss)
heterogeneous_Hyperelastic_response = HyperElasticForm(potential_energy=elastic_potential, source=[0., 0., 0.],
                                                       plane_stress=True)

# Set up the BVP
nl_stretch = BVP(domain=cd, vform=heterogeneous_Hyperelastic_response, bc=stretch)


# ### Solving the Nonlinear System
# 

# In[ ]:


tolerance = 1e-11  # Initial tolerance
nl_stretch.solve(linear_solver='superlu', absolute_tolerance=tolerance, relative_tolerance=tolerance)


# ### Saving the Simulation Results

# In[ ]:


xdmf_save(path='../../data/out/Epi_Surface/Hypocot_surface_softCW.xdmf', solution=nl_stretch.solution, vform=heterogeneous_Hyperelastic_response)

