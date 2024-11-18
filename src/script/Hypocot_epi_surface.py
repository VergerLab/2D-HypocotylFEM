#!/usr/bin/env python
# coding: utf-8
import sys
import os

# Add the directory where function.py is located to the Python path
sys.path.append(os.path.expanduser("./io_function/"))

from lib_fun import *

mesh_file = './www/Hypocot_epi_surface.geo'
initial_scale = 0.01
generate_mesh(mesh_file, initial_scale)
mesh_path = './www/Hypocot_epi_surface.msh'
cd = CustomDomainGmsh.read_and_convert_gmsh(path = mesh_path)
# Locate Boundaries            
ymax = cd.mesh.coordinates()[:, 1].max()
ymin = cd.mesh.coordinates()[:, 1].min()

bottom_border = Boundary(f'near(y, {ymin})')
top_border = Boundary(f'near(y, {ymax})')

poiss = 0.3
TP = 0.5
Force = 2.5

# Boundary conditions
stretch = [dirichlet([0.,-Force,0.], boundary= bottom_border), # applied displ
           dirichlet([0.0,Force,0.], boundary= top_border)] 
# Young's Modulus
young_values_by_labels = {1:10000, 2:5000}
heterogeneous_young = HeterogeneousParameter(cd.cdata, young_values_by_labels)
heterogeneous_Hyperelastic_response = HyperElasticForm(young=heterogeneous_young, poisson = poiss,
                                                   source=[0., 0., 0.],
                                                   plane_stress=True)
# Set up the BVP
nl_stretch = BVP(domain=cd, vform=heterogeneous_Hyperelastic_response, bc=stretch)

tolerance = 1e-11  # Initial tolerance
nl_stretch.solve(linear_solver='superlu', absolute_tolerance=tolerance, relative_tolerance=tolerance)
xdmf_save(path='./out/Epi_Surface/Hypocot_surface_softML.xdmf', solution=nl_stretch.solution, vform=heterogeneous_Hyperelastic_response)

# Young's Modulus
young_values_by_labels = {1:5000, 2:10000}
heterogeneous_young = HeterogeneousParameter(cd.cdata, young_values_by_labels)
heterogeneous_Hyperelastic_response = HyperElasticForm(young=heterogeneous_young, poisson = poiss,
                                                   source=[0., 0., 0.],
                                                   plane_stress=True)
# Set up the BVP
nl_stretch = BVP(domain=cd, vform=heterogeneous_Hyperelastic_response, bc=stretch)

tolerance = 1e-11  # Initial tolerance
nl_stretch.solve(linear_solver='superlu', absolute_tolerance=tolerance, relative_tolerance=tolerance)

xdmf_save(path='./out/Epi_Surface/Hypocot_surface_softCW.xdmf', solution=nl_stretch.solution, vform=heterogeneous_Hyperelastic_response)

