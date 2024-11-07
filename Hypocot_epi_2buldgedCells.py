#!/usr/bin/env python
# coding: utf-8

import sys
import os

# Add the directory where function.py is located to the Python path
sys.path.append(os.path.expanduser("./io_function/"))

from lib_fun import *

mesh_file = './www/2C_buldged.geo'
initial_scale = 0.08
generate_mesh(mesh_file, initial_scale)
mesh_path = './www/2C_buldged.msh'
cd = CustomDomainGmsh.read_and_convert_gmsh(path = mesh_path)

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

poiss = 0.3
TP = 0.5
Force = 0.5

# Boundary conditions
stretch = [NormalNeumann(val=-TP, boundary=inner_border),
           dirichlet([-Force,0.0,0.], boundary= left_border), # applied displ
           dirichlet([Force,0.,0.], boundary= right_border)] 

# Young's Modulus
young_values_by_labels = {1:10000, 2:5000}
heterogeneous_young = HeterogeneousParameter(cd.cdata, young_values_by_labels)
heterogeneous_Hyperelastic_response = HyperElasticForm(young=heterogeneous_young, poisson = poiss,
                                                   source=[0., 0., 0.],
                                                   plane_stress=True)
# Set up the BVP
nl_inflation = BVP(domain=cd, vform=heterogeneous_Hyperelastic_response, bc=stretch)

tolerance = 1e-11  # Initial tolerance
nl_inflation.solve(linear_solver='superlu', absolute_tolerance=tolerance, relative_tolerance=tolerance)

