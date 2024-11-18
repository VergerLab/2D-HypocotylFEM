#!/usr/bin/env python
# coding: utf-8

import sys
import os

# Add the directory where function.py is located to the Python path
sys.path.append(os.path.expanduser("./io_function/"))

from lib_fun import *

mesh_file = './www/2C_SOEWdomain.geo'
initial_scale = 0.16
generate_mesh(mesh_file, initial_scale)
mesh_path = './www/2C_SOEWdomain.msh'
cd = CustomDomainGmsh.read_and_convert_gmsh(path = mesh_path)

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
            filename = f"./out/SOEW/sim_{simu}_{round(scale*100)}.xdmf"
            xdmf_save(path=filename, solution=nl_bulging, vform=heterogeneous_Hyperelastic_response)
            break  # Exit the loop if successful
        except Exception as e:
            error_str = str(e)
            tolerance /= 10  # decrease tolerance on each iteration
