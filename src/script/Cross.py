
import sys
import os

# Add the directory where function.py is located to the Python path
sys.path.append(os.path.expanduser("../io_function/"))

from lib_fun import *

mesh_file = '../../data/in/Orga_cross.geo'
initial_scale = 0.3
generate_mesh(mesh_file, initial_scale)
mesh_path = '../../data/in/Orga_cross.msh'
cd = CustomDomainGmsh(fname=mesh_path)

# Boundary conditions
all_borders = Boundary('all')
out_wall = Boundary(f'near((x)*(x)+(y)*(y), 52.75*52.75, 50)') & all_borders
center = Boundary(f'near((x)*(x)+(y)*(y), 7.2*7.2,2)') & all_borders
inner_border = all_borders  & ~ out_wall

TP = 0.5
#bc = [dirichlet([0., 0. , 0.], boundary = out_wall)]
bc = [NormalNeumann(val=-TP, boundary=inner_border),
     dirichlet([0., 0. , 0.], boundary = center)] #NormalNeumann(val=-TP,

elastic_potential = StVenantKirchoffPotential(young=10000, poisson=0.3)
heterogeneous_Hyperelastic_response = HyperElasticForm(potential_energy=elastic_potential, source=[0., 0., 0.],
                                                       plane_stress=True)

# The Boundary value problem
cross_turgor = BVP(domain=cd, vform=heterogeneous_Hyperelastic_response, bc=bc)

# SOLVE
tolerance = 1e-10  # Initial tolerance
cross_turgor.solve(linear_solver='mumps', krylov_solver={'absolute_tolerance':1e-13}, absolute_tolerance=tolerance, relative_tolerance=tolerance,
          preconditioner = 'none',
          maximum_iterations = 500)

# SAVE
filename = f"../../data/out/cross/xdmf/cross_turgor_15.xdmf"
xdmf_save(path=filename, solution=cross_turgor.solution, vform=heterogeneous_Hyperelastic_response)