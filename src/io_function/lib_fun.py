
from bvpy.domains import CustomDomain, CustomDomainGmsh
from bvpy.utils.pre_processing import HeterogeneousParameter
from bvpy.utils.io import save
import pyvista as pv
from bvpy.utils.visu_pyvista import visualize
from bvpy.vforms import LinearElasticForm, StVenantKirchoffForm
from bvpy.boundary_conditions import dirichlet, NormalDirichlet, NormalNeumann, neumann, Boundary
from bvpy.utils.visu import plot
import fenics as fe
import numpy as np
import sys
import meshio as mio
from bvpy import BVP
from subprocess import call, DEVNULL
import pandas as pd

def xdmf_save(path, solution, vform):
    solution.rename("Displacement Vector", "")
    strain = vform.get_strain(solution)
    strain.rename("Strain", "")
    stress = vform.get_stress(solution)
    stress.rename("Stress", "")

    xdmf_file = fe.XDMFFile(fe.MPI.comm_world, path)
    xdmf_file.parameters["flush_output"] = True
    xdmf_file.parameters["functions_share_mesh"] = True
    xdmf_file.parameters["rewrite_function_mesh"] = False    
    xdmf_file.write(solution, 1)
    xdmf_file.write(strain, 1)
    xdmf_file.write(stress, 1)
    
# Func for change mesh size
def generate_mesh(mesh_file, scale):
    print("gmsh "+ mesh_file+ " -2"+ " -clscale "+ f'{scale}')
    call(["gmsh", mesh_file, "-2", "-clscale", f'{scale}'], stdout=DEVNULL, stderr=DEVNULL)