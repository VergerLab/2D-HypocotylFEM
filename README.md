# 2D-Finite Element Modeling (FEM) of hypocotyl epidermal cells

Özer Erguvan, Adrien Heymans, Asal Atakhani, Elsa Gascone, Olivier Ali, Stéphane Verger.

This repository contains codes to stretch and apply pressure on 2D meshes representing hypocotyl epidermal cells using the Finite Element Method (FEM). The FEM framework uses the [bvpy](https://gitlab.inria.fr/mosaic/bvpy) library, and the mesh generation is performed with [GMSH](https://gmsh.info/).

## Repository content

FEM analysis of stress and strain distribution in 2D models of the hypocotyl epidermis under varied Young's modulus conditions:

- 2D mesh representing the surface of a hypocotyl epidermis with staggered cell files, segmented into two subdomains: the cell and the interface between cells.
- 2D mesh of a hypocotyl longitudinal section with two subdomains: the cell wall and the adhesive layer at the interface.
- 2D mesh of a hypocotyl longitudinal section with four subdomains: the Supracellular Outer Epidermal Wall (SOEW), the Outer Epidermal Edge Filling (OEEF), the Inner Walls, and the Middle Lamella (ML)

All python scripts can be run through a Jupyter Notebook.
The analysis of the results is done with a Rscript (OECW_analysis.R) and all ".xdmf" can be visualise with [paraview](https://www.paraview.org/)  

## Requirements

To run the FEM simulation, you will need the following:

- BvPy

```{bash}
git clone -b dev_gmsh https://gitlab.inria.fr/mosaic/bvpy.git
cd bvpy
mamba env create -f conda/env_OS.yaml -n bvpy-dev_gmsh
mamba activate bvpy-dev_gmsh
python setup.py develop --prefix=$CONDA_PREFIX
```
the recommended way is to use [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) to create a virtual environment and run the FEM script in it.


- [GMSH](https://gmsh.info/#Download)

### Surface mesh

Stress distribution on the epidermis surface if interfaces are stiffer than the cell domains
![Stress distribution on the epidermis surface if interfaces are stiffer than the cell domains](./img/surface.png)

### 2 bulged cells

Stress distribution in the edge region of a longitudinal section of epidermal cells if interface is softer than the cell wall domains
![Stress distribution in the edge region of a longitudinal section of epidermal cells if interface is softer than the cell wall domains](./img/2buldgedcells.png)

### 2 cells with a Supracellular Outer Epidermal Wall

Stress-strain distribution in the edge region of a longitudinal section of epidermal cells if:

- all cell wall subdomain have uniform mechanical properties
- the Middle Lamella is softer than the other cell wall domains
- the Outer Epidermal Edge Filling is softer than the other cell wall domains
- the Supracellular Outer Epidermal Wall is softer than the other cell wall domains
- the SOEW and OEEF are softer than the other cell wall domains

Stress distribution in the edge region of a longitudinal section of epidermal cells if Outer Epidermal Edge Filling is softer than the other cell wall domains
![Stress distribution in the edge region of a longitudinal section of epidermal cells if Outer Epidermal Edge Filling is softer than the other cell wall domains](./img/oeef.png)

## Citation

cite the preprint as soon as it is out.
Meanwhile:

>  Özer Erguvan, Adrien Heymans, Asal Atakhani, Elsa Gascone, Olivier Ali, Stéphane Verger. 2024. "Ultrastructural Characterization of Cell Adhesion in Plants"; Plant Computational Biology Workshop 2024





