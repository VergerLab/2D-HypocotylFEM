> [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14872201.svg)](https://doi.org/10.5281/zenodo.14872201)  v1.0.0

# 2D-Finite Element Modeling (FEM) of hypocotyl epidermal cells

Özer Erguvan<sup>1 </sup>, Adrien Heymans<sup>1 </sup>, Asal Atakhani<sup>2</sup>, Elsa Gascon<sup>3</sup>, Olivier Ali<sup>3 </sup>, Stéphane Verger<sup>1, 2</sup>

<sub>1. Umeå Plant Science Centre (UPSC), Department of Forest Genetics and Plant Physiology, Swedish University of Agricultural Sciences, 901 83 Umeå, Sweden​</sub>

<sub>2. Umeå Plant Science Centre (UPSC), Department of Plant Physiology, Umeå University, 901 87 Umeå, Sweden​</sub>

<sub>3. Laboratoire Reproduction et Développement des Plantes, Université de Lyon, ENS de Lyon, UCB Lyon, CNRS, INRAE, INRIA, F-69342, Lyon, 69364 Cedex 07, France</sub> 

## 1. About

This repository contains the Finite Element Modeling (FEM) analysis code used to explore stress and strain distribution in 2D representations of hypocotyl epidermal cells under varied mechanical conditions. Cell-cell adhesion plays a critical role in tissue integrity and controlled separation during development. Here, we investigate how adhesion mechanisms may vary under tension by simulating different subdomain Young's modulus configurations in the hypocotyl epidermis.

The FEM framework uses the [BvPy](https://gitlab.inria.fr/mosaic/bvpy), which is a python library, based on [FEniCS](https://fenicsproject.org/) and [GMSH](https://gmsh.info/).

The different study cases are:

- 2D mesh representing the surface of a hypocotyl epidermis with staggered cell files, segmented into two subdomains: the cell and the interface between cells.
- 2D mesh of a hypocotyl longitudinal section with two subdomains: the cell wall and the adhesive layer at the interface.
- 2D mesh of a hypocotyl longitudinal section with four subdomains: the Supracellular Outer Epidermal Wall (SOEW), the Outer Epidermal Edge Filling (OEEF), the Inner Walls, and the Middle Lamella (ML)

## 2. Install

You can download the content of the repository using for instance the `git` command line tool

```bash
git clone https://github.com/VergerLab/2D-HypocotylFEM.git
cd 2D-HypocotylFEM
```

#### Requirements

- Python 3.9
- FEniCS 2019.1.0
- GMSH 4.11
- Bvpy-develop
- Paraview 5.11.1
- R 4.3.1

### From mamba/conda

>[!NOTE] 
> We recommend to use [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) to create a virtual environment and run the FEM script in it ([Anaconda](https://www.anaconda.com/download) works also)
>
> For more information on how to set-up conda, please check the [conda user guide](https://conda.io/projects/conda/en/latest/user-guide/install)

```{bash}
mamba env create -f conda/hypocot_env.yaml
mamba activate hypocot_env
```

## 3. Usage and Repository content

FEM analysis of stress and strain distribution in 2D models of the hypocotyl epidermis under varied Young's modulus conditions. We tested our hypothesis in three different configurations:

- 2D mesh representing the surface of a hypocotyl epidermis with staggered cell files, segmented into two subdomains: the cell and the interface between cells.
- 2D mesh of a hypocotyl longitudinal section with two subdomains: the cell wall and the adhesive layer at the interface.
- 2D mesh of a hypocotyl longitudinal section with four subdomains: the Supracellular Outer Epidermal Wall (SOEW), the Outer Epidermal Edge Filling (OEEF), the Inner Walls, and the Middle Lamella (ML)

The analysis can be coducted with the Rscript (OECW_analysis.R) and all ".xdmf" can be visualise with [paraview](https://www.paraview.org/).

Indepenantly, all python scripts can be run through Jupyter Notebook.

![Code workflow](./data/out/img/Workflow.png)

### Surface mesh

Stress distribution on the epidermis surface if interfaces are stiffer than the cell domains
![Stress distribution on the epidermis surface if interfaces are stiffer than the cell domains](./data/out/img/surface.png)

### 2 bulged cells

Stress distribution in the edge region of a longitudinal section of epidermal cells if interface is softer than the cell wall domains
![Stress distribution in the edge region of a longitudinal section of epidermal cells if interface is softer than the cell wall domains](./data/out/img/2buldgedcells.png)

### 2 cells with a Supracellular Outer Epidermal Wall

Stress-strain distribution in the edge region of a longitudinal section of epidermal cells if:

- all cell wall subdomain have uniform mechanical properties
- the Middle Lamella is softer than the other cell wall domains
- the Outer Epidermal Edge Filling is softer than the other cell wall domains
- the Supracellular Outer Epidermal Wall is softer than the other cell wall domains
- the SOEW and OEEF are softer than the other cell wall domains

Stress distribution in the edge region of a longitudinal section of epidermal cells if Outer Epidermal Edge Filling is softer than the other cell wall domains
![Stress distribution in the edge region of a longitudinal section of epidermal cells if Outer Epidermal Edge Filling is softer than the other cell wall domains](./data/out/img/oeef.png)

## Citation

cite the preprint as soon as it is out.
Meanwhile:

>  Özer Erguvan, Adrien Heymans, Asal Atakhani, Elsa Gascon, Olivier Ali, Stéphane Verger. 2024. "Ultrastructural Characterization of Cell Adhesion in Plants"; Plant Computational Biology Workshop 2024





