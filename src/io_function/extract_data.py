
# extract data from xdmf files
# Paraview macro edited by Adrien Heymans - Sept. 2024
# 

import os
import csv
import numpy as np
from paraview.simple import *

# Disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

###################################################################################

# Path to input XDMF files folder
input_folder = "./data/out/SOEW/"
# Path to output CSV files folder
output_folder = "./data/out/SOEW/csv/"
# Define center and radius for selection of the edge
center = (6.15, 14, 0)  
radius = 6

# Iterate over XDMF files in input folder
for file_name in os.listdir(input_folder):
    if file_name.endswith(".xdmf"):
        print(file_name)
        # Create a new 'Xdmf3ReaderS'
        file = os.path.join(input_folder, file_name)
        xdmf_file = Xdmf3ReaderS(registrationName=file_name[:-5], FileName=file)
        print('loading')
        xdmf_file.PointArrays = ['Displacement Vector', 'Strain', 'Stress']
        print('loading completed')

        # Create a new 'Clip' to clip the dataset with a sphere
        clip = Clip(Input=xdmf_file)
        clip.ClipType = 'Sphere'
        clip.ClipType.Center = center
        clip.ClipType.Radius = radius
        clip.Invert = True  # Select cells inside the sphere

        out_name = output_folder+file_name[:-5]+'_a.csv'
        # save data
        SaveData(out_name, proxy=clip, PointDataArrays=['Stress'])
        
print('cropped selection saved for SOEW simulation in csv files')

###################################################################################

# Path to input XDMF files folder
input_folder = "./data/out/Epi_Surface/"
# Path to output CSV files folder
output_folder = "./data/out/Epi_Surface/csv/"

# Iterate over XDMF files in input folder
for file_name in os.listdir(input_folder):
    if file_name.endswith(".xdmf"):
        print(file_name)
        # Create a new 'Xdmf3ReaderS'
        file = os.path.join(input_folder, file_name)
        xdmf_file = Xdmf3ReaderS(registrationName=file_name[:-5], FileName=file)
        print('loading')
        xdmf_file.PointArrays = ['Displacement Vector', 'Strain', 'Stress']
        print('loading completed')

        out_name = output_folder+file_name[:-5]+'_a.csv'
        # save data
        SaveData(out_name, proxy=xdmf_file, PointDataArrays=['Stress'])

print('Surface simulation saved in csv files')

###################################################################################

# Path to input XDMF files folder
input_folder = "./data/out/2BuldgedCells/"
# Path to output CSV files folder
output_folder = "./data/out/2BuldgedCells/csv/"

# Iterate over XDMF files in input folder
for file_name in os.listdir(input_folder):
    if file_name.endswith(".xdmf"):
        print(file_name)
        # Create a new 'Xdmf3ReaderS'
        file = os.path.join(input_folder, file_name)
        xdmf_file = Xdmf3ReaderS(registrationName=file_name[:-5], FileName=file)
        print('loading')
        xdmf_file.PointArrays = ['Displacement Vector', 'Strain', 'Stress']
        print('loading completed')

        out_name = output_folder+file_name[:-5]+'_a.csv'
        # save data
        SaveData(out_name, proxy=xdmf_file, PointDataArrays=['Stress'])

print('2 buldged cells simulation saved in csv files')