import os
from paraview.simple import *

# Path to input XDMF files folder
input_folder = "C:/Users/adhs0001/Documents/GitHub/VergerLab/2D-HypocotylFEM/data/out/SOEW/"

# Define center and radius for selection of the edge
center = (6.15, 14, 0)  
radius = 6

# Iterate over XDMF files in input folder
for file_name in os.listdir(input_folder):
    if file_name.endswith(".xdmf") and "16" in file_name:
        print(file_name)
        file = os.path.join(input_folder, file_name)
        xdmf_file = Xdmf3ReaderS(FileName=file)
        xdmf_file.PointArrays = ['Displacement Vector', 'Strain', 'Stress']
        print('loading completed')

        # Clip the dataset with a sphere
        clip = Clip(Input=xdmf_file)
        clip.ClipType = 'Sphere'
        clip.ClipType.Center = center
        clip.ClipType.Radius = radius
        clip.Invert = True  # Select cells inside the sphere

        # Create WarpByVector for displacement visualization
        warpDisplacement = WarpByVector(Input=clip)
        warpDisplacement.Vectors = ['POINTS', 'Displacement Vector']

        # Create WarpByVector for stress visualization
        warpStress = WarpByVector(Input=clip)
        warpStress.Vectors = ['POINTS', 'Displacement Vector']

        # Apply Color and Opacity Mapping
        displayDisplacement = Show(warpDisplacement, GetActiveView())
        displayDisplacement.DiffuseColor = [0.0, 0.0, 0.] 

        # Setup Stress display with proper LookupTable
        displayStress = Show(warpStress, GetActiveView())
        displayStress.ColorArrayName = ['POINTS', 'Stress']
        stressLUT = GetColorTransferFunction('Stress')  # Get LookupTable
        displayStress.LookupTable = stressLUT
        displayStress.SetRepresentationType('Surface')


        # Get active view
        renderView = GetActiveViewOrCreate('RenderView')
        renderView.Background = [1, 1, 1]  # White background
        #LoadCameraConfiguration(f"{input_folder}3junctions.pvcc")
        
        # Iterate over scaling and opacity values
        for i in range(21):
            scale_factor = i * 0.05
            warpDisplacement.ScaleFactor = scale_factor
            warpStress.ScaleFactor = scale_factor
            displayStress.Opacity = scale_factor  # Gradually increase opacity
            
            # Update the view and save screenshot
            Render()
            screenshot_path = f"{input_folder}/screenshot/{file_name[:-5]}_{i:02d}.png"
            SaveScreenshot(screenshot_path, renderView, ImageResolution=[3026, 1536])
            print(f"Saved: {screenshot_path}")

        Hide(warpDisplacement, renderView)
        Hide(warpStress, renderView)

# Path to input XDMF files folder
input_cross_file = "C:/Users/adhs0001/Documents/GitHub/VergerLab/2D-HypocotylFEM/data/out/cross/xdmf/cross_turgor_05.xdmf"

xdmf_file = Xdmf3ReaderS(FileName=input_cross_file)
xdmf_file.PointArrays = ['Displacement Vector', 'Strain', 'Stress']
# Create WarpByVector for displacement visualization
warpDisplacement = WarpByVector(Input=clip)
        warpDisplacement.Vectors = ['POINTS', 'Displacement Vector']

        # Create WarpByVector for stress visualization
        warpStress = WarpByVector(Input=clip)
        warpStress.Vectors = ['POINTS', 'Displacement Vector']

        # Apply Color and Opacity Mapping
        displayDisplacement = Show(warpDisplacement, GetActiveView())
        displayDisplacement.DiffuseColor = [0.0, 0.0, 0.] 

        # Setup Stress display with proper LookupTable
        displayStress = Show(warpStress, GetActiveView())
        displayStress.ColorArrayName = ['POINTS', 'Stress']
        stressLUT = GetColorTransferFunction('Stress')  # Get LookupTable
        displayStress.LookupTable = stressLUT
        displayStress.SetRepresentationType('Surface')


        # Get active view
        renderView = GetActiveViewOrCreate('RenderView')
        renderView.Background = [1, 1, 1]  # White background
        #LoadCameraConfiguration(f"{input_folder}3junctions.pvcc")
        
        # Iterate over scaling and opacity values
        for i in range(21):
            scale_factor = i * 0.05
            warpDisplacement.ScaleFactor = scale_factor
            warpStress.ScaleFactor = scale_factor
            displayStress.Opacity = scale_factor  # Gradually increase opacity
            
            # Update the view and save screenshot
            Render()
            screenshot_path = f"{input_folder}/screenshot/{file_name[:-5]}_{i:02d}.png"
            SaveScreenshot(screenshot_path, renderView, ImageResolution=[3026, 1536])
            print(f"Saved: {screenshot_path}")

        Hide(warpDisplacement, renderView)
        Hide(warpStress, renderView)