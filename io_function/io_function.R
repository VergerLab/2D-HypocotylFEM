

# IO functions for analysis of .xdmf files converted to .csv with the extract_data.py script

# Works for 2D mesh 
# Adrien Heymans, September 2024
# 


library(tidyverse)
library(data.table)

# Function to construct the 2x2 tensor matrix
df_to_tensor <- function(df_row) {
  sigma_xx <- df_row[5]
  sigma_xy <- df_row[2]
  sigma_yx <- df_row[4]
  sigma_yy <- df_row[1]
  
  # Since the tensor is symmetric, we average the off-diagonal components
  sigma_xy <- (sigma_xy + sigma_yx) / 2
  
  # Construct the 2x2 tensor matrix
  tensor_matrix <- matrix(c(
    sigma_xx, sigma_xy,
    sigma_xy, sigma_yy
  ), nrow = 2, ncol = 2, byrow = TRUE)
  
  return(tensor_matrix)
}

# Function to compute amplitude (Frobenius norm divided by sqrt of dimension)
amplitude <- function(tensor, dim =2) {
  return(norm(tensor, type = "F") / sqrt(dim))
}

# Function to compute amplitude on horizontal axis
amplitude_y <- function(tensor) {
  return(sqrt(tensor[1,1]^2+tensor[2,1]^2))
}
amplitude_x <- function(tensor) {
  return(sqrt(tensor[2,2]^2+tensor[2,1]^2))
}

# Function to compute the spherical (isotropic) part of the tensor
sph <- function(tensor, dim = 2) {
  di = sum(unlist(diag(tensor)))
  sph_matrix <- matrix(c(
    di[1], 0,
    0, di[1]
  ), nrow = 2, ncol = 2, byrow = TRUE)
  return( sph_matrix/ dim)
}

# Function to compute the deviatoric (non-isotropic) part of the tensor
dev <- function(tensor) {
  return(tensor - sph(tensor))
}

# Function to compute anisotropy (ratio of deviatoric to spherical amplitude)
anisotropy <- function(tensor) {
  return(amplitude(dev(tensor)) / amplitude(sph(tensor)))
}

# Function to compute von Mises stress/strain
compute_von_mises <- function(df_row) {
  sigma_xx <- df_row[5]
  sigma_xy <- df_row[2]
  sigma_yx <- df_row[4]
  sigma_yy <- df_row[1]
  
  # Compute von Mises component
  sigma_vm <- sqrt(
    sigma_xx^2 + sigma_yy^2 - sigma_xx * sigma_yy + 3 * sigma_xy^2
  )
  return(sigma_vm)
}

# Function to gather data
stress_strain <- function(data){
  
  data_stress = data%>%select(starts_with("Stress"))
  data_strain = data%>%select(starts_with("Strain"))

  data$VonMises_stress <- apply(data_stress, 1, function(row) {
    compute_von_mises(row)
  })
  data$VonMises_strain <- apply(data_strain, 1, function(row) {
    compute_von_mises(row)
  })

  data$stress_magnitude <- apply(data_stress, 1, function(row) {
    tensor <- df_to_tensor(row)
    amplitude(tensor)
  })
  data$strain_magnitude <- apply(data_strain, 1, function(row) {
    tensor <- df_to_tensor(row)
    amplitude(tensor)
  })
  
  data$stress_x <- apply(data_stress, 1, function(row) {
    tensor <- df_to_tensor(row)
    amplitude_x(tensor)
  })
  data$strain_x <- apply(data_strain, 1, function(row) {
    tensor <- df_to_tensor(row)
    amplitude_x(tensor)
  })
  
  data$stress_y <- apply(data_stress, 1, function(row) {
    tensor <- df_to_tensor(row)
    amplitude_y(tensor)
  })
  data$strain_y <- apply(data_strain, 1, function(row) {
    tensor <- df_to_tensor(row)
    amplitude_y(tensor)
  })
  
  # Apply the principal stress calculation to each row
  for (i in 1:nrow(data_stress)) {
    tensor <- df_to_tensor(data_stress[i,])
    principal <- principal_stress(tensor)
    
    data$principal_stress_1[i] <- principal$stresses[1]
    data$principal_stress_2[i] <- principal$stresses[2]
    data$angle_stress[i] <- principal$angle
  }
  # Apply the principal stress calculation to each row
  for (i in 1:nrow(data_strain)) {
    tensor <- df_to_tensor(data_strain[i,])
    principal <- principal_stress(tensor)
    
    data$principal_strain_1[i] <- principal$stresses[1]
    data$principal_strain_2[i] <- principal$stresses[2]
    data$angle_strain[i] <- principal$angle
  }
  

  return(data)
}

# Function to compute principal stresses and their orientation
principal_stress <- function(tensor) {
  eig <- eigen(matrix(unlist(tensor), ncol = 2))
  
  # Principal stresses (eigenvalues)
  principal_stresses <- eig$values
  
  # Principal directions (eigenvectors)
  principal_directions <- eig$vectors
  
  angle <- atan2(principal_directions[2, 1], principal_directions[1, 1])
  
  return(list(stresses = principal_stresses, angle = angle))
}

# Function to write .obj file containing the axes of the anisotropy glyphs
glyph_to_obj <- function(tmp, path ){
  p = tibble(x = c(tmp$x1,tmp$x2), y = c(tmp$y1,tmp$y2), z = 0, id = 1:(2*nrow(tmp)))
  obj <- '# Simple Wavefront file\n'
  point_data <- paste0('v ',p$x,' ',p$y,' ',p$z, '\n', collapse = "")
  obj <- paste0(obj, point_data)
  line = tibble(id1 =1:nrow(tmp), id2 = (nrow(tmp)+1):(2*nrow(tmp)) )
  l <- paste0('l ',line$id1,' ',line$id2, '\n', collapse = "")
  obj <- paste0(obj, l)
  cat(obj, file= paste0(path))
}

