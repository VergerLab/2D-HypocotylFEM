

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


concatenate_res <- function(data){
  cn = colnames(data)
  cn[which(cn== "X.1" | cn == ":1" | cn == "Points.1")]="y_coor"
  cn[which(cn== "X.0" | cn == ":0" | cn == "Points.0" )]="x_coor"
  cn = str_replace_all(cn, "[: ]", ".")
  colnames(data) = cn
  data%>%transmute(y = y_coor,x = x_coor,
                            displ_y = Displacement.Vector.1,
                            displ_x = Displacement.Vector.0,
                            displ_magnitude = sqrt(Displacement.Vector.0^2+Displacement.Vector.1^2),
                            VonMises_stress = VonMises_stress,
                            VonMises_strain = VonMises_strain,
                            stress_magnitude = stress_magnitude,
                            strain_magnitude = strain_magnitude,
                            stress_x = stress_x,
                            strain_x = strain_x,
                            stress_y = stress_y,
                            strain_y = strain_y,
                            pcstress1 = principal_stress_1,
                            pcstress2 = principal_stress_2,
                            angle_stress = angle_stress,
                            pcstrain1 = principal_strain_1,
                            pcstrain2 = principal_strain_2,
                            angle_strain = angle_strain)%>%
    mutate(pcstress_main = ifelse(abs(pcstress2)>abs(pcstress1), abs(pcstress2), abs(pcstress1)),
           pcstress_secd = ifelse(abs(pcstress2)>abs(pcstress1), abs(pcstress1), abs(pcstress2)),
           angle_stress_mod = ifelse(abs(pcstress2)>abs(pcstress1), angle_stress+pi/2, angle_stress),
           pcstrain_main = ifelse(abs(pcstrain2)>abs(pcstrain1), abs(pcstrain2), abs(pcstrain1)),
           pcstrain_secd = ifelse(abs(pcstrain2)>abs(pcstrain1), abs(pcstrain1), abs(pcstrain2)),
           angle_strain_mod = ifelse(abs(pcstrain2)>abs(pcstrain1), angle_strain+pi/2, angle_strain),
           ani_stress = (abs(pcstress_main)-abs(pcstress_secd))/(abs(pcstress_secd)+abs(pcstress_main)),
           ani_strain = (abs(pcstrain_main)-abs(pcstrain_secd))/(abs(pcstrain_secd)+abs(pcstrain_main)))
  
}


# Function to create hexagon vertices
hexagon <- function(center_x, center_y, size = 1) {
  angles <- seq(0, 2 * pi, length.out = 7)  
  tibble(
    x = center_x + size * cos(angles),
    y = center_y + size * sin(angles)
  )
}

# Function to generate hexagonal grid
hex_grid <- function(rows, cols, size = 1) {
  hex_width <- sqrt(3) * size
  hex_height <- 2 * size
  grid <- expand.grid(
    row = 1:rows,
    col = 1:cols
  ) %>%
    mutate(
      x = col *sqrt(3)* hex_width + ifelse(row %% 2 == 0, (sqrt(3)*hex_width)/2, 0),
      y = row * hex_width/2
    )
  
  # Create hexagon polygons
  hex_data <- grid %>%
    rowwise() %>%
    mutate(hex = list(hexagon(x, y, size))) %>%
    unnest(hex, names_sep = "_") %>%  # Prevents name conflict
    mutate(id = paste(row, col, sep = "-")) %>%
    ungroup()
  
  return(hex_data)
}


# Define the function to snap a single point to the closest circle
snap_point <- function(row, radii) {
  x <- row["x"]
  y <- row["y"]
  
  dist <- sqrt(x^2 + y^2)  # Compute distance from origin
  closest_r <- radii[which.min(abs(radii - dist))]  # Find closest radius
  angle <- atan2(y, x)  # Compute angle
  
  # Compute new coordinates on the closest circle
  snapped_x <- closest_r * cos(angle)
  snapped_y <- closest_r * sin(angle)
  
  return(c(snapped_x, snapped_y, closest_r, angle))
}

# Function to apply snapping transformation
snap_to_circle <- function(data, radii) {
  snapped_coords <- t(apply(data[, c("x", "y")], 1, snap_point, radii = radii))
  colnames(snapped_coords) <- c("snapped_x", "snapped_y", "closest_r", "angle")
  
  # Combine with original data
  data <- cbind(data, snapped_coords)
  return(data)
}

create_cross_section <- function(){
  rows = 20
  cols = 8
  hex_df <- hex_grid(rows, cols, size = 3)
  
  
  center = hex_df%>%
    filter(id == paste0(round(rows/2),"-",round(cols/2))) %>% 
    dplyr::group_by()%>%
    dplyr::summarise(cx = mean(x), cy = mean(y), .groups="drop")
  
  cell = hex_df%>%
    dplyr::group_by(id)%>%
    dplyr::summarise(mx = mean(x), my = mean(y),
                     .groups="drop")%>%
    mutate(euc = sqrt((mx-center$cx)^2+(my-center$cy)^2),
           atan = atan2(my, mx))
  
  hypo_cell = hex_df  %>%
    left_join(cell, by = "id")
  
  names = c("0", rep("1", 6), 
            rep("2", 2*6),
            rep("3", 3*6), 
            rep("4", 4*6))
  
  hypo_cell = hypo_cell %>% 
    filter(euc < 22) %>% 
    arrange(euc, atan, id) %>% 
    mutate(type = sort(rep(names, 7)))
  
  hypo_cell %>% ggplot(aes(hex_x,hex_y))+
    geom_polygon(aes(group = id), colour = "white", size = 3, alpha =0.2, data =hex_df)+
    geom_polygon(aes(group = id), colour = "white", size = 3)+
    geom_point()+
    coord_fixed()+theme_classic()
  
  radii = c(1:17,21)
  # radii = c(21)
  
  snapp = snap_to_circle(hypo_cell %>% mutate(x = hex_x -center$cx,y = hex_y - center$cy ),
                         radii)
  binder = tibble(id = unique(snapp$id), id_cell = 1:length(unique(snapp$id)))
  snapped_data = snapp
  df_hypo = snapped_data %>% mutate(x = snapped_x, y = snapped_y) %>% 
    left_join(binder, by = "id")
  df_hypo %>% ggplot(aes(snapped_x,snapped_y))+
    geom_polygon(aes(group = id_cell), colour = "white", size = 3)+
    geom_point()+
    coord_fixed()
  
  id_cell_vector = unique(df_hypo$id_cell)
  
  hypo_cell = NULL
  for(i in id_cell_vector){
    
    tmp = df_hypo %>% filter(id_cell == i) #%>% 
    #mutate(id_point = paste0(round(closest_r), ";", round(angle*12))) %>% 
    #filter(!duplicated(id_point))
    
    sf_linestring <- st_sfc(st_linestring(as.matrix(rbind(tmp[, c("x", "y")],tmp[1, c("x", "y")]))), crs = 2056)
    my_multilinestring = sf::st_sf(geom = sf::st_sfc(sf_linestring), crs = 2056)
    r_poly_ruff <-  sf::st_union(my_multilinestring)%>% sf::st_polygonize() %>% sf::st_collection_extract()
    r_poly_smooth <- smoothr::smooth(r_poly_ruff, method = "ksmooth", smoothness = 0.5)
    shrunken_polygon <- st_buffer(r_poly_smooth, -0.2)
    
    ggplot() +
      geom_sf(data = r_poly_ruff, fill = "lightblue", color = "black", alpha = 0.5) +
      geom_sf(data = r_poly_smooth, fill = "red", color = "darkred", alpha = 0.5) +
      geom_sf(data = shrunken_polygon, fill = "darkblue", color = "blue", alpha = 0.5) +
      ggtitle("Resized Polygons with Holes") +
      theme_minimal()
    
    # Get the coordinates
    coords <- sf::st_coordinates(shrunken_polygon )[, 1:2]
    
    pol = tibble(x = coords[,1], y = coords[,2], 
                 id_cell = i)
    pol = pol %>% 
      mutate(x1 = x,
             y1 = y,
             x2 = c(pol$x[-1], pol$x[1]),
             y2 = c(pol$y[-1], pol$y[1]))
    
    hypo_cell = rbind(hypo_cell, pol)
    
    
  }
  
  
  center = hypo_cell%>%
    dplyr::group_by()%>%
    dplyr::summarise(cx = mean(x), cy = mean(y), .groups="drop")
  
  cell = hypo_cell%>%
    dplyr::group_by(id_cell)%>%
    dplyr::summarise(mx = mean(x), my = mean(y),
                     .groups="drop")%>%
    mutate(euc = sqrt((mx-center$cx)^2+(my-center$cy)^2),
           atan = atan2(my, mx))
  
  hypo_cell = hypo_cell %>%
    left_join(cell, by = "id_cell")
  hypo_cell%>% 
    ggplot(aes(x,y))+
    #geom_polygon(aes(x,y), data = pol_out)+
    geom_polygon(aes(fill = id_cell, group = factor(id_cell)))+
    geom_point(aes(x,y), data = hypo_cell %>%  mutate(dist= sqrt((x)^2+(y)^2)) %>% 
                 filter(dist > 20.6))+
    #geom_path(data = cir(hypo_cell, 142), aes(x, y), color = "red", linewidth = 1) +  # Add circle
    coord_fixed()+
    theme_classic()+
    viridis::scale_fill_viridis(option = "A")
  
  outer = hypo_cell %>%  mutate(dist= sqrt((x)^2+(y)^2)) %>% 
    filter(dist > 20.6) %>% 
    mutate(atan = atan2(y, x)) %>% arrange(atan)
  snapped_outer <- snap_to_circle(outer, radii = 21.2)
  
  sf_linestring <- st_sfc(st_linestring(as.matrix(rbind(snapped_outer[, c("snapped_x", "snapped_y")],snapped_outer[1, c("snapped_x", "snapped_y")]))), crs = 2056)
  plot(sf_linestring)
  my_multilinestring = sf::st_sf(geom = sf::st_sfc(sf_linestring), crs = 2056)
  r_poly_smooth <-  sf::st_union(my_multilinestring)%>% sf::st_polygonize() %>% sf::st_collection_extract()
  out_wall <- sf::st_coordinates(r_poly_smooth)[, 1:2]
  pol_out = tibble(x = out_wall[,1], y = out_wall[,2], 
                   id_cell = max(hypo_cell$id_cell)+1)
  pol_out = pol_out%>% 
    mutate(x1 = x,
           y1 = y,
           x2 = c(pol_out$x[-1], pol_out$x[1]),
           y2 = c(pol_out$y[-1], pol_out$y[1]), 
           mx = 0, my = 0,
           euc = sqrt(x^2+y^2),
           atan = atan2(y, x))
  
  
  rbind (hypo_cell) %>% 
    ggplot(aes(x,y))+
    geom_polygon(aes(x,y), data = pol_out)+
    geom_polygon(aes(fill = -id_cell, group = factor(id_cell)))+
    #geom_point(colour = "white", size = 5, alpha = 0.5)+
    #geom_path(data = circle_data, aes(x, y), color = "red", linewidth = 1) +  # Add circle
    coord_fixed()+
    theme_classic()+
    viridis::scale_fill_viridis(option = "A")
  
  
  
  write_geo(rbind (pol_out, hypo_cell)%>%mutate(res = 1),dim = 2, path_geo = "~/GitHub/VergerLab/2D-HypocotylFEM/data/in/hypocross05.geo", extrude = 0)
}


write_geo <- function(data, dim = 2, path_geo, averaging = FALSE, extrude=0){
  
  date = Sys.time()
  x1 = paste0('// Gmsh project created on ', date,'\nSetFactory("OpenCASCADE");\n//+\n')
  
  if(dim == 2){
    data$z = 0
    data$z1 = 0
    data$z2 = 0
  }
  
  k = h = j = 1
  txt = x1
  for(i in sort(unique(data$id_cell))){
    tmp = data%>%filter(id_cell == i)
    i_point = unique_point_id(tmp)
    
    i_point$idx = seq(k,-1+k+nrow(i_point),1)
    i_point$idx2 = c(i_point$idx[-1],i_point$idx[1])
    
    if(i < max(data$id_cell)){
      dots = paste0("Point(",i_point$idx,") = {",round(i_point$x,4), ' ,',round(i_point$y,4),' ,',round(i_point$z,4), ',1.0};\n', collapse = "")
      Lines = paste0("Line(",i_point$idx,") = {",i_point$idx, ", ", i_point$idx2,'};\n//+\n', collapse = "")
      Surf = paste0("Curve Loop(",j,") = {",paste0(i_point$idx, collapse = ", "),'};\n//+\nSurface(',
                    h,") = {", j,'};\n//+\n', collapse = "")
      Physical = paste0("//Physical Surface(",h,") = {",h,"};\n")
      
      txt = paste0(txt, dots, Lines, Surf, Physical)
    }else{
      dots = paste0("Point(",i_point$idx,") = {",round(i_point$x,4), ' ,',round(i_point$y,4),' ,',round(i_point$z,4), ',0.5};\n', collapse = "")
      Lines = paste0("Line(",i_point$idx,") = {",i_point$idx, ", ", i_point$idx2,'};\n//+\n', collapse = "")
      Surf = paste0("Curve Loop(",j,") = {",paste0(i_point$idx, collapse = ", "),'};\n//+\nPlane Surface(',
                    h,") = {", paste0(sort(seq(1,j,2), decreasing = T), collapse = ", "),'};\n//+\n', collapse = "")
      Physical = paste0("Physical Surface(0) = {",h,"};\n")
      
      txt = paste0(txt, dots, Lines, Surf, Physical)
    }
    
    
    
    
    h = h+1
    j = j+2
    k = max(i_point$idx)+1
  }
  write(txt, path_geo)
}


