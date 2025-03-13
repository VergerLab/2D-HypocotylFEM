
#' @title 2D-Finite Element Modeling (FEM) of hypocotyl epidermal cells
#' @author Adrien Heymans, Özer Erguvan, Stéphane Verger
#' @description
#' Analysis script of the 3 meshes
#' A) 2D mesh representing the surface of a hypocotyl epidermis with staggered 
#' cell files, segmented into two subdomains: the cell and the interface 
#' between cells.
#' 
#' B) 2D mesh of a hypocotyl longitudinal section with two subdomains: the cell 
#' wall and the adhesive layer at the interface.
#' 
#' c) 2D mesh of a hypocotyl longitudinal section with four subdomains: the 
#' Supracellular Outer Epidermal Wall (SOEW), the Outer Epidermal Edge Filling 
#' (OEEF), the Inner Walls, and the Middle Lamella (ML)
#' @date
#' novembre 2024
#' 

# Library import
library(tidyverse)
library(data.table)
# Folder path
setwd("~/GitHub/VergerLab/2D-HypocotylFEM/")
source("./src/io_function/io_function.R")


################
# -- mesh 0 -- #  
################

# Run python script with the FEM
# 2D mesh representing the cell wall of a cross section of a plant tissue.

orga_cross = create_cross_section()

orga_cross%>%mutate(id_cell = ifelse(id_cell == 62, 0, id_cell)) %>% 
  ggplot(aes(x,y))+
  geom_polygon(aes(fill = -id_cell, group = factor(id_cell)))+
  coord_fixed()+
  theme_classic()+
  viridis::scale_fill_viridis(option = "A")

system("python ./src/script/cross.py")


path_csv = "./data/out/cross/csv/"
fls = list.files(path_csv)

data_cross_raw = read.csv(paste0(path_csv, fls[1]))
data_cross = stress_strain(data_cross_raw)
data_cross = concatenate_res(data_cross)

tmp1 = data_cross %>%
  mutate(x1 = x + pcstress_main/100 * sin(angle_stress_mod),
         x2 = x - pcstress_main/100 * sin(angle_stress_mod),
         y1 = y + pcstress_main/100 * cos(angle_stress_mod),
         y2 = y - pcstress_main/100 * cos(angle_stress_mod))

tmp2 = data_cross%>%
  mutate(x1 = x + pcstress_secd/100 * cos(angle_stress_mod),
         x2 = x - pcstress_secd/100 * cos(angle_stress_mod),
         y1 = y - pcstress_secd/100 * sin(angle_stress_mod),
         y2 = y + pcstress_secd/100 * sin(angle_stress_mod))

glyph_to_obj(tmp1, path = paste0("./data/out/cross/glyphs/cross_stress1.obj"))
glyph_to_obj(tmp2, path = paste0("./data/out/cross/glyphs/cross_stress2.obj"))

data_cross = data_cross%>%
  mutate(euc = sqrt(x^2+y^2))

ggplot()+
  geom_point(aes(x,y), data = data_cross)+
  geom_point(aes(x,y, colour = stress_magnitude), data = data_cross %>% 
               filter(euc >= 22))+
  coord_fixed()+
  viridis::scale_colour_viridis()

ggplot()+
  geom_point(aes(euc,stress_magnitude), data = data_cross)+
  viridis::scale_colour_viridis()+
  labs(x = "distance from center [µm]",
       y = "Stress Magnitude [MPa]")+
  xlim(3,23.5)+
  theme_classic()
ggsave("./data/out/img/Cross_StressQ.svg")


################
# -- mesh 1 -- #  
################

# Run python script with the FEM
# 2D mesh representing the surface of a hypocotyl epidermis with staggered cell files, 
# segmented into two subdomains: the cell and the interface between cells.
# 2 scenarios: 
# - soft ML & stiff CW
# - stiff ML & soft CW

system("python ./src/script/Hypocot_epi_surface.py")


path_csv = "./data/out/Epi_Surface/csv/"
fls = list.files(path_csv)
fls = fls[grepl(pattern = "_a", fls)]

data_surface = NULL
for(fl in fls){
  
  tmp = read.csv(paste0(path_csv, fl))
  config = unlist(str_split(fl, "_"))
  simu_type = ifelse(config[3]== "softML", "Stiff cells, soft interface","Soft cells, stiff interface")
  
  tmp = stress_strain(tmp)
  tmp = concatenate_res(tmp)
  tmp = tmp %>%
    mutate(simu = simu_type)
  
  data_surface = rbind(data_surface, tmp)
}

data_surface%>%filter(simu == "Stiff cells, soft interface")%>%
  filter(((y < 85 & y > 75) &(x >= 20.4 & x <= 20.8)) | ((x > 20.4 & x < 31.4) &(y >= 78.8 & y <= 79.2))|
           (x > 20.8 & x < 20.9 & y > 78.7 & y < 79.3) | ((y < 85 & y > 75) &(x >= 30.7 & x < 31.3)) )%>%
  ggplot(aes(x,y))+
  geom_point(aes(colour =stress_magnitude ))+
  coord_fixed()+theme_classic()
  
  

# Write the main and minor axis of the stress/strain anistropy
for(j in unique(data_surface$simu)){
  tmp1 = data_surface%>%filter(simu == j)%>%
    filter(((y < 85 & y > 75) &(x >= 20.4 & x <= 20.8)) | ((x > 20.4 & x < 31.4) &(y >= 78.8 & y <= 79.2))|
             (x > 20.8 & x < 20.9 & y > 78.7 & y < 79.3) | ((y < 85 & y > 75) &(x >= 30.7 & x < 31.3)) )%>%
    mutate(x1 = x + pcstress_main/3000 * sin(angle_stress_mod),
           x2 = x - pcstress_main/3000 * sin(angle_stress_mod),
           y1 = y + pcstress_main/3000 * cos(angle_stress_mod),
           y2 = y - pcstress_main/3000 * cos(angle_stress_mod))
  
  tmp2 = data_surface%>%filter(simu == j)%>%
    filter(((y < 85 & y > 75) &(x >= 20.4 & x <= 20.8)) | ((x > 20.4 & x < 31.4) &(y >= 78.8 & y <= 79.2))|
             (x > 20.8 & x < 20.9 & y > 78.7 & y < 79.3) | ((y < 85 & y > 75) &(x >= 30.7 & x < 31.3)) )%>%
    mutate(x1 = x + pcstress_secd/3000 * cos(angle_stress_mod),
           x2 = x - pcstress_secd/3000 * cos(angle_stress_mod),
           y1 = y - pcstress_secd/3000 * sin(angle_stress_mod),
           y2 = y + pcstress_secd/3000 * sin(angle_stress_mod))
  
  glyph_to_obj(tmp1, path = paste0("./data/out/Epi_Surface/glyphs/",j,"_stress1.obj"))
  glyph_to_obj(tmp2, path = paste0("./data/out/Epi_Surface/glyphs/",j,"_stress2.obj"))
  
  tmp1 = data_surface%>%filter(simu == j)%>%
    filter(((y < 85 & y > 75) &(x >= 20.4 & x <= 20.8)) | ((x > 20.4 & x < 31.4) &(y >= 78.8 & y <= 79.2))|
             (x > 20.8 & x < 20.9 & y > 78.7 & y < 79.3) | ((y < 85 & y > 75) &(x >= 30.7 & x < 31.3)) )%>%
    mutate(x1 = x + pcstrain_main * sin(angle_strain_mod),
           x2 = x - pcstrain_main * sin(angle_strain_mod),
           y1 = y + pcstrain_main * cos(angle_strain_mod),
           y2 = y - pcstrain_main * cos(angle_strain_mod))
  
  tmp2 = data_surface%>%filter(simu == j)%>%
    filter(((y < 85 & y > 75) &(x >= 20.4 & x <= 20.8)) | ((x > 20.4 & x < 31.4) &(y >= 78.8 & y <= 79.2))|
             (x > 20.8 & x < 20.9 & y > 78.7 & y < 79.3) | ((y < 85 & y > 75) &(x >= 30.7 & x < 31.3)) )%>%
    mutate(x1 = x + pcstrain_secd * cos(angle_strain_mod),
           x2 = x - pcstrain_secd * cos(angle_strain_mod),
           y1 = y - pcstrain_secd * sin(angle_strain_mod),
           y2 = y + pcstrain_secd * sin(angle_strain_mod))
  
  glyph_to_obj(tmp1, path = paste0("./data/out/Epi_Surface/glyphs/",j,"_strain1.obj"))
  glyph_to_obj(tmp2, path = paste0("./data/out/Epi_Surface/glyphs/",j,"_strain2.obj"))
  
}


################
# -- mesh 2 -- #  
################

# Run python script with the FEM
# 2D mesh of a hypocotyl longitudinal section with two subdomains: 
# the cell wall and the adhesive layer at the interface.
# 1 scenarios: 
# - soft ML & stiff CW

system("python ./src/script/Hypocot_epi_2buldgedCells.py")

# Find files names
path_csv = "./data/out/2BuldgedCells/csv/"
fls = list.files(path_csv)
fls = fls[grepl(pattern = "_a", fls)]

# Gather data
data_buldged = NULL
for(fl in fls){
  
  tmp = read.csv(paste0(path_csv, fl))
  # Rep is related to the mesh resolution
  rep =file.info(paste0(path_csv, fl))$size
  
  # Tidy the data
  tmp = stress_strain(tmp)
  tmp = concatenate_res(tmp)
  buldged = tmp%>%mutate(rep = rep)
  data_buldged = rbind(data_buldged, buldged)
}

buldV_data =data_buldged %>%
  filter(x >= 20.2, x <= 20.3, rep == max(data_buldged$rep) )%>%
  mutate(y = round(y*100)/100) %>% 
  dplyr::group_by(y)%>%
  dplyr::summarise(angle_stress = 180-mean(abs(angle_stress_mod)/(pi))*180,
                   pcstress1 = mean(pcstress_main),
                   stress_x = mean(stress_x)/sqrt(2),
                   stress = mean(stress_magnitude),
                   ani_stress = mean(ani_stress),
                   angle_strain = 180-mean(abs(angle_strain_mod)/(pi))*180,
                   pcstrain1 = mean(pcstrain_main),
                   strain_x = mean(strain_x)/sqrt(2),
                   strain = mean(strain_magnitude),
                   ani_strain = mean(ani_strain),.groups = "drop")

buldV = rbind(tibble(s = "Stress",y = buldV_data$y, unit = "[°]",value = buldV_data$angle_stress, type = "Main tensor direction [°]"),
              # tibble(s = "Stress",y = buldV_data$y, unit = "[MPa]",value = buldV_data$stress, type = 'Magnitude [MPa]'),
              tibble(s = "Stress",y = buldV_data$y, unit = " [MPa]",value = buldV_data$stress, type = 'Magnitude [MPa]'),
              #tibble(s = "Stress",y = buldV_data$y, unit = "[-]",value = buldV_data$ani_stress, type = 'Anisotropy [-]'),
              tibble(s = "Strain",y = buldV_data$y, unit = "[°]",value = buldV_data$angle_strain, type = "Main tensor direction [°]"),
              #  tibble(s = "Strain",y = buldV_data$y, unit = "[.]",value = buldV_data$strain, type = 'Magnitude'),
              tibble(s = "Strain",y = buldV_data$y, unit = "[.]",value = buldV_data$strain, type = 'Magnitude [-]')
              #tibble(s = "Strain",y = buldV_data$y, unit = "[-]",value = buldV_data$ani_strain, type = 'Anisotropy [-]')
)

buldV%>%
  ggplot(aes(y, value, colour = type))+
  geom_line(size = 1.2)+
  facet_grid(~unit, scale = 'free')+
  theme_classic()+coord_flip()+viridis::scale_colour_viridis(discrete = T, option = "D")
ggsave("./data/out/img/2BuldgedCells_anticlinal_stress-strain-angle.png")

buldV%>%
  filter(y > 9.2) %>% 
  ggplot(aes(y, value, colour = type))+
  geom_line(size = 1.2)+
  facet_grid(~unit, scale = 'free')+
  theme_classic()+coord_flip()+viridis::scale_colour_viridis(discrete = T, option = "D")
ggsave("./data/out/img/2BuldgedCells_anticlinal_Zoom_stress-strain-angle.png")


# Write the main and minor axis of the stress anistropy
tmp_1 = data_buldged%>%
  filter(rep == min(data_buldged$rep)) %>% 
  mutate(pcstress_main = ifelse(pcstress_main > 150, 150, pcstress_main))%>%
  mutate(x1 = x + pcstress_main/3000* sin(angle_stress),
         x2 = x - pcstress_main/3000* sin(angle_stress),
         y1 = y + pcstress_main/3000 * cos(angle_stress),
         y2 = y - pcstress_main/3000 * cos(angle_stress))
tmp_2 = data_buldged%>%
  filter(rep == min(data_buldged$rep)) %>% 
  mutate(x1 = x + pcstress_secd/3000 * cos(angle_stress_mod),
         x2 = x - pcstress_secd/3000 * cos(angle_stress_mod),
         y1 = y - pcstress_secd/3000 * sin(angle_stress_mod),
         y2 = y + pcstress_secd/3000 * sin(angle_stress_mod))

glyph_to_obj(tmp_1, path = "./data/out/2BuldgedCells/glyphs/buldged_cells_stress_1.obj")
glyph_to_obj(tmp_2, path = "./data/out/2BuldgedCells/glyphs/buldged_cells_stress_2.obj")

# Write the main and minor axis of the strain anistropy
tmp_1 = data_buldged%>%
  filter(rep == min(data_buldged$rep)) %>% 
  mutate(x1 = x + pcstrain_main * sin(angle_strain_mod),
         x2 = x - pcstrain_main * sin(angle_strain_mod),
         y1 = y + pcstrain_main * cos(angle_strain_mod),
         y2 = y - pcstrain_main * cos(angle_strain_mod))
tmp_2 = data_buldged%>%
  filter(rep == min(data_buldged$rep)) %>% 
  mutate(x1 = x + pcstrain_secd * cos(angle_strain_mod),
         x2 = x - pcstrain_secd * cos(angle_strain_mod),
         y1 = y - pcstrain_secd * sin(angle_strain_mod),
         y2 = y + pcstrain_secd * sin(angle_strain_mod))

glyph_to_obj(tmp_1, path = "./data/out/2BuldgedCells/glyphs/buldged_cells_strain_1.obj")
glyph_to_obj(tmp_2, path = "./data/out/2BuldgedCells/glyphs/buldged_cells_strain_2.obj")


################
# -- mesh 3 -- #  
################

# Run python script with the FEM
# 2D mesh of a hypocotyl longitudinal section with four subdomains: the 
# Supracellular Outer Epidermal Wall (SOEW), the Outer Epidermal Edge Filling 
# (OEEF), the Inner Walls, and the Middle Lamella (ML)
# 5 scenarios: 
# - Uniform
# - soft ML
# - soft OEEF
# - soft SOEW
# - soft SOEW & OEEF

# Input data frame for the simulation related to this 2D mesh
path_sim = "./data/in/sim_paper.csv"
sim1 = tibble(simu = "Uniform", Cuti = 10000, IW = 10000, SOEW = 10000, OEEF = 10000, ML = 10000)
sim2 = tibble(simu = "Weak_ML", Cuti = 10000, IW = 10000, SOEW = 10000, OEEF = 10000, ML = 1000)
sim3 = tibble(simu = "Weak_OEEF", Cuti = 10000, IW = 10000, SOEW = 10000, OEEF = 1000, ML = 10000)
sim4 = tibble(simu = "Weak_SOEW", Cuti = 1000, IW = 10000, SOEW = 1000, OEEF = 10000, ML = 10000)
sim5 = tibble(simu = "Weak_SOEW_OEEF", Cuti = 1000, IW = 10000, SOEW = 1000, OEEF = 1000, ML = 10000)
sim = rbind(sim5, sim2,sim3, sim4, sim1)%>%
  mutate(SOEW_thickness = 0.3,
         TP = 0.2, Force = 0.25)
# Write the input data frame
write.csv(sim, path_sim)

system("python ./src/script/Hypocot_epi_SOEW.py")
############################################

# Find files names
path_csv = "./data/out/SOEW/csv/"
fls = list.files(path_csv)
fls = fls[grepl(pattern = "_a", fls)]
fls = fls[grepl(pattern = "sim", fls)]

# Gather data
data = NULL
for(fl in fls){
  
  tmp = read.csv(paste0(path_csv, fl))
  config = unlist(str_split(fl, "_"))
  simu_type = ifelse(config[2]== "Uniform", "Uniform",
                     ifelse(config[3] != "SOEW", paste0("Weak ",config[3]),
                            ifelse(config[4] == "OEEF", "Weak SOEW and OEEF", 
                                   "Weak SOEW")))
  # Rep is related to the mesh resolution
  rep =parse_integer(config)
  rep = rep[!is.na(rep)]
  # Tidy the data
  tmp = stress_strain(tmp)
  tmp = concatenate_res(tmp)%>%
    mutate(simu = simu_type,
           rep = rep)
  data = rbind(data, tmp)
}

# make strain main glyphs for the low res meshes
for(j in unique(data$simu)){
  tmp = data%>%filter(simu == j, rep == 30)%>%
    mutate(x1 = x + 10*pcstrain_main * sin(angle_strain_mod),
           x2 = x - 10*pcstrain_main * sin(angle_strain_mod),
           y1 = y + 10*pcstrain_main * cos(angle_strain_mod),
           y2 = y - 10*pcstrain_main * cos(angle_strain_mod))
  glyph_to_obj(tmp, path = paste0("./data/out/SOEW/glyphs/",j,"_strain.obj"))
}
# make stress main glyphs for the low res meshes
for(j in unique(data$simu)){
  tmp = data%>%filter(simu == j, rep == 30)%>%
    mutate(x1 = x + pcstress_main/1000 * sin(angle_stress_mod),
           x2 = x - pcstress_main/1000 * sin(angle_stress_mod),
           y1 = y + pcstress_main/1000 * cos(angle_stress_mod),
           y2 = y - pcstress_main/1000 * cos(angle_stress_mod))
  glyph_to_obj(tmp, path = paste0("./data/out/SOEW/glyphs/",j,"_stress.obj"))
  
}

# Stress & strain data within the cell junction area.
# This area encompass nodes within 0.1625 µm away from the median line that 
# separate the two cells (Inner walls and ML).

data_2C = data  %>%filter(x >= 6, x <= 6.375)%>%
  dplyr::group_by(simu, y = round(y*10)/10)%>%
  dplyr::summarise(stress_x = mean(stress_x/sqrt(2)),
                   strain_x = mean(strain_x/sqrt(2)),
                   Vstress = mean(VonMises_stress),
                   stress = mean(stress_magnitude),
                   Vstrain = mean(VonMises_strain),
                   strain = mean(strain_magnitude),.groups = "drop")

two_C = rbind(tibble(y = data_2C$y,
                     simu = data_2C$simu,
                     value = data_2C$stress_x, 
                     what = "Separating stress [MPa]"),
              tibble(y = data_2C$y,
                     simu = data_2C$simu,
                     value = data_2C$strain_x, 
                     what = "Separating strain [-]"),
              tibble(y = data_2C$y,
                     simu = data_2C$simu,
                     value = data_2C$stress, 
                     what = "Stress magnitude [MPa]"),
              tibble(y = data_2C$y,
                     simu = data_2C$simu,
                     value = data_2C$strain, 
                     what = "Strain magnitude [-]"))
# Along the anticlinal axis

two_C %>%
  ggplot( aes(x = y))+
  geom_line(aes(y = value, colour = simu), size = 1.2)+
  geom_vline(aes(xintercept = 12.3), linetype = 2)+
  geom_vline(aes(xintercept = 9), linetype = 2)+
  geom_vline(aes(xintercept = 12.6), linetype = 2)+
  labs(x = "Anticlinal axis [µm]",
       colour = "Simulated scenarios:") +
  facet_wrap(~what, scales = "free")+
  coord_flip()+
  theme_dark()+
  viridis::scale_colour_viridis(discrete = T)

ggsave("./data/out/img/SOEW_anticlinal_stress-strain.png")

# Stress & strain data within the SOEW region over the cell junction.
data_2C_x = data  %>%filter(y >= 12.3)%>%
  dplyr::group_by(simu, x = round(x*5)/5-6.1625)%>% # center over the junction
  dplyr::summarise(stress_x = mean(stress_x/sqrt(2)),
                   strain_x = mean(strain_x/sqrt(2)),
                   Vstress = mean(VonMises_stress),
                   stress = mean(stress_magnitude),
                   Vstrain = mean(VonMises_strain),
                   strain = mean(strain_magnitude),.groups = "drop")

two_C_x = rbind(tibble(x = data_2C_x$x,
                     simu = data_2C_x$simu,
                     value = data_2C_x$stress_x, 
                     what = "Separating stress [MPa]"),
              tibble(x = data_2C_x$x,
                     simu = data_2C_x$simu,
                     value = data_2C_x$strain_x, 
                     what = "Separating strain [-]"),
              tibble(x = data_2C_x$x,
                     simu = data_2C_x$simu,
                     value = data_2C_x$stress, 
                     what = "Stress magnitude [MPa]"),
              tibble(x = data_2C_x$x,
                     simu = data_2C_x$simu,
                     value = data_2C_x$strain, 
                     what = "Stress magnitude [-]"))

two_C_x %>%
  ggplot( aes(x = x))+
  geom_line(aes(y = value, colour = simu), size = 1)+
  labs(x = "SOEW, Periclinal axis [µm]",
       colour = "Simulated scenarios:") +
  facet_wrap(~what, scales = "free")+
  theme_dark()+
  viridis::scale_colour_viridis(discrete = T)

ggsave("./data/out/img/SOEW_periclinal_stress-strain.png")






