####################################################################################################################################
# SCRIPT 8

# Getting started with ABC simulations
# 2024-12-04
# Barbara Walkowiak bw18

# Script to run simulations (Henry wants to do this on a grid so I guess this is what we are going to try)

###################################################################################################################################
# LIBRARIES
library(data.table)
library(ggplot2)
library(plyr)

###################################################################################################################################
# WHAT I WANT TO DO

# initialize a large grid
# the grid (2D) represents the developing embryo

# set size of the grid
gridX = 30
gridY = 30

# Number of starting cells 
start_cells = 1 # how many cells to start with (zygote so 1 cell)
start_div = 4 # number of divisions before cells selected to the ICM
cells_to_icm = 3 # how many cells to select to ICM
icm_div = 4 # number of times each ICM cell divides once the ICM has formed

# Start with one cell on the grid

# Divide the cell n times, place the descendants on the grid (not completely random)

# Select cells to the ICM and keep dividing those


adjacells <- getAdjacentCellsDataFrame()


# Make multiple runs (Replication of simulation) and take the average of stats
st <- data.frame()
st_row<- vector()
for(i in 1:kNumReplications) {
  area_df <- resetIteration()
  
  seedAreaWithPioneers(numPioneers,seeding.opt)
  simstats <- accommodateSettlers(kNumSettlers, settling.option)
  
  found.home <- simstats[1]
  max.look.around <- simstats[2]
  #compute for this iterations
  st_row <- store_iteration_stats(i, kNumSettlers, found.home, max.look.around)
  st <- rbind(st,st_row)
}

#Render
p <- drawArea(area_df)
p
names(st) <- c("Iter", "FoundHome", "NumSettlers", "Percent")
st











