####################################################################################################################################
# SCRIPT 8

# Getting started with ABC simulations
# 2024-12-04
# Barbara Walkowiak bw18

# Script to run simulations (Henry wants to do this on a grid so I guess this is what we are going to try)
# I found some useful code here: https://gist.github.com/Ram-N/5019098 

###################################################################################################################################
# LIBRARIES
library(data.table)
library(ggplot2)
library(plyr)

###################################################################################################################################
# SIMULATION SUMMARY 

# PARAMETERS 
# Number of cells selected to the ICM
# Stage of cell division when cells selected to the ICM
# Number of cell divisions in the ICM
# Asymmetry of split
# Degree of cell mixing 

# ASSUMPTIONS
# A lot
# Cells selected to the ICM once 
# Each cell divides the same number of times 
# No cell death, no cell arrest etc. 
# Spatial structure accurately approx by a 2D grid 

# SUMMARY STATISTICS 
# TBD, we don't know 

###################################################################################################################################
# SIMULATION: FUNCTIONS

# Define necessary functions

# creates an empty dataframe and selects at random the first cell to start 
start_sim = function(x, y, center = TRUE){ # create empty dataframe with 0s that will be further populated 
  
  # create an empty grid 
  grid = data.frame( matrix(0, nrow=x, ncol=y))
  names(grid) = c(1:x)  
  
  # initialize the starting cell in the middle of the grid 
  if (center == TRUE){
    xNew = x/2
    yNew = y/2
    grid[xNew, yNew] = 1 # 1 indicates there is a cell there 
    return(grid)
  }
  
  else{
    # initialize the starting cell wherever
    xNew = sample(1:x, 1) # get new coordinates (x direction)
    yNew = sample(1:y, 1) # get new coordinates (y direction) 
    grid[xNew, yNew] = 1 # 1 indicates there is a cell there 
    return(grid)
  }
  
}

# visualise grid 
draw_area = function(grid){  
  
  df = reshape2::melt(as.matrix(grid, nrow=dim(grid)[1])) 
  df$value = as.factor(df$value)
  square = 15
  p = ggplot(df, aes(Var1, Var2, color=value)) + 
    geom_point(shape=square, size=6) + 
    theme_classic() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    scale_colour_manual(values = c("black","blue", "#10c8c5", "purple", "darkred")) # black if absent, red if present 
  return(p)
}

# divide all cells on the grid (for each existing cell, create two daughter cells, then remove the existing one)
# use this function when dividing cells to contribute to the ICM 

# for each cell existing on the grid, this function selects coordinates of two new cells
# check that coordinates are not outside of grid bounds and do not overlap 

divide_cells_pre_icm = function(grid, sd = 5){
  
  # count the nr of cells on the grid
  nr_cells = sum(grid==1) 
  print(paste('Number of cells to divide:', nr_cells))
  
  # find grid dimensions
  x = dim(grid)[1]
  y = dim(grid)[2]
  
  # for each cell, select coordinates of the new daughter cell 
  for (i in 1:nr_cells){

    # coordinates of the parent cell 
    xParent = which(grid==1, arr.ind=TRUE)[i, 1]
    yParent = which(grid==1, arr.ind=TRUE)[i, 2]
    
    # sample coordinates for new cells 
    # the idea is to sample around the same row and around the same column
    # this creates a spatial structure on the grid 
    # ensure that cells do not overlap 
    
    # sample coordinates for the first cell
    repeat{ # repeat until coordinate sampled 
      xNew1 = rnorm(1, mean = xParent, sd = sd) 
      yNew1 = rnorm(1, mean = yParent, sd = sd)
      if (xNew1 > 1 & xNew1 < x & yNew1 > 1 & yNew1 < y){ # check coordinates not outside the grid
        if (grid[xNew1, yNew1] == 0){ # add the cell if the position is not occupied
          grid[xNew1, yNew1] = 1 
          break 
        } 
      }
    }
    
    # sample coordinates for the second cell
    repeat{ # repeat until coordinate sampled 
      xNew2 = rnorm(1, mean = xParent, sd = sd) 
      yNew2 = rnorm(1, mean = yParent, sd = sd) 
      if (xNew1 > 1 & xNew1 < x & yNew1 > 1 & yNew1 < y){ # check coordinates not outside the grid 
        if (grid[xNew2, yNew2] == 0){ # add the cell if the position is not occupied 
          grid[xNew2, yNew2] = 1 
          break 
        }
      }
    }
    
    # remove the mother cell from the grid 
    grid[xParent, yParent] = 0
    
  } 
  
  nr_cells_added = sum(grid==1) 
  expected_nr_cells_added = nr_cells * 2
  
  print(paste('Number of added cells:', nr_cells_added))
  print(paste('Expected number of added cells:', expected_nr_cells_added))
  
  return(grid)
  draw_area(grid)
  
}

# select cells from existing ones (select to ICM)
select_to_icm = function(grid, n){ # select n cells from all existing cells 
  
  # count the nr of cells on the grid
  nr_cells = sum(grid==1) 
  print(paste0('Number of cells available to go to the ICM:', nr_cells))
  
  # obtain coordinates of all cells which are present 
  coords_available = which(grid==1, arr.ind=TRUE)
  
  coords_to_icm = sample(nrow(coords_available), size = n) # sample n random rows from the matrix (= select of n cells to icm)
  coords_icm = coords_available[c(coords_to_icm),]
  
  # change the value of those cells to 2 
  # will allow to show which cells got selected
  grid[coords_icm] = 2
  return(grid)
  
}

divide_cells_post_icm = function(grid, sd = 5){
  
  # count the nr of cells on the grid
  nr_cells = sum(grid==2) 
  print(paste('Number of ICM cells to divide:', nr_cells))
  
  # find grid dimensions
  x = dim(grid)[1]
  y = dim(grid)[2]
  
  # for each cell, select coordinates of the new daughter cell 
  for (i in 1:nr_cells){
    
    # coordinates of the parent cell 
    xParent = which(grid==2, arr.ind=TRUE)[i, 1]
    yParent = which(grid==2, arr.ind=TRUE)[i, 2]
    
    # sample coordinates for new cells 
    # the idea is to sample around the same row and around the same column
    # this creates a spatial structure on the grid 
    # ensure that cells do not overlap 
    
    # sample coordinates for the first cell
    repeat{ # repeat until coordinate sampled 
      xNew1 = rnorm(1, mean = xParent, sd = sd) 
      yNew1 = rnorm(1, mean = yParent, sd = sd)
      if (xNew1 > 1 & xNew1 < x & yNew1 > 1 & yNew1 < y){ # check coordinates not outside the grid
        if (grid[xNew1, yNew1] == 0){ # add the cell if the position is not occupied
          grid[xNew1, yNew1] = 3 
          break 
        } 
      }
    }
    
    # sample coordinates for the second cell
    repeat{ # repeat until coordinate sampled 
      xNew2 = rnorm(1, mean = xParent, sd = sd) 
      yNew2 = rnorm(1, mean = yParent, sd = sd) 
      if (xNew1 > 1 & xNew1 < x & yNew1 > 1 & yNew1 < y){ # check coordinates not outside the grid 
        if (grid[xNew2, yNew2] == 0){ # add the cell if the position is not occupied 
          grid[xNew2, yNew2] = 3
          break 
        }
      }
    }
    
  } 
  
  nr_cells_added = sum(grid==3) 
  expected_nr_cells_added = nr_cells * 2
  
  print(paste('Number of added cells:', nr_cells_added))
  print(paste('Expected number of added cells:', expected_nr_cells_added))
  
  return(grid)
  
}

# I guess you'd like to ensure that by the end of the cells dividing you have enough cells to split the grid and it's ~full
# so you can just draw a line through the grid or something of that kind 

# split the grid (cells to twin1 vs to twin2 )
split_twins = function(grid, p){ # specify proportion of cells that should go to twin1 
  
  # available cells will have the highest number
  total_nr_cells = sum(grid==max(grid))
  nr_cells_twin1 = round(p * total_nr_cells)
  sum_cells = 0
  
  for (i in 1:nrow(grid)){
    cells_in_row = which(grid[i,]==1)
    if (length(cells_in_row > 0)){
      for (j in cells_in_row){
        grid[i, j] = max(grid)
        sum_cells = sum_cells + 1
      } 
    }
    
    if (sum_cells >= nr_cells_twin1){
      return(grid)
    }
  }
}

###################################################################################################################################
# SIMULATION: parameters 

# specify parameters
x = 100
y = 100
sd = 5
s = 5 # rounds of cell division before ICM cells selected
n = 4 # number of cells to select to the ICM 
s2 = 3 # number of cell divisions of the ICM cells
p = 0.7

# run the simulation 
grid0 = start_sim(x, y)
draw_area(grid0)
grid1 = Reduce(divide_cells_pre_icm, 1:s, init = grid0)
draw_area(grid1)
grid2 = select_to_icm(grid1, n)
draw_area(grid2)
grid3 = Reduce(divide_cells_post_icm, 1:s2, init = grid2)
draw_area(grid3)
grid4 = split_twins(grid3, p)
draw_area(grid4)


###################################################################################################################################
# SIMULATION: run 

# we will be doing simulations so set seed so we can reproduce this 
set.seed(124)
start_df = create_empty_df(x, y)
draw_area(start_df)  
start_df1 = start_sim(start_df)
draw_area(start_df1)  
start_df2 = divide_cells_pre_icm(start_df1)
draw_area(start_df2) 
start_df3 = divide_cells_pre_icm(start_df2)
draw_area(start_df3) 
start_df4 = divide_cells(start_df3)
draw_area(start_df4) 
start_df5 = select_to_icm(3, start_df4) # select 3 cells 
draw_area(start_df4)  
draw_area(start_df5)
start_df6 = divide_icm_cells(start_df5, 1) # select 3 cells 
draw_area(start_df6)  

###################################################################################################################################
# SIMULATION: output / summary statistics 

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
p <- draw_area(area_df)
p
names(st) <- c("Iter", "FoundHome", "NumSettlers", "Percent")
st











