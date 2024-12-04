####################################################################################################################################
# SCRIPT 8

# Getting started with ABC simulations
# 2024-12-04
# Barbara Walkowiak bw18

# Script to run simulations (Henry wants to do this on a grid so I guess this is what we are going to try)
# I am trying to google how to do this and found some useful code here: https://gist.github.com/Ram-N/5019098 

###################################################################################################################################
# LIBRARIES
library(data.table)
library(ggplot2)
library(plyr)

###################################################################################################################################
# WHAT I WANT TO DO

# we will be doing simulations so set seed so this thing reproduces
set.seed(124)

# initialize a large grid
# the grid (2D) represents the developing embryo

# set size of the grid
x = 100
y = 100

# Number of starting cells 
start_cells = 1 # how many cells to start with (zygote so 1 cell)
start_div = 4 # number of divisions before cells selected to the ICM
cells_to_icm = 3 # how many cells to select to ICM
icm_div = 4 # number of times each ICM cell divides once the ICM has formed

# Define necessary functions

create_empty_df = function(x, y){ # create empty dataframe with 0s to populate it 
  area_df = data.frame( matrix(0, nrow=x, ncol=y))
  names(area_df) = c(1:x)  
  return(area_df)
}

start_df = create_empty_df(x, y)
drawArea(start_df) # okay so this is just showing you have an empty grid for now 

# start simulation by randomly placing a cell on the grid 
start_sim = function(grid){ # get new coordinates by sampling x and y coordinates at random 
  gridX = dim(grid)[1]
  gridY = dim(grid)[2]
  xNew = sample(1:gridX, 1) # get new coordinates (x direction)
  yNew = sample(1:gridY, 1) # get new coordinates (y direction) 
  grid[xNew, yNew] = 1 # 1 indicates there is a cell there 
  return(grid)
}

start_df1 = start_sim(start_df)
drawArea(start_df1) # okay so this is just showing you have an empty grid for now 

find_cell_coordinates = function(grid){
  x_coord = which(start_df1==1, arr.ind=TRUE)[1]
  y_coord = which(start_df1==1, arr.ind=TRUE)[2]
  coords = c(x_coord, y_coord)
  return(coords) # indices of where the cell is 
}

# divide all cells on the grid (create 2x that many cells and place them in a spatially-coherent manner)
# count how many cells you have and for each cell, create two new cells to show 

divide_cells = function(grid){
  
  # available coordinates 
  gridX = dim(grid)[1]
  gridY = dim(grid)[2]
  
  # count the nr of cells on the grid
  nr_cells = sum(grid==1) 
  print(paste0('Number of cells to divide:', nr_cells))
  
  # find coordinates of the existing cells

  # select coordinates for each of the new daughters 
  # each cell on the grid needs to give rise to two daughter cells
  
  coordinates = c() # list to store all coordinates in to make sure they are unique for each cell 
  
  for (i in 1:nr_cells){
    
    # coordinates of the parent cell 
    present_cell_x = which(grid==1, arr.ind=TRUE)[i, 1]
    present_cell_y = which(grid==1, arr.ind=TRUE)[i, 2]
    
    # we don't care if the cell is on its parent or not 
    # gridX = setdiff(1:gridX, present_cell_x)
    # gridY = setdiff(1:gridY, present_cell_x)
    
    # sample coordinates for new cells 
    # the idea is to sample around the same row and around the same column
    # this creates a spatial structure on the grid 
    
    # sample coordinates for the first cell 
    xNew1 = rnorm(1, mean = present_cell_x, sd = 1) 
    yNew1 = rnorm(1, mean = present_cell_y, sd = 1) 
    grid[xNew1, yNew1] = 1
    cNew1 = c(xNew1, yNew1) # save coordinates 
    
    # sample coordinates for the second cell 
    gridX = setdiff(1:gridX, xNew1) # exclude previously chosen X
    gridY = setdiff(1:gridY, yNew1) # exclude previously chosen Y
    
    xNew2 = rnorm(1, mean = present_cell_x, sd = 5) # I guess one can tune SD so cells are more or less clustered
    yNew2 = rnorm(1, mean = present_cell_y, sd = 5) # I guess one can tune SD so cells are more or less clustered
    grid[xNew2, yNew2] = 1 
    cNew2 = c(xNew2, yNew2)
    
  } 
  
  nr_cells_added = sum(grid==1) - nr_cells
  expected_nr_cells_added = nr_cells * 2
  
  if (nr_cells_added == expected_nr_cells_added){
    print('Added the correct number of cells')
  }
  else{
    print('Something weird happened when adding cells')
  }
  
  return(grid)
  
}

start_df2 = divide_cells(start_df1)
drawArea(start_df2) # okay so this is just showing you have an empty grid for now 

start_df3 = divide_cells(start_df2)
drawArea(start_df3) # okay so this is just showing you have an empty grid for now 

start_df4 = divide_cells(start_df3)
drawArea(start_df4) # okay so this is just showing you have an empty grid for now 

# select cells from existing ones (select to ICM)
select_to_icm = function(n, grid){ # select n cells from all existing cells 
  
  # available coordinates 
  gridX = dim(grid)[1]
  gridY = dim(grid)[2]
  
  # count the nr of cells on the grid
  nr_cells = sum(grid==1) 
  print(paste0('Number of cells available to go to the ICM:', nr_cells))
  
  # obtain coordinates of all cells which are present 
  coords_available = which(grid==1, arr.ind=TRUE)
  
  coords_to_icm = sample(nrow(coords_available), size = n) # sample n random rows from the matrix (= select of n cells to icm)
  coords_icm = coords_available[c(coords_to_icm),]
  
  print(coords_to_icm)
  print(coords_icm)
  
  # change the value of those cells to 2 
  # will allow to show which cells got selected
  grid[coords_icm] = 2
  return(grid)
}

start_df5 = select_to_icm(3, start_df4) # select 3 cells 
drawArea(start_df4)  
drawArea(start_df5)

divide_icm_cells = function(grid, s){
  
  nr_icm_cells = sum(grid==2) # count nr of icm cells
  
  for (i in 1:nr_icm_cells){
    
    # coordinates of the parent cell 
    present_cell_x = which(grid==1, arr.ind=TRUE)[i, 1]
    present_cell_y = which(grid==1, arr.ind=TRUE)[i, 2]
    
    # we don't care if the cell is on its parent or not 
    # gridX = setdiff(1:gridX, present_cell_x)
    # gridY = setdiff(1:gridY, present_cell_x)
    
    # sample coordinates for new cells 
    # the idea is to sample around the same row and around the same column
    # this creates a spatial structure on the grid 
    
    # sample coordinates for the first cell 
    xNew1 = rnorm(1, mean = present_cell_x, sd = 1) 
    yNew1 = rnorm(1, mean = present_cell_y, sd = 1) 
    grid[xNew1, yNew1] = 3 # add new cells as 3 to be able to tell on the plot
    cNew1 = c(xNew1, yNew1) # save coordinates 
    
    # sample coordinates for the second cell 
    gridX = setdiff(1:gridX, xNew1) # exclude previously chosen X
    gridY = setdiff(1:gridY, yNew1) # exclude previously chosen Y
    
    xNew2 = rnorm(1, mean = present_cell_x, sd = 5) # I guess one can tune SD so cells are more or less clustered
    yNew2 = rnorm(1, mean = present_cell_y, sd = 5) # I guess one can tune SD so cells are more or less clustered
    grid[xNew2, yNew2] = 3 # add new cells as 3 to be able to tell on the plot  
    cNew2 = c(xNew2, yNew2)
    
  } 
  
  nr_cells_added = sum(grid==1) - nr_icm_cells
  expected_nr_cells_added = nr_icm_cells * 2
  
  if (nr_cells_added == expected_nr_cells_added){
    print('Added the correct number of cells')
  }
  else{
    print('Something weird happened when adding cells')
  }
  
  return(grid)
  
}
  
start_df6 = divide_icm_cells(start_df5, 1) # select 3 cells 
drawArea(start_df6)  

# I guess you'd like to ensure that by the end of the cells dividing you have enough cells to split the grid and it's ~full
# so you can just draw a line through the grid or something of that kind 

# split the grid (cells to twin1 vs to twin2 )
split_twins = function(p){ # specify proportion of cells that should go to twin1 
  
}


# plot the output of the simulation
drawArea = function(area_df){  
  
  df = reshape2::melt(as.matrix(area_df, nrow=30)) 
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
    scale_colour_manual(values = c("black","red", "blue", "purple")) # black if absent, red if present 
  return(p)
}



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











