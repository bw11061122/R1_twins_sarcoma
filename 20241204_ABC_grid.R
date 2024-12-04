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

# Define necessary functions

create_empty_df = function(){ # create empty dataframe with 0s to populate it 
  area_df <- data.frame( matrix(0, nrow=gridX, ncol=gridY))
  names(area_df) <- c(1:gridX)  
  return(area_df)
}

start_df = create_empty_df()

# start simulation by randomly placing a cell on the grid 
start_sim = function(){ # get new coordinates by sampling x and y coordinates at random 
    xNew = sample(1:gridX, 1) # get new coordinates (x direction)
    yNew = sample(1:gridY, 1) # get new coordinates (y direction) 
    return( c(xNew,yNew))
}

# divide all cells on the grid (create 2x that many cells and place them in a spatially-coherent manner)
divide_cells = function(){
  
  
  
}

# select cells from existing ones (select to ICM)
select_to_icm = function(n, cells){ # select n cells from all existing cells 
  
  
}


# split the grid (cells to twin1 vs to twin2 )
split_twins = function(p){ # specify proportion of cells that should go to twin1 
  
}


# plot the output of the simulation
drawArea <- function(area_df){  
  
  df <- melt(as.matrix(area_df, nrow=30))  
  brk = c(-1, 0, 1, 2, 24, 1000)
  df$valBucket =cut(df$value, breaks=brk) 
  
  square=15
  p <-NULL
  p <- ggplot(df, aes(X1,X2, color=valBucket)) + geom_point(shape=square, size=6)
  p <- p + scale_colour_manual(values = c("black","blue","orange","yellow","red", "white"))
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











