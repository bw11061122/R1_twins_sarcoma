####################################################################################################################################
# SCRIPT 5

# Getting started with ABC simulations
# December 2024
# Barbara Walkowiak bw18

# Script to run simulations (Henry wants to do this on a grid so I guess this is what we are going to try)
# I found some useful code here: https://gist.github.com/Ram-N/5019098 

###################################################################################################################################
# LIBRARIES
library(data.table)
library(ggplot2)
library(plyr)
library(abc)
library(dplyr)
library(tidyverse)
library(glue)

###################################################################################################################################
# Set working directory 
setwd('/Users/bw18/Desktop/1SB')

######################################################################################################
# PLOT SETTINGS

# Specify colors for plotting 
col_background = "black"
col_TE = "#c81048"
col_PD62341 = "#0ac368"
col_PD63383 = "#a249e8"

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
# TBD, we don't know (number of mutations and their VAFs for now)

###################################################################################################################################
# SIMULATION: ADDING MUTATIONS (IDs = RANDOM STRINGS)

# create ~20k strings that can be used as mutation IDs to assign to cells
# generate all possible 2-letter combination of letters
muts_1id = letters 
muts_2id = c()
for (i in 1:length(letters)){
  for (j in 1:length(letters)){
    muts_2id = c(muts_2id, paste0(letters[i], letters[j])) 
  }
}
muts_3id = c()
for (i in 1:length(letters)){
  for (j in 1:length(letters)){
    for (z in 1:length(letters)){
      muts_3id = c(muts_3id, paste0(letters[i], letters[j], letters[z]))  
    }
  }
}
muts_ids = c(muts_1id, muts_2id, muts_3id)
paste('Number of mutations available to sample from:', length(muts_ids)) # 18,278 so definitely enough

###################################################################################################################################
# SIMULATION: FUNCTIONS

# Define necessary functions

# visualize grid 
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
    coord_equal(ratio = 1)+
    scale_colour_manual(labels = c("no cells", "TE cells", "ICM cells (twin 2)", "ICM cells (twin 1)"),
                        values = c(col_background, col_TE, col_PD63383, col_PD62341, "#dba20a")) # black if absent, red if present 
  return(p)
}


# Helper function to check indices of selected cells
check_indices = function(idx) {
  sorted_indices = sort(idx)
  for (i in 1:(length(idx) - 1)) {
    diff = sorted_indices[i + 1] - sorted_indices[i]
    if (diff == 1 && sorted_indices[i] %% 2 == 1) {
      return(FALSE) # indices differ by 1 and the first is odd, means that you sampled two cells from the same parent
    }
  }
  return(TRUE) # fine, can accept 
}

# Function to simulate twinning 
simulate_twinning = function(x = 100, y = 100, mu1 = 1, mu2 = 1, sd1 = 10, sd2 = 3, s1 = 4, n = 3, s2 = 4, p = 0.5, center = TRUE, cell_mixing = FALSE){
  
  # start the simulation
  grid = data.frame( matrix(0, nrow=x, ncol=y))
  names(grid) = c(1:x)  
  
  # initialize the starting cell in the middle of the grid 
  if (center == TRUE){
    xNew = x/2
    yNew = y/2
    grid[xNew, yNew] = 1 # 1 indicates there is a cell there 
  }
  
  else{
    # initialize the starting cell wherever
    xNew = sample(1:x, 1) # get new coordinates (x direction)
    yNew = sample(1:y, 1) # get new coordinates (y direction) 
    grid[xNew, yNew] = 1 # 1 indicates there is a cell there 
  }
  
  # add mutation to the first cell (1 mutation so we can keep track of who's whose parent, will be present in all cells)
  
  # select mutations from the list
  muts_acquired = sample(muts_ids, 1)
  
  # update muts_ids to avoid sampling the same mutation twice
  muts_ids = setdiff(muts_ids, muts_acquired)
  
  # add mutations to the dt 
  muts_row = t(c(xNew, yNew, 0, muts_acquired)) # first row = 1st cell in the embryo 
  muts_sim = rbind(muts_sim, muts_row, use.names=FALSE)
  
  # simulate s cell divisions 
  
  for (iter in 1:s1){
    
    # count the nr of cells on the grid
    nr_cells = sum(grid==1) 
    
    # find grid dimensions
    x = dim(grid)[1]
    y = dim(grid)[2]
    
    coordParents = which(grid==1, arr.ind=TRUE)
    
    # for each cell, select coordinates of the new daughter cell 
    for (i in 1:nr_cells){
      
      # coordinates of the parent cell 
      xParent = coordParents[i, 1]
      yParent = coordParents[i, 2]
      
      # identify mutations present in the parent cell
      muts_parent = muts_sim[cell_x == xParent & cell_y == yParent, cell_muts] %>% unlist()
      
      # sample coordinates for new cells 
      # the idea is to sample around the same row and around the same column
      # this creates a spatial structure on the grid 
      # ensure that cells do not overlap 
      
      # sample coordinates for the first cell
      repeat{ # repeat until coordinate sampled 
        xNew1 = round(rnorm(1, mean = xParent, sd = sd1))
        yNew1 = round(rnorm(1, mean = yParent, sd = sd1))
        
        if (xNew1 >= 1 & xNew1 <= x & yNew1 >= 1 & yNew1 <= y){ # check coordinates not outside the grid
          if (grid[xNew1, yNew1] == 0){ # add the cell if the position is not occupied
            grid[xNew1, yNew1] = 1 
            
            # add mutations
            # choose nr of mutations according to the mutation rate
            if (s1 <= 2){
              nr_muts = rpois(lambda=mu1, n=1)
              } # before ZGA (until 4-cell stage), use mu1
            else{
              nr_muts = rpois(lambda=mu2, n=1)
              } # after ZGA (past 4-cell stage), use mu2

            # select mutations from the list
            muts_acquired1 = sample(muts_ids, nr_muts)
            muts_all1 = c(muts_parent, muts_acquired1)
            muts_all1 = paste(muts_all1, collapse = ',')
            
            # update muts_ids to avoid sampling the same mutation twice
            muts_ids = setdiff(muts_ids, muts_all1)
            
            # add mutations to the dt 
            muts_row = t(c(xNew1, yNew1, iter, muts_all1))
            muts_sim = rbind(muts_sim, muts_row, use.names=FALSE)
            
            break 
          } 
        }
      }
      
      # sample coordinates for the second cell
      repeat{ # repeat until coordinate sampled 
        xNew2 = round(rnorm(1, mean = xParent, sd = sd1))
        yNew2 = round(rnorm(1, mean = yParent, sd = sd1))
        
        if (xNew2 >= 1 & xNew2 <= x & yNew2 >= 1 & yNew2 <= y){ # check coordinates not outside the grid 
          if (grid[xNew2, yNew2] == 0){ # add the cell if the position is not occupied 
            grid[xNew2, yNew2] = 1 
            
            # add mutations
            # choose nr of mutations according to the mutation rate
            if (s1 <= 2){
              nr_muts = rpois(lambda=mu1, n=1)} # before ZGA (until 4-cell stage), use mu1
            else{
              nr_muts = rpois(lambda=mu2, n=1)} # after ZGA (past 4-cell stage), use mu2
            
            # select mutations from the list
            muts_acquired2 = sample(muts_ids, nr_muts)
            muts_all2 = c(muts_parent, muts_acquired2)
            muts_all2 = paste(muts_all2, collapse = ',')
            
            # update muts_ids to avoid sampling the same mutation twice
            muts_ids = setdiff(muts_ids, muts_all2)
            
            # add mutations to the dt 
            muts_row = t(c(xNew2, yNew2, iter, muts_all2))
            muts_sim = rbind(muts_sim, muts_row, use.names=FALSE)
            
            break 
          }
        }
      }
      
      # remove the mother cell from the grid 
      grid[xParent, yParent] = 0
      
    } 
    
    nr_cells_added = sum(grid==1) 
    expected_nr_cells_added = nr_cells * 2
  
  }
  
  # select coordinates of cells existing after the last cell division
  last_gen = s1
  coords_cells = muts_sim[cell_gen==last_gen, c('cell_x', 'cell_y'), with=FALSE]
  coords_cells[, names(coords_cells) := lapply(.SD, as.numeric), .SDcols = names(coords_cells)]
  
  # cells ingressed through asymmetric divisions contribute to the ICM
  # therefore, allow only one of the two cells generated from each parent to be allocated to the ICM
  
  # check the ID of the first row of the generation you are sampling from
  repeat{
    indices_to_icm = sample(nrow(coords_cells), size = n) # sample n random row indexes from the matrix (= select of n cells to icm)
    
    # check that you don't have 1 and 2; 3 and 4; 5 and 6 and so on 
    if (check_indices(indices_to_icm)) {
      accepted_indices = indices_to_icm
      break 
    }
  }
  
  coords_icm = coords_cells[accepted_indices,c('cell_x', 'cell_y'), with=FALSE] # select coordinates for ICM cells

  # add column in muts_sim to indicate which cells went to the ICM
  # calculate the indices in the actual dt
  sum_cells_before = dim(muts_sim[cell_gen<=s1-1])[1]
  accepted_indices_2 = accepted_indices + sum_cells_before
  muts_sim[, makes_ICM := ifelse(.I %in% accepted_indices_2, 1, 0)]
  
  # change value of ICM cells to 2
  for (i in 1:nrow(coords_icm)){
   x_icm = as.numeric(coords_icm[i,1])
   y_icm = as.numeric(coords_icm[i,2])
   grid[x_icm, y_icm] = 2
  }

  cells_icm = sum(grid==2)
  
  # now divide ICM cells
  for (iter2 in 1:s2){ # repeat cell division s2 times 
    
    coordParents = which(grid==2, arr.ind=TRUE)
    nr_cells_divide = nrow(coordParents)
    
    for (c in 1:nr_cells_divide){ # for each cell in the ICM 
      
      # coordinates of the parent cell 
      xParent = coordParents[c, 1]
      yParent = coordParents[c, 2]
      
      # identify mutations present in the parent cell
      muts_parent = muts_sim[cell_x == xParent & cell_y == yParent, cell_muts] %>% unlist()
      
      # sample coordinates for new cells 
      # the idea is to sample around the same row and around the same column
      # this creates a spatial structure on the grid 
      # ensure that cells do not overlap 
      
      # sample coordinates for the first cell
      repeat{ # repeat until coordinate sampled 
        xNew1 = round(rnorm(1, mean = xParent, sd = sd2))
        yNew1 = round(rnorm(1, mean = yParent, sd = sd2))
        
        if (xNew1 >= 1 & xNew1 <= x & yNew1 >= 1 & yNew1 <= y){ # check coordinates not outside the grid
          if (grid[xNew1, yNew1] == 0){ # add the cell if the position is not occupied
            grid[xNew1, yNew1] = 2 
            
            # add mutations
            # choose nr of mutations according to the mutation rate
            nr_muts = rpois(lambda=mu2, n=1) # only return one number 
            
            # select mutations from the list
            muts_acquired1 = sample(muts_ids, nr_muts)
            muts_all1 = c(muts_parent, muts_acquired1)
            muts_all1 = paste(muts_all1, collapse = ',')
            
            # update muts_ids to avoid sampling the same mutation twice
            muts_ids = setdiff(muts_ids, muts_all1)
            # 2 to label that this is a daughter cell of one of the ICM cells
            
            # add mutations to the dt 
            muts_row = t(c(xNew1, yNew1, iter2+s1, muts_all1, 2))
            muts_sim = rbind(muts_sim, muts_row, use.names=FALSE)
            
            break 
          } 
        }
      }
      
      # sample coordinates for the second cell
      repeat{ # repeat until coordinate sampled 
        xNew2 = round(rnorm(1, mean = xParent, sd = sd2))
        yNew2 = round(rnorm(1, mean = yParent, sd = sd2))
        
        if (xNew2 >= 1 & xNew2 <= x & yNew2 >= 1 & yNew2 <= y){ # check coordinates not outside the grid 
          if (grid[xNew2, yNew2] == 0){ # add the cell if the position is not occupied 
            grid[xNew2, yNew2] = 2 
            
            # add mutations
            # choose nr of mutations according to the mutation rate
            nr_muts = rpois(lambda=mu2, n=1) # only return one number 
            
            # select mutations from the list
            muts_acquired2 = sample(muts_ids, nr_muts)
            muts_all2 = c(muts_parent, muts_acquired2)
            muts_all2 = paste(muts_all2, collapse = ',')
            
            # update muts_ids to avoid sampling the same mutation twice
            muts_ids = setdiff(muts_ids, muts_all2)
            
            # add mutations to the dt 
            muts_row = t(c(xNew2, yNew2, iter2+s1, muts_all2, 2))
            muts_sim = rbind(muts_sim, muts_row, use.names=FALSE)
            
            break 
          }
        }
      }
      
      # remove the mother cell from the grid 
      grid[xParent, yParent] = 0
      
    } 
    
  }
 
  # finished all ICM divisions until embryo split
  # select cells to twin 1 and twin 2
  
  # determine nr of cells to be allocated to twin1
  nr_cells_twins = sum(grid==2)
  nr_cells_twin1 = round(p * nr_cells_twins)
  
  # if cell mixing not allowed, bisect the grid 
  if (cell_mixing == FALSE){
    
    sum_cells = 0
  
    for (i in 1:nrow(grid)){
      cells_in_row = which(grid[i,]==2)
      if (sum_cells <= nr_cells_twin1){
        if (length(cells_in_row > 0)){
          for (j in cells_in_row){
            grid[i, j] = 3 # label twin 1 cells with a different number
            sum_cells = sum_cells + 1
          } 
        }
      }
    }
  }
  
  # if cell mixing is allowed, select cells at random
  else {

    coords_avail = which(grid==2, arr.ind = TRUE)
    coords_twin1 = coords_avail[sample(nrow(coords_avail), size = nr_cells_twin1),]
    for (i in 1:nrow(coords_twin1)){
        a = as.numeric(coords_twin1[i, 1])
        b = as.numeric(coords_twin1[i, 2])
        grid[a, b] = 3
      }
    }

  nr_cells_twin1 = sum(grid==3)
  nr_cells_twin2 = sum(grid==2)

  # add column to specify if cell allocated to twin1 or twin2
  coords_twin1 = data.frame(which(grid==3, arr.ind=TRUE))
  coords_twin2 = data.frame(which(grid==2, arr.ind=TRUE))
  coords_twin1$coord = paste0(coords_twin1$row, '_', coords_twin1$col)
  ctwin1 = coords_twin1$coord %>% unlist()
  coords_twin2$coord = paste0(coords_twin2$row, '_', coords_twin2$col)
  ctwin2 = coords_twin2$coord %>% unlist()
  
  muts_sim[, coord := paste0(cell_x, '_', cell_y)]
  muts_sim[, twin := factor(fcase(
    coord %in% ctwin1, 'twin 1',
    coord %in% ctwin2, 'twin 2',
    !coord %in% c(ctwin1, ctwin2) & makes_ICM == '1', 'ICM',
    !coord %in% c(ctwin1, ctwin2) & makes_ICM == '0' & cell_gen==s1, 'TE',
    !coord %in% c(ctwin1, ctwin2) & makes_ICM == '0' & cell_gen!=s1, 'divided',
    !coord %in% c(ctwin1, ctwin2) & cell_gen>s1, 'divided'
  ))]
  
  # before you return muts_sim dataframe, I want to also add info on which parameters this was obtained with
  muts_sim[, mu1 := mu1]
  muts_sim[, mu2 := mu2]
  muts_sim[, s1 := s1]
  muts_sim[, s2 := s2]
  muts_sim[, sd1 := sd1]
  muts_sim[, sd2 := sd2]
  muts_sim[, n := n]
  muts_sim[, p := p]

  out = list(grid, muts_sim, muts_ids)
  return(out)
  
}

# write a function to extract interesting output from the simulation
get_output = function(muts_sim){  
  
  nr_cells_twin1 = dim(muts_sim[twin == 'twin 1'])[1]
  nr_cells_twin2 = dim(muts_sim[twin == 'twin 2'])[1]

  # extract IDs of mutations present in cells from either twin 
  muts_twin1 = muts_sim[twin == 'twin 1', cell_muts] %>% unlist()
  muts_twin2 = muts_sim[twin == 'twin 2', cell_muts] %>% unlist()
  muts_twin1 = unlist(lapply(muts_twin1, function(x) strsplit(x, split = ",")))
  muts_twin2 = unlist(lapply(muts_twin2, function(x) strsplit(x, split = ",")))
  
  # count mutations present in each twin
  counts_muts_twin1 = data.table(table(muts_twin1))
  counts_muts_twin2 = data.table(table(muts_twin2))
  setnames(counts_muts_twin1, c('muts_twin1', 'N'), c('mut_ID', 'N1'))
  setnames(counts_muts_twin2, c('muts_twin2', 'N'), c('mut_ID', 'N2'))
  
  # add info on mutation VAF and fraction of cells where mutation is present 
  counts_muts_twin1[, ccf1 := N1/nr_cells_twin1]
  counts_muts_twin2[, ccf2 := N2/nr_cells_twin2]
  counts_muts_twin1[, vaf1 := ccf1/2]
  counts_muts_twin2[, vaf2 := ccf2/2]
  
  # create a shared df with muts in either twin 
  counts_muts_both = merge(counts_muts_twin1, counts_muts_twin2, by = 'mut_ID', all = TRUE)
  counts_muts_both[is.na(counts_muts_both)] = 0 # replace NA with 0
  counts_muts_both[, twin1_specific := as.numeric(vaf1 >= 0.1 & vaf2 == 0)]
  counts_muts_both[, twin2_specific := as.numeric(vaf1 == 0 & vaf2 >= 0.1)]
  counts_muts_both[, shared := as.numeric(vaf1 >= 0.1 & vaf2 >= 0.1)]
  
  # extract parameters to write to new dt
  mu1 = muts_sim[, mu1] %>% unlist() %>% unique()
  mu2 = muts_sim[, mu2] %>% unlist() %>% unique()
  sd1 = muts_sim[, sd1] %>% unlist() %>% unique()
  sd2 = muts_sim[, sd2] %>% unlist() %>% unique()
  s1 = muts_sim[, s1] %>% unlist() %>% unique()
  s2 = muts_sim[, s2] %>% unlist() %>% unique()
  n = muts_sim[, n] %>% unlist() %>% unique()
  p = muts_sim[, p] %>% unlist() %>% unique()

  # extract summary stats 
  nr_shared = dim(counts_muts_both[shared==1])[1]
  nr_twin1 = dim(counts_muts_both[twin1_specific==1])[1]
  nr_twin2 = dim(counts_muts_both[twin2_specific==1])[1]
  
  # extract vafs for each category of mutation
  vafs_shared_twin1 = paste(counts_muts_both[shared==1, vaf1] %>% unlist(), collapse = ',')
  vafs_shared_twin2 = paste(counts_muts_both[shared==1, vaf2] %>% unlist(), collapse = ',')
  vafs_spec_twin1 = paste(counts_muts_both[twin1_specific==1, vaf1] %>% unlist(), collapse = ',')
  vafs_spec_twin2 = paste(counts_muts_both[twin2_specific==1, vaf2] %>% unlist(), collapse = ',')
  
  out_row = t(c(mu1, mu2, sd1, sd2, s1, s2, n, p,
                nr_shared, nr_twin1, nr_twin2, vafs_shared_twin1, vafs_shared_twin2, vafs_spec_twin1, vafs_spec_twin2))
  out_sim = rbind(out_sim, out_row, use.names = FALSE)
  return(out_sim)
  
}

# Simulate sequencing

# Convert the VAF character to a list
char_to_list = function(vafs, lambda = 30){
  vafs = vafs %>% strsplit(',') %>% unlist() %>% as.numeric() # create a list of VAF values from the character 
  return(vafs)
}

# Take each VAF stored in the column, sample reads from Poisson w/ lambda 30, determine sampled VAF
sim_sequencing = function(vafs, lambda = 30){
  vafs = vafs %>% strsplit(',') %>% unlist() %>% as.numeric() # create a list of VAF values from the character 
  if (length(vafs)==0){
    return(0)
  }
  else {
    vafs_seq = c()
    for (vaf in vafs){ 
      depth = rpois(1, lambda) # sample nr of all reads sequenced from Poisson distribution, lambda = 30 
      mtr = rbinom(1, size = depth, prob = vaf) # sample nr of mutant reads given coverage and true VAF 
      vaf_seq = mtr / depth # calculate the VAF that would be estimated through sequencing 
      vafs_seq = c(vafs_seq, vaf_seq)
    } 
  return(vafs_seq)
  }
}

# calculate difference b/n simulated and observed VAFs (summary statistics) 
get_vaf_diff = function(my_list, sum_observed_vafs) {
  sim_vafs = sum(my_list) # sum simulated (sequenced) vafs
  result = sim_vafs - sum_observed_vafs  
  return(result)
}

###################################################################################################################################
# SIMULATION: run with 1 varying parameter (s2)

# what are the reasonable parameter values?
# x, y, sd1, sd2 = for x and y, anything big will do
# sd1, sd2 = does Henry have any intuition for this? I was thinking sd1 > sd2, but not sure otherwise
# s1 = 3-5?
# s2 = 2-6?
# n = 2-6
# p = 0.25, 0.5, 0.75 (can try more fine-grained but realistically 0.2-0.8)
# mu1, mu2 = 2 for before ZGA and 0.5 for after ZGA

# set fixed parameters 
x = 100
y = 100
sd1 = 10
sd2 = 3
s1 = 3
n = 3
p = 0.5
mu1 = 2
mu2 = 0.5

# observed statistics 
nr_shared_observed = 1
nr_spec_twin1_observed = 6
nr_spec_twin2_observed = 12
vafs_shared_twin1_observed = 0.2
vafs_shared_twin2_observed = 0.5
vafs_spec_twin1_observed = c(0.48, 0.4, 0.4, 0.4, 0.4, 0.2)
vafs_spec_twin2_observed = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
sum_vafs_spec_twin1_observed = sum(vafs_spec_twin1_observed)
sum_vafs_spec_twin2_observed = sum(vafs_spec_twin2_observed)

# set seed
seed = (runif(1, min=0, max=100) + runif(1, min=50, max=100)) * runif(1, min=0, max=100)
print(paste0("Random seed used in the simulation: ", seed))
set.seed(seed)

# create an empty dt to store mutations for each cell on the phylogeny
muts_sim=data.table(cell_x=numeric(), cell_y=numeric(), cell_gen=numeric(), cell_muts=character())
# create an empty dt to store simulation results in 
out_sim=data.table(mu1=numeric(), mu2=numeric(), sd1=numeric(), sd2=numeric(), s1=numeric(), s2=numeric(), n=numeric(), p=numeric(),
                   nr_shared=numeric(), nr_twin1=numeric(), nr_twin2 = numeric(),
                   vafs_shared_twin1=character(), vafs_shared_twin2=character(), vafs_spec_twin1=character(), vafs_spec_twin2=character())

# test out from s2 in 2 to 6 (2-6 ICM cell divisions)
for (s2 in 2:6){
  
  # run 100 simulations for each parameter value 
  for (rep in 1:1000){
  
    if (rep == 1000){
      print(rep)
    }
    
    out = simulate_twinning(x, y, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, n = n, s1 = s1, s2 = s2, p = p)
    
    out_muts = out[[2]]
    out_sim = get_output(out_muts)
    
    out_grid = out[[1]]
    if (rep == 100){
      draw_area(out_grid)
      ggsave(glue('FiguresAdd/F5/F5_sim_output_mu1_{mu1}_mu2_{mu2}_sd1_{sd1}_sd2_{sd2}_s1_{s1}_s2_{s2}_n_{n}_p_{p}.pdf'), width = 4, height = 4)
    }
    
  }
}

out_sim_s2_bisect = data.table(out_sim)
out_sim_s2_bisect[, mixing := FALSE]

# make sure that relevant columns in out_sim are numeric
num_cols = c('mu1', 'mu2', 'sd1', 'sd2', 's1', 's2', 'n', 'p', 'nr_shared', 'nr_twin1', 'nr_twin2')
out_sim_s2_bisect[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

# ADD SEQUENCING SIMULATION TO DEPTHS 30x (1 SAMPLE) AND 200x (AGGREGATED)

# add VAFs observed "after sequencing" to 30x coverage (1 WGS sample)
out_sim_s2_bisect[, vafs_shared_twin1_seq30 := lapply(vafs_shared_twin1, sim_sequencing, lambda = 30)] # vafs of shared mutations in twin1 
out_sim_s2_bisect[, vafs_shared_twin2_seq30 := lapply(vafs_shared_twin2, sim_sequencing, lambda = 30)] # vafs of shared mutations in twin1
out_sim_s2_bisect[, vafs_spec_twin1_seq30 := lapply(vafs_spec_twin1, sim_sequencing, lambda = 30)]
out_sim_s2_bisect[, vafs_spec_twin2_seq30 := lapply(vafs_spec_twin2, sim_sequencing, lambda = 30)]
# add how many mutations would be observed based on the sequenced VAF (from sequencing to 30x)
out_sim_s2_bisect[, nr_shared_twin1_seq30 := sapply(vafs_shared_twin1_seq30, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_s2_bisect[, nr_shared_twin2_seq30 := sapply(vafs_shared_twin2_seq30, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_s2_bisect[, nr_spec_twin1_seq30 := sapply(vafs_spec_twin1_seq30, function(x) length(x[x > 0.1]))]
out_sim_s2_bisect[, nr_spec_twin2_seq30 := sapply(vafs_spec_twin2_seq30, function(x) length(x[x > 0.1]))]

# add VAFs observed "after sequencing" to 200x coverage (all WGS samples aggregated)
out_sim_s2_bisect[, vafs_shared_twin1_seq200 := lapply(vafs_shared_twin1, sim_sequencing, lambda = 200)] # vafs of shared mutations in twin1 
out_sim_s2_bisect[, vafs_shared_twin2_seq200 := lapply(vafs_shared_twin2, sim_sequencing, lambda = 200)] # vafs of shared mutations in twin1
out_sim_s2_bisect[, vafs_spec_twin1_seq200 := lapply(vafs_spec_twin1, sim_sequencing, lambda = 200)]
out_sim_s2_bisect[, vafs_spec_twin2_seq200 := lapply(vafs_spec_twin2, sim_sequencing, lambda = 200)]
# add how many mutations would be observed based on the sequenced VAF 
out_sim_s2_bisect[, nr_shared_twin1_seq200 := sapply(vafs_shared_twin1_seq200, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_s2_bisect[, nr_shared_twin2_seq200 := sapply(vafs_shared_twin2_seq200, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_s2_bisect[, nr_spec_twin1_seq200 := sapply(vafs_spec_twin1_seq200, function(x) length(x[x > 0.1]))]
out_sim_s2_bisect[, nr_spec_twin2_seq200 := sapply(vafs_spec_twin2_seq200, function(x) length(x[x > 0.1]))]

# CALCULATE DIFFERENCE TO SUMMARY STATISTICS 
# determine difference between observed and simulated data (raw - w/o "sequencing")
out_sim_s2_bisect[, nr_shared_diff_raw := nr_shared - nr_shared_observed]
out_sim_s2_bisect[, nr_spec_twin1_diff_raw := nr_twin1 - nr_spec_twin1_observed]
out_sim_s2_bisect[, nr_spec_twin2_diff_raw := nr_twin2 - nr_spec_twin2_observed]
out_sim_s2_bisect[, vafs_shared_twin1 := lapply(vafs_shared_twin1, char_to_list)]
out_sim_s2_bisect[, vafs_shared_twin2 := lapply(vafs_shared_twin2, char_to_list)]
out_sim_s2_bisect[, vafs_spec_twin1 := lapply(vafs_spec_twin1, char_to_list)]
out_sim_s2_bisect[, vafs_spec_twin2 := lapply(vafs_spec_twin2, char_to_list)]
out_sim_s2_bisect[, vaf_shared_twin1_diff_raw := mapply(get_vaf_diff, vafs_shared_twin1, vafs_shared_twin1_observed)]
out_sim_s2_bisect[, vaf_shared_twin2_diff_raw := mapply(get_vaf_diff, vafs_shared_twin2, vafs_shared_twin2_observed)]
out_sim_s2_bisect[, vaf_spec_twin1_diff_raw := mapply(get_vaf_diff, vafs_spec_twin1, sum_vafs_spec_twin1_observed)]
out_sim_s2_bisect[, vaf_spec_twin2_diff_raw := mapply(get_vaf_diff, vafs_spec_twin2, sum_vafs_spec_twin2_observed)]

# difference between simulated (sequenced) and observed statistics (sequencing to 30x)
out_sim_s2_bisect[, nr_shared_twin1_diff_seq30 := nr_shared_twin1_seq30 - nr_shared_observed]
out_sim_s2_bisect[, nr_shared_twin2_diff_seq30 := nr_shared_twin2_seq30 - nr_shared_observed]
out_sim_s2_bisect[, nr_spec_twin1_diff_seq30 := nr_spec_twin1_seq30 - nr_spec_twin1_observed]
out_sim_s2_bisect[, nr_spec_twin2_diff_seq30 := nr_spec_twin2_seq30 - nr_spec_twin2_observed]
out_sim_s2_bisect[, vaf_shared_twin1_diff_seq30 := mapply(get_vaf_diff, vafs_shared_twin1_seq30, vafs_shared_twin1_observed)]
out_sim_s2_bisect[, vaf_shared_twin2_diff_seq30 := mapply(get_vaf_diff, vafs_shared_twin2_seq30, vafs_shared_twin2_observed)]
out_sim_s2_bisect[, vaf_spec_twin1_diff_seq30 := mapply(get_vaf_diff, vafs_spec_twin1_seq30, sum_vafs_spec_twin1_observed)]
out_sim_s2_bisect[, vaf_spec_twin2_diff_seq30 := mapply(get_vaf_diff, vafs_spec_twin2_seq30, sum_vafs_spec_twin2_observed)]

# difference between simulated (sequenced) and observed statistics (sequencing to 200x)
out_sim_s2_bisect[, nr_shared_twin1_diff_seq200 := nr_shared_twin1_seq200 - nr_shared_observed]
out_sim_s2_bisect[, nr_shared_twin2_diff_seq200 := nr_shared_twin2_seq200 - nr_shared_observed]
out_sim_s2_bisect[, nr_spec_twin1_diff_seq200 := nr_spec_twin1_seq200 - nr_spec_twin1_observed]
out_sim_s2_bisect[, nr_spec_twin2_diff_seq200 := nr_spec_twin2_seq200 - nr_spec_twin2_observed]
out_sim_s2_bisect[, vaf_shared_twin1_diff_seq200 := mapply(get_vaf_diff, vafs_shared_twin1_seq200, vafs_shared_twin1_observed)]
out_sim_s2_bisect[, vaf_shared_twin2_diff_seq200 := mapply(get_vaf_diff, vafs_shared_twin2_seq200, vafs_shared_twin2_observed)]
out_sim_s2_bisect[, vaf_spec_twin1_diff_seq200 := mapply(get_vaf_diff, vafs_spec_twin1_seq200, sum_vafs_spec_twin1_observed)]
out_sim_s2_bisect[, vaf_spec_twin2_diff_seq200 := mapply(get_vaf_diff, vafs_spec_twin2_seq200, sum_vafs_spec_twin2_observed)]

# repeat, but now allow for cell mixing 

# set seed
seed = (runif(1, min=0, max=100) + runif(1, min=50, max=100)) * runif(1, min=0, max=100)
print(paste0("Random seed used in the simulation: ", seed))
set.seed(seed)

# create an empty dt to store mutations for each cell on the phylogeny
muts_sim=data.table(cell_x=numeric(), cell_y=numeric(), cell_gen=numeric(), cell_muts=character())
# create an empty dt to store simulation results in 
out_sim=data.table(mu1=numeric(), mu2=numeric(), sd1=numeric(), sd2=numeric(), s1=numeric(), s2=numeric(), n=numeric(), p=numeric(),
                   nr_shared=numeric(), nr_twin1=numeric(), nr_twin2 = numeric(),
                   vafs_shared_twin1=character(), vafs_shared_twin2=character(), vafs_spec_twin1=character(), vafs_spec_twin2=character())

# test out from s2 in 2 to 6 (2-6 ICM cell divisions)
for (s2 in 2:6){
  
  # run 100 simulations for each parameter value 
  for (rep in 1:1000){
    
    if (rep == 1000){
      print(rep)
    }
    out = simulate_twinning(x, y, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, n = n, s1 = s1, s2 = s2, p = p, cell_mixing = TRUE)
    
    out_muts = out[[2]]
    out_sim = get_output(out_muts)
    
    out_grid = out[[1]]
    if (rep == 100){
      draw_area(out_grid)
      ggsave(glue('FiguresAdd/F5/F5_sim_output_mu1_{mu1}_mu2_{mu2}_sd1_{sd1}_sd2_{sd2}_s1_{s1}_s2_{s2}_n_{n}_p_{p}_mixing.pdf'), width = 4, height = 4)
    }
    
  }
}

# make sure that relevant columns in out_sim are numeric
num_cols = c('mu1', 'mu2', 'sd1', 'sd2', 's1', 's2', 'n', 'p', 'nr_shared', 'nr_twin1', 'nr_twin2')
out_sim[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

out_sim_mixing = data.table(out_sim)
out_sim_mixing[, mixing := TRUE]

# add VAFs observed "after sequencing" to 30x coverage (1 WGS sample)
out_sim_mixing[, vafs_shared_twin1_seq30 := lapply(vafs_shared_twin1, sim_sequencing, lambda = 30)] # vafs of shared mutations in twin1 
out_sim_mixing[, vafs_shared_twin2_seq30 := lapply(vafs_shared_twin2, sim_sequencing, lambda = 30)] # vafs of shared mutations in twin1
out_sim_mixing[, vafs_spec_twin1_seq30 := lapply(vafs_spec_twin1, sim_sequencing, lambda = 30)]
out_sim_mixing[, vafs_spec_twin2_seq30 := lapply(vafs_spec_twin2, sim_sequencing, lambda = 30)]
# add how many mutations would be observed based on the sequenced VAF (from sequencing to 30x)
out_sim_mixing[, nr_shared_twin1_seq30 := sapply(vafs_shared_twin1_seq30, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_mixing[, nr_shared_twin2_seq30 := sapply(vafs_shared_twin2_seq30, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_mixing[, nr_spec_twin1_seq30 := sapply(vafs_spec_twin1_seq30, function(x) length(x[x > 0.1]))]
out_sim_mixing[, nr_spec_twin2_seq30 := sapply(vafs_spec_twin2_seq30, function(x) length(x[x > 0.1]))]

# add VAFs observed "after sequencing" to 200x coverage (all WGS samples aggregated)
out_sim_mixing[, vafs_shared_twin1_seq200 := lapply(vafs_shared_twin1, sim_sequencing, lambda = 200)] # vafs of shared mutations in twin1 
out_sim_mixing[, vafs_shared_twin2_seq200 := lapply(vafs_shared_twin2, sim_sequencing, lambda = 200)] # vafs of shared mutations in twin1
out_sim_mixing[, vafs_spec_twin1_seq200 := lapply(vafs_spec_twin1, sim_sequencing, lambda = 200)]
out_sim_mixing[, vafs_spec_twin2_seq200 := lapply(vafs_spec_twin2, sim_sequencing, lambda = 200)]
# add how many mutations would be observed based on the sequenced VAF 
out_sim_mixing[, nr_shared_twin1_seq200 := sapply(vafs_shared_twin1_seq200, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_mixing[, nr_shared_twin2_seq200 := sapply(vafs_shared_twin2_seq200, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_mixing[, nr_spec_twin1_seq200 := sapply(vafs_spec_twin1_seq200, function(x) length(x[x > 0.1]))]
out_sim_mixing[, nr_spec_twin2_seq200 := sapply(vafs_spec_twin2_seq200, function(x) length(x[x > 0.1]))]

# CALCULATE DIFFERENCE TO SUMMARY STATISTICS 
# determine difference between observed and simulated data (raw - w/o "sequencing")
out_sim_mixing[, nr_shared_diff_raw := nr_shared - nr_shared_observed]
out_sim_mixing[, nr_spec_twin1_diff_raw := nr_twin1 - nr_spec_twin1_observed]
out_sim_mixing[, nr_spec_twin2_diff_raw := nr_twin2 - nr_spec_twin2_observed]
out_sim_mixing[, vafs_shared_twin1 := lapply(vafs_shared_twin1, char_to_list)]
out_sim_mixing[, vafs_shared_twin2 := lapply(vafs_shared_twin2, char_to_list)]
out_sim_mixing[, vafs_spec_twin1 := lapply(vafs_spec_twin1, char_to_list)]
out_sim_mixing[, vafs_spec_twin2 := lapply(vafs_spec_twin2, char_to_list)]
out_sim_mixing[, vaf_shared_twin1_diff_raw := mapply(get_vaf_diff, vafs_shared_twin1, vafs_shared_twin1_observed)]
out_sim_mixing[, vaf_shared_twin2_diff_raw := mapply(get_vaf_diff, vafs_shared_twin2, vafs_shared_twin2_observed)]
out_sim_mixing[, vaf_spec_twin1_diff_raw := mapply(get_vaf_diff, vafs_spec_twin1, sum_vafs_spec_twin1_observed)]
out_sim_mixing[, vaf_spec_twin2_diff_raw := mapply(get_vaf_diff, vafs_spec_twin2, sum_vafs_spec_twin2_observed)]

# difference between simulated (sequenced) and observed statistics (sequencing to 30x)
out_sim_mixing[, nr_shared_twin1_diff_seq30 := nr_shared_twin1_seq30 - nr_shared_observed]
out_sim_mixing[, nr_shared_twin2_diff_seq30 := nr_shared_twin2_seq30 - nr_shared_observed]
out_sim_mixing[, nr_spec_twin1_diff_seq30 := nr_spec_twin1_seq30 - nr_spec_twin1_observed]
out_sim_mixing[, nr_spec_twin2_diff_seq30 := nr_spec_twin2_seq30 - nr_spec_twin2_observed]
out_sim_mixing[, vaf_shared_twin1_diff_seq30 := mapply(get_vaf_diff, vafs_shared_twin1_seq30, vafs_shared_twin1_observed)]
out_sim_mixing[, vaf_shared_twin2_diff_seq30 := mapply(get_vaf_diff, vafs_shared_twin2_seq30, vafs_shared_twin2_observed)]
out_sim_mixing[, vaf_spec_twin1_diff_seq30 := mapply(get_vaf_diff, vafs_spec_twin1_seq30, sum_vafs_spec_twin1_observed)]
out_sim_mixing[, vaf_spec_twin2_diff_seq30 := mapply(get_vaf_diff, vafs_spec_twin2_seq30, sum_vafs_spec_twin2_observed)]

# difference between simulated (sequenced) and observed statistics (sequencing to 200x)
out_sim_mixing[, nr_shared_twin1_diff_seq200 := nr_shared_twin1_seq200 - nr_shared_observed]
out_sim_mixing[, nr_shared_twin2_diff_seq200 := nr_shared_twin2_seq200 - nr_shared_observed]
out_sim_mixing[, nr_spec_twin1_diff_seq200 := nr_spec_twin1_seq200 - nr_spec_twin1_observed]
out_sim_mixing[, nr_spec_twin2_diff_seq200 := nr_spec_twin2_seq200 - nr_spec_twin2_observed]
out_sim_mixing[, vaf_shared_twin1_diff_seq200 := mapply(get_vaf_diff, vafs_shared_twin1_seq200, vafs_shared_twin1_observed)]
out_sim_mixing[, vaf_shared_twin2_diff_seq200 := mapply(get_vaf_diff, vafs_shared_twin2_seq200, vafs_shared_twin2_observed)]
out_sim_mixing[, vaf_spec_twin1_diff_seq200 := mapply(get_vaf_diff, vafs_spec_twin1_seq200, sum_vafs_spec_twin1_observed)]
out_sim_mixing[, vaf_spec_twin2_diff_seq200 := mapply(get_vaf_diff, vafs_spec_twin2_seq200, sum_vafs_spec_twin2_observed)]

# SAVE SIMULATION RESULTS
# save simulation results to out
out_sim_s2 = rbind(out_sim_s2_bisect, out_sim_mixing)

# define which summary statistics are of interest
summary_stats_raw = c("nr_shared_diff_raw", "nr_spec_twin1_diff_raw",   
                  "nr_spec_twin2_diff_raw",  "vaf_shared_twin1_diff_raw", "vaf_shared_twin2_diff_raw",
                  "vaf_spec_twin1_diff_raw",   "vaf_spec_twin2_diff_raw")
summary_stats_seq30 = c("nr_shared_twin1_diff_seq30",  "nr_shared_twin2_diff_seq30",  "nr_spec_twin1_diff_seq30",   
                        "nr_spec_twin2_diff_seq30",  "vaf_shared_twin1_diff_seq30", "vaf_shared_twin2_diff_seq30",
                        "vaf_spec_twin1_diff_seq30",   "vaf_spec_twin2_diff_seq30")
summary_stats_seq200 = c("nr_shared_twin1_diff_seq200",  "nr_shared_twin2_diff_seq200",  "nr_spec_twin1_diff_seq200",   
                         "nr_spec_twin2_diff_seq200",  "vaf_shared_twin1_diff_seq200", "vaf_shared_twin2_diff_seq200",
                         "vaf_spec_twin1_diff_seq200",   "vaf_spec_twin2_diff_seq200")

# melt results dt
out_sim_s2_sub_raw = out_sim_s2[, c('s2', summary_stats_raw, 'mixing'), with=FALSE]
out_sim_s2_raw_melt = data.table::melt(out_sim_s2_sub_raw, id.vars = c('s2', 'mixing'))
out_sim_s2_raw_melt[, variable := as.factor(variable)]
out_sim_s2_raw_melt[, s2:= as.factor(s2)]

# compare values of summary statistics for different summary stats
ggplot(out_sim_s2_raw_melt, aes(x = s2, y = value, col = mixing))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 3)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data', col = 'Cell mixing')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cfMixing_uncorrectedSeq.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2_raw_melt[mixing==FALSE], aes(x = s2, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 3)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_noMixing_uncorrectedSeq.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2_raw_melt[mixing==TRUE], aes(x = s2, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 3)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data', col)
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cellMixing_uncorrectedSeq.pdf'), height = 8, width = 8)  

# results (sequencing to 30x coverage)
out_sim_s2_sub_seq30 = out_sim_s2[, c('s2', summary_stats_seq30, 'mixing'), with=FALSE]
out_sim_s2_seq30_melt = data.table::melt(out_sim_s2_sub_seq30, id.vars = c('s2', 'mixing'))
out_sim_s2_seq30_melt[, variable := as.factor(variable)]
out_sim_s2_seq30_melt[, s2:= as.factor(s2)]

# compare values of summary statistics for different summary stats
ggplot(out_sim_s2_seq30_melt, aes(x = s2, y = value, col = mixing))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data', col = 'Cell mixing')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cfMixing_30xSeq.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2_seq30_melt[mixing==FALSE], aes(x = s2, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_noMixing_30xSeq.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2_seq30_melt[mixing==TRUE], aes(x = s2, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data', col)
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cellMixing_30xSeq.pdf'), height = 8, width = 8)  

# results cf sequencing to depth 200
out_sim_s2_sub_seq200 = out_sim_s2[, c('s2', summary_stats_seq200, 'mixing'), with=FALSE]
out_sim_s2_seq200_melt = data.table::melt(out_sim_s2_sub_seq200, id.vars = c('s2', 'mixing'))
out_sim_s2_seq200_melt[, variable := as.factor(variable)]
out_sim_s2_seq200_melt[, s2:= as.factor(s2)]

# compare values of summary statistics for different summary stats
ggplot(out_sim_s2_seq200_melt, aes(x = s2, y = value, col = mixing))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data', col = 'Cell mixing')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cfMixing_200xSeq.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2_seq200_melt[mixing==FALSE], aes(x = s2, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_noMixing_200xSeq.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2_seq200_melt[mixing==TRUE], aes(x = s2, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data', col)
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cellMixing_200xSeq.pdf'), height = 8, width = 8)  

# save the simulation output to a file
# convert columns to character or it screams otherwise
out_sim_s2 = apply(out_sim_s2, 2, as.character)
write.table(out_sim_s2, 'Out/F5/F5_20241215_out_sim_s2_2to6_mixingVsNomixing.csv', sep = ',', quote=F, row.names=F)

# does this work for a combo s2 = 2 and p = 0.3
muts_sim=data.table(cell_x=numeric(), cell_y=numeric(), cell_gen=numeric(), cell_muts=character())
# create an empty dt to store simulation results in 
out_sim=data.table(mu1=numeric(), mu2=numeric(), sd1=numeric(), sd2=numeric(), s1=numeric(), s2=numeric(), n=numeric(), p=numeric(),
                   nr_shared=numeric(), nr_twin1=numeric(), nr_twin2 = numeric(),
                   vafs_shared_twin1=character(), vafs_shared_twin2=character(), vafs_spec_twin1=character(), vafs_spec_twin2=character())

# test out from s2 in 2 to 6 (2-6 ICM cell divisions)
for (p in c(0.3, 0.5, 0.7)){
    
    # run 100 simulations for each parameter value 
    for (rep in 1:1000){
      
      if (rep == 1){
        print(c(p, rep))
      }
      
      out = simulate_twinning(x, y, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, n = n, s1 = s1, s2 = 2, p = p)
      
      out_muts = out[[2]]
      out_sim = get_output(out_muts)
      
      out_grid = out[[1]]
      if (rep == 1000){
        draw_area(out_grid)
        ggsave(glue('FiguresAdd/F5/F5_sim_output_mu1_{mu1}_mu2_{mu2}_sd1_{sd1}_sd2_{sd2}_s1_{s1}_s2_{s2}_n_{n}_p_{p}.pdf'), width = 4, height = 4)
      }
      
    }
  }

###################################################################################################################################
# SIMULATION: run with 2 varying parameter (s2, p)

# what are the reasonable parameter values?
# x, y, sd1, sd2 = for x and y, anything big will do
# sd1, sd2 = does Henry have any intuition for this? I was thinking sd1 > sd2, but not sure otherwise
# s1 = set to 4
# s2 = 2-6?
# n = 2-6
# p = 0.25, 0.5, 0.75 (can try more fine-grained but realistically 0.2-0.8)
# mu1, mu2 = 2 for before ZGA and 0.5 for after ZGA

# set fixed parameters 
x = 100
y = 100
sd1 = 10
sd2 = 3
s1 = 3
n = 3
mu1 = 2
mu2 = 0.5

# observed statistics 
nr_shared_observed = 1
nr_spec_twin1_observed = 6
nr_spec_twin2_observed = 12
vafs_shared_twin1_observed = 0.2
vafs_shared_twin2_observed = 0.5
vafs_spec_twin1_observed = c(0.48, 0.4, 0.4, 0.4, 0.4, 0.2)
vafs_spec_twin2_observed = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
sum_vafs_spec_twin1_observed = sum(vafs_spec_twin1_observed)
sum_vafs_spec_twin2_observed = sum(vafs_spec_twin2_observed)

# set seed
seed = (runif(1, min=0, max=100) + runif(1, min=50, max=100)) * runif(1, min=0, max=100)
print(paste0("Random seed used in the simulation: ", seed)) # 2198.61151301869 run with 
set.seed(seed) 

# create an empty dt to store mutations for each cell on the phylogeny
muts_sim=data.table(cell_x=numeric(), cell_y=numeric(), cell_gen=numeric(), cell_muts=character())
# create an empty dt to store simulation results in 
out_sim=data.table(mu1=numeric(), mu2=numeric(), sd1=numeric(), sd2=numeric(), s1=numeric(), s2=numeric(), n=numeric(), p=numeric(),
                   nr_shared=numeric(), nr_twin1=numeric(), nr_twin2 = numeric(),
                   vafs_shared_twin1=character(), vafs_shared_twin2=character(), vafs_spec_twin1=character(), vafs_spec_twin2=character())

# test out from s2 in 2 to 6 (2-6 ICM cell divisions)
for (s2 in c(2:7)){
  
  for (p in c(0.3, 0.5, 0.7)){
  
    # run 100 simulations for each parameter value 
    for (rep in 1:1000){
      
      if (rep == 1){
        print(c(s2, p, rep))
      }
      
      out = simulate_twinning(x, y, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, n = n, s1 = s1, s2 = s2, p = p)
      
      out_muts = out[[2]]
      out_sim = get_output(out_muts)
      
      out_grid = out[[1]]
      if (rep == 1000){
        draw_area(out_grid)
        ggsave(glue('FiguresAdd/F5/F5_sim_output_mu1_{mu1}_mu2_{mu2}_sd1_{sd1}_sd2_{sd2}_s1_{s1}_s2_{s2}_n_{n}_p_{p}.pdf'), width = 4, height = 4)
      }
    }
  }
}

out_sim_s2p_bisect = data.table(out_sim)
out_sim_s2p_bisect[, mixing := FALSE]

# make sure that relevant columns in out_sim are numeric
num_cols = c('mu1', 'mu2', 'sd1', 'sd2', 's1', 's2', 'n', 'p', 'nr_shared', 'nr_twin1', 'nr_twin2')
out_sim_s2p_bisect[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

# ADD SEQUENCING SIMULATION TO DEPTHS 30x (1 SAMPLE) AND 200x (AGGREGATED)

# add VAFs observed "after sequencing" to 30x coverage (1 WGS sample)
out_sim_s2p_bisect[, vafs_shared_twin1_seq30 := lapply(vafs_shared_twin1, sim_sequencing, lambda = 30)] # vafs of shared mutations in twin1 
out_sim_s2p_bisect[, vafs_shared_twin2_seq30 := lapply(vafs_shared_twin2, sim_sequencing, lambda = 30)] # vafs of shared mutations in twin1
out_sim_s2p_bisect[, vafs_spec_twin1_seq30 := lapply(vafs_spec_twin1, sim_sequencing, lambda = 30)]
out_sim_s2p_bisect[, vafs_spec_twin2_seq30 := lapply(vafs_spec_twin2, sim_sequencing, lambda = 30)]
# add how many mutations would be observed based on the sequenced VAF (from sequencing to 30x)
out_sim_s2p_bisect[, nr_shared_twin1_seq30 := sapply(vafs_shared_twin1_seq30, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_s2p_bisect[, nr_shared_twin2_seq30 := sapply(vafs_shared_twin2_seq30, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_s2p_bisect[, nr_spec_twin1_seq30 := sapply(vafs_spec_twin1_seq30, function(x) length(x[x > 0.1]))]
out_sim_s2p_bisect[, nr_spec_twin2_seq30 := sapply(vafs_spec_twin2_seq30, function(x) length(x[x > 0.1]))]

# add VAFs observed "after sequencing" to 200x coverage (all WGS samples aggregated)
out_sim_s2p_bisect[, vafs_shared_twin1_seq200 := lapply(vafs_shared_twin1, sim_sequencing, lambda = 200)] # vafs of shared mutations in twin1 
out_sim_s2p_bisect[, vafs_shared_twin2_seq200 := lapply(vafs_shared_twin2, sim_sequencing, lambda = 200)] # vafs of shared mutations in twin1
out_sim_s2p_bisect[, vafs_spec_twin1_seq200 := lapply(vafs_spec_twin1, sim_sequencing, lambda = 200)]
out_sim_s2p_bisect[, vafs_spec_twin2_seq200 := lapply(vafs_spec_twin2, sim_sequencing, lambda = 200)]
# add how many mutations would be observed based on the sequenced VAF 
out_sim_s2p_bisect[, nr_shared_twin1_seq200 := sapply(vafs_shared_twin1_seq200, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_s2p_bisect[, nr_shared_twin2_seq200 := sapply(vafs_shared_twin2_seq200, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_s2p_bisect[, nr_spec_twin1_seq200 := sapply(vafs_spec_twin1_seq200, function(x) length(x[x > 0.1]))]
out_sim_s2p_bisect[, nr_spec_twin2_seq200 := sapply(vafs_spec_twin2_seq200, function(x) length(x[x > 0.1]))]

# CALCULATE DIFFERENCE TO SUMMARY STATISTICS 
# determine difference between observed and simulated data (raw - w/o "sequencing")
out_sim_s2p_bisect[, nr_shared_diff_raw := nr_shared - nr_shared_observed]
out_sim_s2p_bisect[, nr_spec_twin1_diff_raw := nr_twin1 - nr_spec_twin1_observed]
out_sim_s2p_bisect[, nr_spec_twin2_diff_raw := nr_twin2 - nr_spec_twin2_observed]
out_sim_s2p_bisect[, vafs_shared_twin1 := lapply(vafs_shared_twin1, char_to_list)]
out_sim_s2p_bisect[, vafs_shared_twin2 := lapply(vafs_shared_twin2, char_to_list)]
out_sim_s2p_bisect[, vafs_spec_twin1 := lapply(vafs_spec_twin1, char_to_list)]
out_sim_s2p_bisect[, vafs_spec_twin2 := lapply(vafs_spec_twin2, char_to_list)]
out_sim_s2p_bisect[, vaf_shared_twin1_diff_raw := mapply(get_vaf_diff, vafs_shared_twin1, vafs_shared_twin1_observed)]
out_sim_s2p_bisect[, vaf_shared_twin2_diff_raw := mapply(get_vaf_diff, vafs_shared_twin2, vafs_shared_twin2_observed)]
out_sim_s2p_bisect[, vaf_spec_twin1_diff_raw := mapply(get_vaf_diff, vafs_spec_twin1, sum_vafs_spec_twin1_observed)]
out_sim_s2p_bisect[, vaf_spec_twin2_diff_raw := mapply(get_vaf_diff, vafs_spec_twin2, sum_vafs_spec_twin2_observed)]

# difference between simulated (sequenced) and observed statistics (sequencing to 30x)
out_sim_s2p_bisect[, nr_shared_twin1_diff_seq30 := nr_shared_twin1_seq30 - nr_shared_observed]
out_sim_s2p_bisect[, nr_shared_twin2_diff_seq30 := nr_shared_twin2_seq30 - nr_shared_observed]
out_sim_s2p_bisect[, nr_spec_twin1_diff_seq30 := nr_spec_twin1_seq30 - nr_spec_twin1_observed]
out_sim_s2p_bisect[, nr_spec_twin2_diff_seq30 := nr_spec_twin2_seq30 - nr_spec_twin2_observed]
out_sim_s2p_bisect[, vaf_shared_twin1_diff_seq30 := mapply(get_vaf_diff, vafs_shared_twin1_seq30, vafs_shared_twin1_observed)]
out_sim_s2p_bisect[, vaf_shared_twin2_diff_seq30 := mapply(get_vaf_diff, vafs_shared_twin2_seq30, vafs_shared_twin2_observed)]
out_sim_s2p_bisect[, vaf_spec_twin1_diff_seq30 := mapply(get_vaf_diff, vafs_spec_twin1_seq30, sum_vafs_spec_twin1_observed)]
out_sim_s2p_bisect[, vaf_spec_twin2_diff_seq30 := mapply(get_vaf_diff, vafs_spec_twin2_seq30, sum_vafs_spec_twin2_observed)]

# difference between simulated (sequenced) and observed statistics (sequencing to 200x)
out_sim_s2p_bisect[, nr_shared_twin1_diff_seq200 := nr_shared_twin1_seq200 - nr_shared_observed]
out_sim_s2p_bisect[, nr_shared_twin2_diff_seq200 := nr_shared_twin2_seq200 - nr_shared_observed]
out_sim_s2p_bisect[, nr_spec_twin1_diff_seq200 := nr_spec_twin1_seq200 - nr_spec_twin1_observed]
out_sim_s2p_bisect[, nr_spec_twin2_diff_seq200 := nr_spec_twin2_seq200 - nr_spec_twin2_observed]
out_sim_s2p_bisect[, vaf_shared_twin1_diff_seq200 := mapply(get_vaf_diff, vafs_shared_twin1_seq200, vafs_shared_twin1_observed)]
out_sim_s2p_bisect[, vaf_shared_twin2_diff_seq200 := mapply(get_vaf_diff, vafs_shared_twin2_seq200, vafs_shared_twin2_observed)]
out_sim_s2p_bisect[, vaf_spec_twin1_diff_seq200 := mapply(get_vaf_diff, vafs_spec_twin1_seq200, sum_vafs_spec_twin1_observed)]
out_sim_s2p_bisect[, vaf_spec_twin2_diff_seq200 := mapply(get_vaf_diff, vafs_spec_twin2_seq200, sum_vafs_spec_twin2_observed)]

# repeat, but now allow for cell mixing 

# set seed
seed = (runif(1, min=0, max=100) + runif(1, min=50, max=100)) * runif(1, min=0, max=100)
print(paste0("Random seed used in the simulation: ", seed))
set.seed(seed)

# create an empty dt to store mutations for each cell on the phylogeny
muts_sim=data.table(cell_x=numeric(), cell_y=numeric(), cell_gen=numeric(), cell_muts=character())
# create an empty dt to store simulation results in 
out_sim=data.table(mu1=numeric(), mu2=numeric(), sd1=numeric(), sd2=numeric(), s1=numeric(), s2=numeric(), n=numeric(), p=numeric(),
                   nr_shared=numeric(), nr_twin1=numeric(), nr_twin2 = numeric(),
                   vafs_shared_twin1=character(), vafs_shared_twin2=character(), vafs_spec_twin1=character(), vafs_spec_twin2=character())

# test out from s2 in 2 to 6 (2-6 ICM cell divisions)
for (s2 in 2:7){
  
  for (p in c(0.3, 0.5, 0.7)){
    
    # run 100 simulations for each parameter value 
    for (rep in 1:1000){
      
      if (rep == 1000){
        print(rep)
      }
      out = simulate_twinning(x, y, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, n = n, s1 = s1, s2 = s2, p = p, cell_mixing = TRUE)
      
      out_muts = out[[2]]
      out_sim = get_output(out_muts)
      
      out_grid = out[[1]]
      if (rep == 1000){
        draw_area(out_grid)
        ggsave(glue('FiguresAdd/F5/F5_sim_output_mu1_{mu1}_mu2_{mu2}_sd1_{sd1}_sd2_{sd2}_s1_{s1}_s2_{s2}_n_{n}_p_{p}_mixing.pdf'), width = 4, height = 4)
      }
    }
  }
}

# make sure that relevant columns in out_sim are numeric
num_cols = c('mu1', 'mu2', 'sd1', 'sd2', 's1', 's2', 'n', 'p', 'nr_shared', 'nr_twin1', 'nr_twin2')
out_sim[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

out_sim_s2p_mixing = data.table(out_sim)
out_sim_s2p_mixing[, mixing := TRUE]

# add VAFs observed "after sequencing" to 30x coverage (1 WGS sample)
out_sim_s2p_mixing[, vafs_shared_twin1_seq30 := lapply(vafs_shared_twin1, sim_sequencing, lambda = 30)] # vafs of shared mutations in twin1 
out_sim_s2p_mixing[, vafs_shared_twin2_seq30 := lapply(vafs_shared_twin2, sim_sequencing, lambda = 30)] # vafs of shared mutations in twin1
out_sim_s2p_mixing[, vafs_spec_twin1_seq30 := lapply(vafs_spec_twin1, sim_sequencing, lambda = 30)]
out_sim_s2p_mixing[, vafs_spec_twin2_seq30 := lapply(vafs_spec_twin2, sim_sequencing, lambda = 30)]
# add how many mutations would be observed based on the sequenced VAF (from sequencing to 30x)
out_sim_s2p_mixing[, nr_shared_twin1_seq30 := sapply(vafs_shared_twin1_seq30, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_s2p_mixing[, nr_shared_twin2_seq30 := sapply(vafs_shared_twin2_seq30, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_s2p_mixing[, nr_spec_twin1_seq30 := sapply(vafs_spec_twin1_seq30, function(x) length(x[x > 0.1]))]
out_sim_s2p_mixing[, nr_spec_twin2_seq30 := sapply(vafs_spec_twin2_seq30, function(x) length(x[x > 0.1]))]

# add VAFs observed "after sequencing" to 200x coverage (all WGS samples aggregated)
out_sim_s2p_mixing[, vafs_shared_twin1_seq200 := lapply(vafs_shared_twin1, sim_sequencing, lambda = 200)] # vafs of shared mutations in twin1 
out_sim_s2p_mixing[, vafs_shared_twin2_seq200 := lapply(vafs_shared_twin2, sim_sequencing, lambda = 200)] # vafs of shared mutations in twin1
out_sim_s2p_mixing[, vafs_spec_twin1_seq200 := lapply(vafs_spec_twin1, sim_sequencing, lambda = 200)]
out_sim_s2p_mixing[, vafs_spec_twin2_seq200 := lapply(vafs_spec_twin2, sim_sequencing, lambda = 200)]
# add how many mutations would be observed based on the sequenced VAF 
out_sim_s2p_mixing[, nr_shared_twin1_seq200 := sapply(vafs_shared_twin1_seq200, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_s2p_mixing[, nr_shared_twin2_seq200 := sapply(vafs_shared_twin2_seq200, function(x) length(x[x > 0.1]))] # this accounts for shared mutations being mis-classified as twin specific due to sampling 
out_sim_s2p_mixing[, nr_spec_twin1_seq200 := sapply(vafs_spec_twin1_seq200, function(x) length(x[x > 0.1]))]
out_sim_s2p_mixing[, nr_spec_twin2_seq200 := sapply(vafs_spec_twin2_seq200, function(x) length(x[x > 0.1]))]

# CALCULATE DIFFERENCE TO SUMMARY STATISTICS 
# determine difference between observed and simulated data (raw - w/o "sequencing")
out_sim_s2p_mixing[, nr_shared_diff_raw := nr_shared - nr_shared_observed]
out_sim_s2p_mixing[, nr_spec_twin1_diff_raw := nr_twin1 - nr_spec_twin1_observed]
out_sim_s2p_mixing[, nr_spec_twin2_diff_raw := nr_twin2 - nr_spec_twin2_observed]
out_sim_s2p_mixing[, vafs_shared_twin1 := lapply(vafs_shared_twin1, char_to_list)]
out_sim_s2p_mixing[, vafs_shared_twin2 := lapply(vafs_shared_twin2, char_to_list)]
out_sim_s2p_mixing[, vafs_spec_twin1 := lapply(vafs_spec_twin1, char_to_list)]
out_sim_s2p_mixing[, vafs_spec_twin2 := lapply(vafs_spec_twin2, char_to_list)]
out_sim_s2p_mixing[, vaf_shared_twin1_diff_raw := mapply(get_vaf_diff, vafs_shared_twin1, vafs_shared_twin1_observed)]
out_sim_s2p_mixing[, vaf_shared_twin2_diff_raw := mapply(get_vaf_diff, vafs_shared_twin2, vafs_shared_twin2_observed)]
out_sim_s2p_mixing[, vaf_spec_twin1_diff_raw := mapply(get_vaf_diff, vafs_spec_twin1, sum_vafs_spec_twin1_observed)]
out_sim_s2p_mixing[, vaf_spec_twin2_diff_raw := mapply(get_vaf_diff, vafs_spec_twin2, sum_vafs_spec_twin2_observed)]

# difference between simulated (sequenced) and observed statistics (sequencing to 30x)
out_sim_s2p_mixing[, nr_shared_twin1_diff_seq30 := nr_shared_twin1_seq30 - nr_shared_observed]
out_sim_s2p_mixing[, nr_shared_twin2_diff_seq30 := nr_shared_twin2_seq30 - nr_shared_observed]
out_sim_s2p_mixing[, nr_spec_twin1_diff_seq30 := nr_spec_twin1_seq30 - nr_spec_twin1_observed]
out_sim_s2p_mixing[, nr_spec_twin2_diff_seq30 := nr_spec_twin2_seq30 - nr_spec_twin2_observed]
out_sim_s2p_mixing[, vaf_shared_twin1_diff_seq30 := mapply(get_vaf_diff, vafs_shared_twin1_seq30, vafs_shared_twin1_observed)]
out_sim_s2p_mixing[, vaf_shared_twin2_diff_seq30 := mapply(get_vaf_diff, vafs_shared_twin2_seq30, vafs_shared_twin2_observed)]
out_sim_s2p_mixing[, vaf_spec_twin1_diff_seq30 := mapply(get_vaf_diff, vafs_spec_twin1_seq30, sum_vafs_spec_twin1_observed)]
out_sim_s2p_mixing[, vaf_spec_twin2_diff_seq30 := mapply(get_vaf_diff, vafs_spec_twin2_seq30, sum_vafs_spec_twin2_observed)]

# difference between simulated (sequenced) and observed statistics (sequencing to 200x)
out_sim_s2p_mixing[, nr_shared_twin1_diff_seq200 := nr_shared_twin1_seq200 - nr_shared_observed]
out_sim_s2p_mixing[, nr_shared_twin2_diff_seq200 := nr_shared_twin2_seq200 - nr_shared_observed]
out_sim_s2p_mixing[, nr_spec_twin1_diff_seq200 := nr_spec_twin1_seq200 - nr_spec_twin1_observed]
out_sim_s2p_mixing[, nr_spec_twin2_diff_seq200 := nr_spec_twin2_seq200 - nr_spec_twin2_observed]
out_sim_s2p_mixing[, vaf_shared_twin1_diff_seq200 := mapply(get_vaf_diff, vafs_shared_twin1_seq200, vafs_shared_twin1_observed)]
out_sim_s2p_mixing[, vaf_shared_twin2_diff_seq200 := mapply(get_vaf_diff, vafs_shared_twin2_seq200, vafs_shared_twin2_observed)]
out_sim_s2p_mixing[, vaf_spec_twin1_diff_seq200 := mapply(get_vaf_diff, vafs_spec_twin1_seq200, sum_vafs_spec_twin1_observed)]
out_sim_s2p_mixing[, vaf_spec_twin2_diff_seq200 := mapply(get_vaf_diff, vafs_spec_twin2_seq200, sum_vafs_spec_twin2_observed)]

# SAVE SIMULATION RESULTS
# save simulation results to out
out_sim_s2p = rbind(out_sim_s2p_bisect, out_sim_s2p_mixing)

# define which summary statistics are of interest
summary_stats_raw = c("nr_shared_diff_raw", "nr_spec_twin1_diff_raw",   
                      "nr_spec_twin2_diff_raw",  "vaf_shared_twin1_diff_raw", "vaf_shared_twin2_diff_raw",
                      "vaf_spec_twin1_diff_raw",   "vaf_spec_twin2_diff_raw")
summary_stats_seq30 = c("nr_shared_twin1_diff_seq30",  "nr_shared_twin2_diff_seq30",  "nr_spec_twin1_diff_seq30",   
                        "nr_spec_twin2_diff_seq30",  "vaf_shared_twin1_diff_seq30", "vaf_shared_twin2_diff_seq30",
                        "vaf_spec_twin1_diff_seq30",   "vaf_spec_twin2_diff_seq30")
summary_stats_seq200 = c("nr_shared_twin1_diff_seq200",  "nr_shared_twin2_diff_seq200",  "nr_spec_twin1_diff_seq200",   
                         "nr_spec_twin2_diff_seq200",  "vaf_shared_twin1_diff_seq200", "vaf_shared_twin2_diff_seq200",
                         "vaf_spec_twin1_diff_seq200",   "vaf_spec_twin2_diff_seq200")

# melt results dt
out_sim_s2p_sub_raw = out_sim_s2p[, c('s2', 'p', summary_stats_raw, 'mixing'), with=FALSE]
out_sim_s2p_raw_melt = data.table::melt(out_sim_s2p_sub_raw, id.vars = c('s2', 'p', 'mixing'))
out_sim_s2p_raw_melt[, variable := as.factor(variable)]
out_sim_s2p_raw_melt[, s2:= as.factor(s2)]
out_sim_s2p_raw_melt[, p:= as.factor(p)]

# compare values of summary statistics for different summary stats
ggplot(out_sim_s2p_raw_melt, aes(x = s2, y = value, col = mixing))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 3)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data', col = 'Cell mixing')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cfMixing_uncorrectedSeq_s2.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_raw_melt[mixing==FALSE], aes(x = s2, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 3)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_noMixing_uncorrectedSeq_s2.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_raw_melt[mixing==TRUE], aes(x = s2, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 3)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data', col)
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cellMixing_uncorrectedSeq_s2.pdf'), height = 8, width = 8)  

# plots to test the importance of asymmetry 
ggplot(out_sim_s2p_raw_melt, aes(x = p, y = value, col = mixing))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 3)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Proportion of cells allocated to twin1', y = 'Difference to observed data', col = 'Cell mixing')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cfMixing_uncorrectedSeq_p.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_raw_melt[mixing==FALSE], aes(x = p, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 3)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Proportion of cells allocated to twin1', y = 'Difference to observed data')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_noMixing_uncorrectedSeq_p.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_raw_melt[mixing==TRUE], aes(x = p, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 3)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Proportion of cells allocated to twin1', y = 'Difference to observed data', col)
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cellMixing_uncorrectedSeq_p.pdf'), height = 8, width = 8)  

# results (sequencing to 30x coverage)
out_sim_s2p_sub_seq30 = out_sim_s2p[, c('s2', 'p', summary_stats_seq30, 'mixing'), with=FALSE]
out_sim_s2p_seq30_melt = data.table::melt(out_sim_s2p_sub_seq30, id.vars = c('s2', 'p', 'mixing'))
out_sim_s2p_seq30_melt[, variable := as.factor(variable)]
out_sim_s2p_seq30_melt[, s2:= as.factor(s2)]
out_sim_s2p_seq30_melt[, p:= as.factor(p)]

# compare values of summary statistics for different summary stats
ggplot(out_sim_s2p_seq30_melt, aes(x = s2, y = value, col = mixing))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data', col = 'Cell mixing')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cfMixing_30xSeq_s2.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_seq30_melt[mixing==FALSE], aes(x = s2, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_noMixing_30xSeq_s2.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_seq30_melt[mixing==TRUE], aes(x = s2, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data', col)
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cellMixing_30xSeq_s2.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_seq30_melt, aes(x = p, y = value, col = mixing))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Proportion of cells allocated to twin1', y = 'Difference to observed data', col = 'Cell mixing')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cfMixing_30xSeq_p.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_seq30_melt[mixing==FALSE], aes(x = p, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Proportion of cells allocated to twin1', y = 'Difference to observed data')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_noMixing_30xSeq_p.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_seq30_melt[mixing==TRUE], aes(x = p, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Proportion of cells allocated to twin1', y = 'Difference to observed data', col)
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cellMixing_30xSeq_p.pdf'), height = 8, width = 8)  

# results cf sequencing to depth 200
out_sim_s2p_sub_seq200 = out_sim_s2p[, c('s2', 'p', summary_stats_seq200, 'mixing'), with=FALSE]
out_sim_s2p_seq200_melt = data.table::melt(out_sim_s2p_sub_seq200, id.vars = c('s2', 'p', 'mixing'))
out_sim_s2p_seq200_melt[, variable := as.factor(variable)]
out_sim_s2p_seq200_melt[, s2:= as.factor(s2)]
out_sim_s2p_seq200_melt[, p:= as.factor(s2)]

# compare values of summary statistics for different summary stats
ggplot(out_sim_s2p_seq200_melt, aes(x = s2, y = value, col = mixing))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data', col = 'Cell mixing')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cfMixing_200xSeq_s2.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_seq200_melt[mixing==FALSE], aes(x = s2, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_noMixing_200xSeq_s2.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_seq200_melt[mixing==TRUE], aes(x = s2, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Number of post-ICM divisions', y = 'Difference to observed data', col)
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cellMixing_200xSeq_s2.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_seq200_melt, aes(x = p, y = value, col = mixing))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Proportion of cells allocated to twin1', y = 'Difference to observed data', col = 'Cell mixing')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cfMixing_200xSeq_p.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_seq200_melt[mixing==FALSE], aes(x = p, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Proportion of cells allocated to twin1', y = 'Difference to observed data')
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_noMixing_200xSeq_p.pdf'), height = 8, width = 8)  

ggplot(out_sim_s2p_seq200_melt[mixing==TRUE], aes(x = p, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, nrow = 4)+
  theme_classic(base_size = 13)+
  geom_hline(yintercept = 0, color = 'black')+
  labs(x = 'Proportion of cells allocated to twin1', y = 'Difference to observed data', col)
ggsave(glue('FiguresAdd/F5/F5_diff_to_stat_cellMixing_200xSeq_p.pdf'), height = 8, width = 8)  

# save the simulation output to a file
# convert columns to character or it screams otherwise
out_sim_s2p = rbind(out_sim_s2p_bisect, out_sim_s2p_mixing)
out_sim_s2p_char = apply(out_sim_s2p, 2, as.character)
write.table(out_sim_s2p_char, 'Out/F5/F5_20241215_out_sim_s2_2to7_p_3values_mixingVsNomixing.txt', sep = '\t', quote=F, row.names=F)

###################################################################################################################################
# Heatmap 

# create a matrix s2 vs p with difference 
summary_stats_diffs = c("nr_shared_diff_raw", "nr_spec_twin1_diff_raw",   
                        "nr_spec_twin2_diff_raw",  "vaf_shared_twin1_diff_raw", "vaf_shared_twin2_diff_raw",
                        "vaf_spec_twin1_diff_raw",   "vaf_spec_twin2_diff_raw", "nr_shared_twin1_diff_seq30", 
                        "nr_shared_twin2_diff_seq30",  "nr_spec_twin1_diff_seq30","nr_spec_twin2_diff_seq30",  
                        "vaf_shared_twin1_diff_seq30", "vaf_shared_twin2_diff_seq30", "vaf_spec_twin1_diff_seq30", 
                        "vaf_spec_twin2_diff_seq30", "nr_shared_twin1_diff_seq200",  "nr_shared_twin2_diff_seq200",  
                        "nr_spec_twin1_diff_seq200",  "nr_spec_twin2_diff_seq200",  "vaf_shared_twin1_diff_seq200", 
                        "vaf_shared_twin2_diff_seq200", "vaf_spec_twin1_diff_seq200",   "vaf_spec_twin2_diff_seq200")

# heatmap for each summary statistic difference to observed 
for (stat_diff in summary_stats_diffs){
  
  # subset the dt to only include s2, p and the specific statistic
  mat = out_sim_s2p[, c('s2', 'p', stat_diff), with=FALSE]
  mat_abs = abs(mat)
  
  mat = dcast(mat, s2 ~ p, fun = mean) # mean difference from 1000 simulations
  mat_abs = dcast(mat_abs, s2 ~ p, fun = mean) # mean difference from 1000 simulations
  
  mat = mat %>% remove_rownames %>% column_to_rownames(var="s2")
  mat_abs = mat_abs %>% remove_rownames %>% column_to_rownames(var="s2")
  
  # Plot output of the simulation as a heatmap (with distance to observed in color)
  pdf(glue("FiguresAdd/F5/F5_heatmap_s2_vs_p_distToObserved_{stat_diff}.pdf"), width=4, height=4)
  pheatmap(mat, 
           main=glue("{stat_diff}"), 
           treeheight_row = 0,
           cluster_rows = F, cluster_cols = F, 
           show_rownames = T, show_colnames = T)
  dev.off()
  
  pdf(glue("FiguresAdd/F5/F5_heatmap_s2_vs_p_distToObserved_{stat_diff}_absValues.pdf"), width=4, height=4)
  pheatmap(mat, 
           main=glue("{stat_diff}"), 
           treeheight_row = 0,
           cluster_rows = F, cluster_cols = F, 
           show_rownames = T, show_colnames = T)
  dev.off()
  
}




