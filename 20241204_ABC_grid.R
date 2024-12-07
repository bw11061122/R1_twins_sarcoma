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
# Set working directory 
setwd('/Users/bw18/Desktop/1SB/ABC')

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
# SIMULATION: ADDING MUTATIONS (IDs = RANDOM STRINGS)

# create 10k strings that can be used as mutation IDs to assign to cells
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
                        values = c("black", "#c81048", "#10c8c5", "purple", "#dba20a")) # black if absent, red if present 
  return(p)
}


# Helper function to check indices selected
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

# Divide cells pre-ICM
# divide all cells on the grid (for each existing cell, create two daughter cells, then remove the existing one)
# use this function when dividing cells to contribute to the ICM 

# for each cell existing on the grid, this function selects coordinates of two new cells
# check that coordinates are not outside of grid bounds and do not overlap 

# each daughter cell acquired a number of mutations selected from the Poisson distribution according to the mutation rate

simulate_twinning = function(x, y, mu1 = 1, mu2 = 1, sd1 = 10, sd2 = 3, s1 = 4, n = 3, s2 = 4, p = 0.5, center = TRUE){
  
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
  
  # add mutation to the first cell (give it 1 mutation so we can keep track of who's whose parent)
  
  # select mutations from the list
  muts_acquired = sample(muts_ids, 1)
  
  # update muts_ids to avoid sampling the same mutation twice
  muts_ids = setdiff(muts_ids, muts_acquired)
  
  # add mutations to the dt 
  muts_row = t(c(xNew, yNew, 1, muts_acquired)) # first row = 1st cell in the embryo 
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
            nr_muts = rpois(lambda=mu1, n=1) # only return one number 
            
            # select mutations from the list
            muts_acquired1 = sample(muts_ids, nr_muts)
            muts_all1 = c(muts_parent, muts_acquired1)
            muts_all1 = paste(muts_all1, collapse = ',')
            
            # update muts_ids to avoid sampling the same mutation twice
            muts_ids = setdiff(muts_ids, muts_all1)
            
            # add mutations to the dt 
            muts_row = t(c(xNew1, yNew1, iter+1, muts_all1))
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
            nr_muts = rpois(lambda=mu1, n=1) # only return one number 
            
            # select mutations from the list
            muts_acquired2 = sample(muts_ids, nr_muts)
            muts_all2 = c(muts_parent, muts_acquired2)
            muts_all2 = paste(muts_all2, collapse = ',')
            
            # update muts_ids to avoid sampling the same mutation twice
            muts_ids = setdiff(muts_ids, muts_all2)
            
            # add mutations to the dt 
            muts_row = t(c(xNew2, yNew2, iter+1, muts_all2))
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
  last_gen = s1+1
  coords_cells = muts_sim[cell_gen==last_gen, c('cell_x', 'cell_y'), with=FALSE]
  coords_cells[, names(coords_cells) := lapply(.SD, as.numeric), .SDcols = names(coords_cells)]
  
  # it would be good to add a check to prevent cells generated from the same parent cell to be both ICM allocated
  # based on the idea that cells ingressed through asymmetric divisions contribute to the ICM
  # > if the division is asymmetric, only one of the two cells will be ICM-allocated 
  
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
  sum_cells_before = dim(muts_sim[cell_gen<=s1])[1]
  accepted_indices_2 = accepted_indices + sum_cells_before
  muts_sim[, makes_ICM := ifelse(.I %in% accepted_indices_2, 1, 0)]
  
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
 
  nr_cells_twins = sum(grid==2)
  nr_cells_twin1 = round(p * nr_cells_twins)
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
    !coord %in% c(ctwin1, ctwin2), 'not ICM'
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

###################################################################################################################################
# SIMULATION: output 

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

###################################################################################################################################
# SIMULATION: parameters 

# what are the reasonable parameter values?
# x, y, sd1, sd2 = for x and y, anything big will do
# sd1, sd2 = does Henry have any intuition for this? I was thinking sd1 > sd2, but not sure otherwise
# s1 = 3-5?
# s2 = 2-6?
# n = 2-6
# p = 0.25, 0.5, 0.75 (can try more fine-grained but realistically 0.2-0.8)
# mu1, mu2 = for now set both to 1, can think about being more fancy later

# specify parameters and test for s2 = 2, 3, 4, 5, 6
x = 100
y = 100
mu1 = 1
mu2 = 1
sd1 = 10
sd2 = 3
s1 = 3
n = 3
p = 0.5

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
  for (rep in 1:100){
    out = simulate_twinning(x, y, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, n = n, s1 = s1, s2 = s2, p = p)
    
    out_muts = out[[2]]
    out_sim = get_output(out_muts)
    
    out_grid = out[[1]]
    if (rep == 100){
      draw_area(out_grid)
      ggsave(glue('Results/20241206_sim_output_mu1_{mu1}_mu2_{mu2}_sd1_{sd1}_sd2_{sd2}_s1_{s1}_s2_{s2}_n_{n}_p_{p}.pdf'), width = 4, height = 4)
    }
    
  }
}

# make sure that relevant columns in out_sim are numeric
num_cols = c('mu1', 'mu2', 'sd1', 'sd2', 's1', 's2', 'n', 'p', 'nr_shared', 'nr_twin1', 'nr_twin2')
out_sim[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

# how many mutations (shared / twin-specific) do we get for each s?
for (s in 2:6){
  
  out_sim_sub = out_sim[s2 == s]
  nr_shared_muts = out_sim_sub[, nr_shared] %>% unlist() %>% as.numeric()
  nr_twin1_muts = out_sim_sub[, nr_twin1] %>% unlist() %>% as.numeric()
  nr_twin2_muts = out_sim_sub[, nr_twin2] %>% unlist() %>% as.numeric()
  
  print(hist(nr_shared_muts, xlab = 'Number of shared mutations', main = glue('s2 = {s}')))
  print(hist(nr_twin1_muts, xlab = 'Number of twin1-specific mutations', main = glue('s2 = {s}')))
  print(hist(nr_twin2_muts, xlab = 'Number of twin2-specific mutations', main = glue('s2 = {s}')))
  
}

# empty dt to write results to
muts_sim=data.table(cell_x=numeric(), cell_y=numeric(), cell_gen=numeric(), cell_muts=character())
out_sim=data.table(mu1=numeric(), mu2=numeric(), sd1=numeric(), sd2=numeric(), s1=numeric(), s2=numeric(), n=numeric(), p=numeric(),
                   nr_shared=numeric(), nr_twin1=numeric(), nr_twin2 = numeric(),
                   vafs_shared_twin1=character(), vafs_shared_twin2=character(), vafs_spec_twin1=character(), vafs_spec_twin2=character())


# test s2 and p (0.25, 0.5, 0.75) (parameter screen kind of)
for (s2 in 2:6){
  
  for (pf in c(0.3, 0.5, 0.7)){
  
    for (rep in 1:100){
      
      # run 100 simulations for each parameter value 
      out = simulate_twinning(x, y, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, n = n, s1 = s1, s2 = s2, p = pf)
      
      out_muts = out[[2]]
      out_sim = get_output(out_muts)
      
      if (rep==1){
        print(c(s2, pf))
        out_grid = out[[1]]
        draw_area(out_grid)
        ggsave('Results/20241207_out_sim')
        }
      
      
      
    }
  }
}
  
num_cols = c('mu1', 'mu2', 'sd1', 'sd2', 's1', 's2', 'n', 'p', 'nr_shared', 'nr_twin1', 'nr_twin2')
out_sim[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
  
# how many mutations (shared / twin-specific) do we get for each s?
for (s_value in 2:6){
  for (p_value in c(0.3, 0.5, 0.7)){
    
    out_sim_sub = out_sim[s2 == s_value & p == p_value]
    nr_shared_muts = out_sim_sub[, nr_shared] %>% unlist() %>% as.numeric()
    nr_twin1_muts = out_sim_sub[, nr_twin1] %>% unlist() %>% as.numeric()
    nr_twin2_muts = out_sim_sub[, nr_twin2] %>% unlist() %>% as.numeric()
    
    print(hist(nr_shared_muts, xlab = 'Number of shared mutations', main = glue('s2 = {s_value}, p = {p_value}'), xlim = c(0, 20), ylim = c(0, 80), breaks = 10))
    print(hist(nr_twin1_muts, xlab = 'Number of twin1-specific mutations', main = glue('s2 = {s_value}, p = {p_value}'), xlim = c(0, 20), ylim = c(0, 80), breaks = 10))
    print(hist(nr_twin2_muts, xlab = 'Number of twin2-specific mutations', main = glue('s2 = {s_value}, p = {p_value}'), xlim = c(0, 20), ylim = c(0, 80), breaks = 10))
    
  }  
}

# ngl this is really slow but this is the best I think I can do coding-wise

###################################################################################################################################
# SIMULATION: target data 

###################################################################################################################################
# SIMULATION: ABC figuring this out 
observed_data=data.frame(nr_shared = 1, 
                         nr_twin1 = 7,
                         nr_twin2 = 11)  

abc_results=abc(target=observed_data, 
                param=all_sim[select,c("s1","s2","n","p")], 
                sumstat=all_sim[select,c("nr_shared","nr_twin1","nr_twin2")], 
                tol=0.01, hcorr=F,method="rejection")





