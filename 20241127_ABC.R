####################################################################################################################################
# SCRIPT 7

# Getting started with ABC simulations
# 2024-11-27
# Barbara Walkowiak bw18

# Script to run simulations (based on embryonic_bottleneck_ABC.R from Tim Coorens' Extensive phylogenies paper + Henry Lee-Six HSC)

###################################################################################################################################
# LIBRARIES

# Load needed libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(glue)
library(stringr)
library(pheatmap)
library(ggrepel)
library(abc)
library(ape)
library(ggtree)
library(phytools)

###################################################################################################################################
# INPUT FILES 

setwd('/Users/bw18/Desktop/1SB/ABC')

options(stringsAsFactors = F)
iter_num=as.numeric(commandArgs(T)[1])

# what you need to do 

# okay how I imagine this simulation
# I simulate a tree of the ICM
# simulate the embryo (10 cell divisions, so 2^10 cells)
# select n cells to the ICM at stage s 
# but after that, you simulate the ICM cells for a bit more (say to get max 100 cells)
# starting from the baseline ICM cells to each cell division, split cells into two groups (w/ asymmetry defined by p)
# if you take 4 cells to the ICM and split w/ 1/4, this means 1 cells out of 4 makes twin 1 and the rest makes twin 2
# if you take say 4 cells to the ICM, let them divide 4 times (so you have 64 cells) and split w/ 1/2, this means each twin starts with 32 cells 
# you basically want to conceptualize the different options and what each one returns 
# I am not sure how to do this in a way that "remembers" which cell has which mutation and smhw splits the tree spatially (cause there is the funky cell mixing parameter I coudl also test, but would be easier to say there is no cell mixing at all for the start)
# Henry suggested simulating a grid, which is really smart but not sure how to start implementing this

# parameters
# 1 mutation rate before ZGA (mu1)
# 2 mutation rate after ZGA (mu2)
# 3 number of cells taken to the ICM 
# 4 stage at which cells selected to the ICM (nr of cell divisions until cells selected to the ICM)
# 5 number of cell divisions of the ICM cells to form the embryo (start with parameter 3 value and divide your cells x times, assume each cell divides the same nr of times)
# 6 split asymmetry = fraction of cells (obtained in stage 5) which are allocated to twin 1; the rest goes to twin 2
# 7 cell death rate = 0 for now

# assumptions
# only the ICM splits to form twins (ie twins separate after the ICM forms)
# cells that split stay as ICM, don't become TE cells or sth 

# summary statistics 
# say number of mutations specific to twin 1; number of mutations specific to twin 2; number of shared mutations
# VAF of twin1-specific, VAF of twin2-specific, VAF of shared mutations 

# how do I simulate the acquisition of mutations on a tree? slightly confused by this atm

# define needed functions
odd = function(x) x%%2 != 0 # return TRUE if odd
even = function(x) x%%2 == 0 # return TRUE if even 

# create a tree representation of a given stage (= nr cell divisions) in the embryo 
make_tree=function(stage){
  vec=1:(2^stage) # this is the nr of cells you need 
  for(n in 1:stage){
    vec=paste0("(",vec[odd(1:length(vec))],
               ",",vec[even(1:length(vec))],")")
  }
  tree=read.tree(text=paste0(vec,";"))
  return(tree)
}

# find the numbers of the parent nodes given the nr of cell divisions and the tree 
find_nodes_gen=function(stage, tree){
  parent_nodes=Ntip(tree)+1 # number of tips + 1
  for(s in 1:stage){
    parent_nodes=tree$edge[(tree$edge[,1]%in%parent_nodes),2]
  }
  return(parent_nodes)
}

tree=make_tree(stage=3) # generate a tree of the embryo at stage 3 (8 cells)
find_nodes_gen(stage=3, tree) # find numbers of nodes at stage 3 (the last 8 cells)
find_nodes_gen(stage=2, tree) # find numbers of nodes of stage 2 (4 cells)
find_nodes_gen(stage=1, tree) # find numbers of nodes of stage 1 (2 cells)

# find children of nodes that you specify as input 
find_children = function(nodes, tree=tree_df){
  all_tips=c()
  for(node in nodes){
    child_nodes = tree$node[tree$parent==node]
    tips = c()
    for (a in 1:length(child_nodes)){
      if (tree$isTip[tree$node==child_nodes[a]]){
        tips=c(tips,child_nodes[a])
      }
      else{
        tips=c(tips,find_children(child_nodes[a],tree))
      }
    }
    all_tips=c(all_tips,tips)
  }
  return(all_tips)
}

# how many cells to allocate to the ICM? 
prob_selection=c() # vector of probabilities that a given nr of cell is allocated to the LCM

# what is the difference between number and picked here??
for(stage in 2:6){ # we allow cells to be selected to the ICM from stage 2 (4 cells) to stage 6 (64 cells)
  for(number in 2:(2^(stage)-1)){ # we have from 2 to n-1 cells available (n = total nr of cells at that stage)
    probs=c()
    for(k in 1:min(number,2^(stage-1))){ # do we select cells to the ICM here? (ma 1/2 cells can be allocated to the ICM)
      probs=c(probs,binom.test(x=342,n=500,p = k/number)$p.value) 
      # given the nr of trials and the nr of successes, how unreasonable is the probability?
      # where did we take the nr of successes and the nr of trials from?
    }
    prob_selection=rbind(prob_selection,c(stage,number,max(probs),which.max(probs)))
  }
}

colnames(prob_selection)=c("stage","number","pval","picked") # picked refers to the nr of cells picked to the ICM
prob_selection=as.data.frame(prob_selection)
prob_selection_final=prob_selection[prob_selection$pval>0.01,]

# simulate a tree where you select cells to the ICM and then divide them for several more times 
sim_tree=function(stage,stage_twinning,number,mu1,mu2, split_asymmetry){
  
  tree=bi_tree(11) # create a tree (until 10 cell divisions) 
  # the logic here is that we set the boundary of 100 ICM cells to split into twins 
  # the min nr of cells allocated to the ICM is 3 (selected at stage 4); the max nr of cells allocated to the ICM is 30 (at stage 6)
  # in either case, at stage 11, we should have at least 96 cells (if 3 selected at stage 6) so we can do the split 
  
  tree$edge.length = rpois(lambda=mu2,n=nrow(tree$edge)) # add mutations onto each branch of the tree (sample from Poisson distribution with lambda = mutation rate after ZGA)
  tree$edge.length[tree$edge[,2]%in%c(find_nodes_gen(1, tree=tree),find_nodes_gen(2, tree=tree))]=rpois(lambda=mu1,n=6) 
  # ZGA at the 8 cell stage, so for the first 6 branches (in stage 1 and stage 2), replace with mutations generated with mu1 (basically higher, so we can have more mutations at this stage)
  tree_df=as.data.frame(fortify(tree)) # convert to a dt so this is nice to work with
  
  nodes_gen=find_nodes_gen(stage, tree=tree) # find nodes of the stage of interest 
  drop_nodes=sample(nodes_gen,size=2^stage-number) # drop nodes at this stage (all cells - nr of cells allocated to the LCM)
  drop_tips=find_children(drop_nodes,tree_df) # drop tips of the nodes that you dropped 
  
  tree_subset=drop.tip(tree,tip=drop_tips) # create the subsetted tree
  tree_df=as.data.frame(fortify(tree_subset)) # dt of the current tree
  
  # now we have this tree, we can select the nr of cells in the ICM we want to split 
  cells_to_split=find_nodes_gen(stage_twinning, tree=tree)
  drop_tips=find_children(cells_to_split,tree_df)
  # drop any cells after this stage
  tree_sub=drop.tip(tree,tip=drop_tips)
  
  # select nodes to be allocated to twin 1
  # 1 get nr of cells to go to twin 1 given the required split asymmetry 
  nr_cells_twin1 = round(split_asymmetry * length(find_nodes_gen(5, tree=tree)))
  nr_cells_twin2 = length(find_nodes_gen(5, tree=tree)) - nr_cells_twin1
  # 2 select cells to go to twin 1
  # here you can kind of have a mixing parameter = do you want to sample cells at random or preferentially from the same branch?
  # you could do sth like take sisters w prob x, and reduce the prob at each step back on the tree (assuming higher chance cells mixed)
  tips= # identify tips
  cells_twin1 = sample(tips, size=nr_cells_twin1) # these are the cells / nodes that you allocated to twin 1
  # here when you sample, you just sample at random, so this assumes no spatial structure in the embryo
  
  # is there a way to sample but only from one branch?
  # you can calculate the 

  # get the summary statistics out 
  
  # Number of mutations 
  
  # sum up mutations present on twin1 tree
  
  # sum up mutations present on twin2 tree
  
  # how many mutations are shared  
  
  # VAF of mutations 
  
  # calculate summary statistics for this simulation 
  orig_furcation=sum(tree_df$parent%in%which(tree_df$x==0)&tree_df$x>0)
  multiforcation_score=sum(tree_df$branch.length[tree_df$x<=5&tree_df$node!=Ntip(tree)+1]<1)/sum(tree_df$x<=5&tree_df$node!=Ntip(tree)+1)
  if(is.nan(multiforcation_score))multiforcation_score=0
  d1=length(find_children(nodes=find_nodes_gen(1, tree=tree_subset)[1],tree=tree_df))
  asymmetry=max(d1,Ntip(tree_subset)-d1)/Ntip(tree_subset)
  burden_in_first_branches=tree_df$branch.length[tree_subset$edge[(tree_subset$edge[,1]==Ntip(tree_subset)+1),2]]
  return(c(orig_furcation,asymmetry,max(burden_in_first_branches),min(burden_in_first_branches),multiforcation_score, mu1, mu2, stage, number))
}


# 
tree=bi_tree(10) # create a tree (until 10 cell divisions)
tree$edge.length = rpois(lambda=mu2,n=nrow(tree$edge)) # add mutations onto each branch of the tree (sample from Poisson distribution with lambda = mutation rate after ZGA)
tree$edge.length[tree$edge[,2]%in%c(find_nodes_gen(1, tree=tree),find_nodes_gen(2, tree=tree))]=rpois(lambda=mu1,n=6) 
# ZGA at the 8 cell stage, so for the first 6 branches, replace with mutations generated with mu1 (basically higher, so we can have more mutations at this stage)
tree_df=as.data.frame(fortify(tree)) # convert to a dt so this is nice to work with

nodes_gen=find_nodes_gen(4, tree=tree) # find nodes of the stage of interest 
drop_nodes=sample(nodes_gen,size=2^4-3) # drop nodes at this stage (all cells - nr of cells allocated to the LCM)
drop_tips=find_children(drop_nodes,tree_df) # drop tips of the nodes that you dropped 

tree_subset=drop.tip(tree,tip=drop_tips) # create the subsetted tree
tree_df=as.data.frame(fortify(tree_subset)) # dt of the current tree


# run the simulation and save summary statistics 
stage_number_mat=c()
for(stage in 2:6){
  num_range=2:min(30,max(2,2^(stage)-1)) # create a matrix for each stage and nr of cells (cap at 30)
  stage_number_mat=rbind(stage_number_mat,data.frame(stage=stage,number=num_range))
}

n_iter=5000
sim_results=matrix(0,ncol=9,nrow=n_iter) # save results for each of the simulations to a matrix 
for(k in 1:n_iter){ # for each iteration 
  mu1=sample(seq(0.5,6,by=0.01),1) # select mu1 from uniform distribution 
  mu2=sample(seq(0.1,1.5,by=0.01),1) # select mu2 from uniform distribution 
  select=sample(1:nrow(stage_number_mat),1) # randomly select a row of the matrix for stages and cells 
  stage=stage_number_mat$stage[select] # get the stage from this row  
  number=stage_number_mat$number[select] # get the number from this row 
  sim_results[k,]=sim_tree(mu1=mu1,mu2=mu2,stage=stage,number=number) # simulate tree for the numbers you've selected 
  
  if (k%%1000==0){ # print each 1000 simulations to see how things are going 
      print(k)}
}


colnames(sim_results)=c("orig_furcation","asymmetry","max_burden_branch","min_burden_branch","multifurcation","mu1","mu2","stage","number")
sim_results=as.data.frame(sim_results)
# write.table(sim_results,paste0("sim_results_",iter_num,"_round2.txt"))

# create a vector of observed summary statistics for your dataset 
twin_stats=data.frame(orig_furcation=2,
                         asymmetry=0.685,
                         min_burden_branch=1,
                         max_burden_branch=5,
                         multifurcation=0.8064516)  

# run ABC to find reasonable parameters
abc_results=abc(target=twin_stats, # vector of OBSERVED summary stats (your data)
                param=all_sim[1:1000,c("mu1","mu2","stage","number")], # these are your parameters with which you simulated the thing 
                sumstat=all_sim[1:1000,c("orig_furcation","asymmetry","min_burden_branch","max_burden_branch","multifurcation")], # these are the summary stats that you get if you run the simulation with those parameters 
                tol=0.01, hcorr=F,method="rejection") # hcorr = do you want to apply the conditional heterodcedastic model? > need to read up on this but the vignette does not explain why we would apply this as a default

# Plot the output 
stage_number_mat=c()
for(stage in 2:6){
  num_range=2:min(30,max(2,2^(stage)-1))
  stage_number_mat=rbind(stage_number_mat,data.frame(stage=stage,number=num_range))
}

cell_alloc_mat=matrix(ncol=6, nrow=30)

for(k in 1:nrow(stage_number_mat)){
  
  n=stage_number_mat[k,"number"]
  s=stage_number_mat[k,"stage"]
  
  cell_alloc_mat[n,s]=sum(abc_results$unadj.values[,c("number")]==n&abc_results$unadj.values[,"stage"]==s)
  
}

library(RColorBrewer)
library(gplots)
pdf("Results/twins_abc.pdf",width=3.5,height=8)
heatmap.2(cell_alloc_mat,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          dendrogram="none",
          na.color = "gray",
          col=colorRampPalette(brewer.pal(9,"Reds"))(100),
          Rowv = F, Colv = F, scale='none',
          key=F, cexRow=0.8, cexCol=1, notecex=1.0,
          mar=c(10,6.5))
dev.off()


###################################################################################################################################
# ME TRYING TO UNDERSTAND WHAT IS GOING ON

# Method section from the "extensive phylogenies" paper
# use a simple model of TE / ICM segregation and ABC to estimate likely parameter values

# Parameters of interest:
# mu2, mutation rate per cell division before ZGA 
# mu1, mutation rate per cell division after ZGA
# s, stage at which cells are selected to contribute to ICM (number of cells in the embryo at that time = 2^s)
# n, number of cells forming the ICM 

# You can also add the following for your simulation
# n, number of cells at the time of embryo split
# a, asymmetry of split (proportion of cells contributing to twin 1)
# cell death (but this only enters after 64-cell stage) > for this, you can perhaps just assume a value or sth
# cell mixing after division (no assume none, and can think about this later)

# Assumptions
# ZGA occurs at the 4-cell stage 
# some specific parameter (or the way to select parameters) for mutation rate / how to choose it 

# Summary statistics
# number of mutations shared between twin 1 and twin 2 (+ VAF)
# number of mutations specific to twin 1 (+ VAF)
# number of mutations specific to twin 2 (+ VAF)

# what they did
# simulate the embryo until 1,024 cells 
# select n cells at stage s 
# simulations with ABC package
# posterior distribution for parameter estimates generated from the 1% best simulations
# determined by the Euclidian distance of the summary statistics from what was observed in the individual 

# Script from Tim's paper 
# https://github.com/TimCoorens/PanBody_Phylogenies/blob/main/Analysis/Embryonic_bottleneck_ABC.R

# ABC model to estimate embryonic bottleneck
# Tim Coorens - February 2021 

# Without death rate parameter

options(stringsAsFactors = F)
iter_num=as.numeric(commandArgs(T)[1])

prob_selection=c() # calculate vector of probabilities for cell selection 
for(stage in 2:6){ # check values from stage 2 to stage 6; stage 2 corresponds to 4 cells (so we are testing 4-64 cell stages)
  for(number in 2:(2^(stage)-1)){ # take the number (of cells I presume), from 2 to all but 1 cell available at that stage (ie at least one cell contributes to each TE / ICM)
    probs=c()
    for(k in 1:min(number,2^(stage-1))){ # from 1 to whatever number we took 
      probs=c(probs,binom.test(x=342,n=500,p = k/number)$p.value) # binomial test: nr of successes, nr of trials, 
    } # where on earth do the 342 and 500 come from???
    prob_selection=rbind(prob_selection,c(stage,number,max(probs),which.max(probs)))
  }
}

colnames(prob_selection)=c("stage","number","pval","picked")
prob_selection=as.data.frame(prob_selection)
prob_selection_final=prob_selection[prob_selection$pval>0.01,] # select only cases where the probability of selection > 0.01

odd <- function(x) x%%2 != 0 # return TRUE if x is odd and FALSE if x is even 
even <- function(x) x%%2 == 0 # return TRUE if x is even and FALSE if x is odd  

# create a symmetric, easy phylo tree for the stage (basically just codes the nr of cell divisions, really easy)
bi_tree=function(stage){ # create a tree as function of stage 
  vec=1:(2^stage)
  for(n in 1:stage){
    vec=paste0("(",vec[odd(1:length(vec))],
               ",",vec[even(1:length(vec))],")")
  }
  tree=read.tree(text=paste0(vec,";"))
  return(tree)
}

# if you want to visualize the output of this, you can do
plot.phylo(bi_tree(3)) # displays the phylogenetic tree for stage 3, or 8-cell embryo 

find_children = function(nodes, tree=tree_df){
  all_tips=c()
  for(node in nodes){
    child_nodes = tree$node[tree$parent==node]
    tips = c()
    for (a in 1:length(child_nodes)){
      if (tree$isTip[tree$node==child_nodes[a]]){
        tips=c(tips,child_nodes[a])
      }
      else{
        tips=c(tips,find_children(child_nodes[a],tree))
      }
    }
    all_tips=c(all_tips,tips)
  }
  return(all_tips)
}

find_nodes_gen=function(stage, tree){
  parent_nodes=Ntip(tree)+1
  for(s in 1:stage){
    parent_nodes=tree$edge[(tree$edge[,1]%in%parent_nodes),2]
  }
  return(parent_nodes)
}

# function to simulate a tree given parameters (stage to select cells to ICM, number of cells, mutation rate before ZGA, mutation rate after ZGA)
sim_tree=function(stage,number,mu1,mu2){
  tree=bi_tree(10) # get a tree until stage 10 (1,024 cells)
  tree$edge.length = rpois(lambda=mu2,n=nrow(tree$edge)) # sample n numbers from a Poisson distribution of parameter lambda 
  tree$edge.length[tree$edge[,2]%in%c(find_nodes_gen(1, tree=tree),find_nodes_gen(2, tree=tree))]=rpois(lambda=mu1,n=6)
  tree_df=as.data.frame(fortify(tree)) # fortify basically converts a graph to a dataframe 
  
  nodes_gen=find_nodes_gen(stage, tree=tree) # find nodes 
  drop_nodes=sample(nodes_gen,size=2^stage-number) # which nodes to get rid of 
  drop_tips=find_children(drop_nodes,tree_df) # remove children of nodes that were dropped  
  
  tree_subset=drop.tip(tree,tip=drop_tips) # drop childre from the tree
  tree_df=as.data.frame(fortify(tree_subset)) # create a new dataframe with some tips removed 
  
  orig_furcation=sum(tree_df$parent%in%which(tree_df$x==0)&tree_df$x>0) # I guess which cell / stage the bifurcation started in?
  multiforcation_score=sum(tree_df$branch.length[tree_df$x<=5&tree_df$node!=Ntip(tree)+1]<1)/sum(tree_df$x<=5&tree_df$node!=Ntip(tree)+1)
  if(is.nan(multiforcation_score))multiforcation_score=0
  d1=length(find_children(nodes=find_nodes_gen(1, tree=tree_subset)[1],tree=tree_df))
  asymmetry=max(d1,Ntip(tree_subset)-d1)/Ntip(tree_subset)
  burden_in_first_branches=tree_df$branch.length[tree_subset$edge[(tree_subset$edge[,1]==Ntip(tree_subset)+1),2]]
  return(c(orig_furcation,asymmetry,max(burden_in_first_branches),min(burden_in_first_branches),multiforcation_score, mu1, mu2, stage, number))
}

# what does this return 
sim_tree(2, 2, 1, 1)
# 2.0000000 0.5000000 5.0000000 1.0000000 0.5646259 1.0000000 1.0000000 2.0000000 2.0000000

sim_tree(6, 3, 1, 1) # you simulate a tree where you select 3 cells at the 64 cell stage to make your ICM, 1 mutation per cell division before and after ZGA 
# I don't get what the first output means 
# the second is the asymmetry you'd observe on the tree 
# next is the range of mutation burden (from max to min) on the early tree (before split)?

stage_number_mat=c()
for(stage in 2:6){
  num_range=2:min(30,max(2,2^(stage)-1)) # create a matrix for each stage and nr of cells (cap at 30)
  stage_number_mat=rbind(stage_number_mat,data.frame(stage=stage,number=num_range))
}

n_iter=1000
sim_results=matrix(0,ncol=9,nrow=n_iter) # save results for each of the simulations to a matrix 
for(k in 1:n_iter){ # for each iteration 
  mu1=sample(seq(0.5,6,by=0.01),1) # select mu1 from uniform distribution 
  mu2=sample(seq(0.1,1.5,by=0.01),1) # select mu2 from uniform distribution 
  select=sample(1:nrow(stage_number_mat),1) # randomly select a row of the matrix for stages and cells 
  stage=stage_number_mat$stage[select] # get the stage from this row  
  number=stage_number_mat$number[select] # get the number from this row 
  sim_results[k,]=sim_tree(mu1=mu1,mu2=mu2,stage=stage,number=number) # simulate tree for the numbers you've selected 
  
#  if (k%%1000==0){
#    print(k)}
  }


colnames(sim_results)=c("orig_furcation","asymmetry","max_burden_branch","min_burden_branch","multifurcation","mu1","mu2","stage","number")
sim_results=as.data.frame(sim_results)
# write.table(sim_results,paste0("sim_results_",iter_num,"_round2.txt"))


#--------
library(abc)
library(data.table)

# all_sim=fread("all_sim_results.txt",data.table=F) #Concatenate and read in results
select=all_sim_stage_number_vec%in%stage_number_mat_vec # okay this was never defined so not like I am going to find out what this is meant to be 
all_sim = sim_results
select=5 # I am just going to pick a random number, I don't know co autor mial na mysli w tym wypadku x

# get the stats for observed data (they had samples from 3 adults and so could compare the summary stats to the observed thing)
PD28690_stats=data.frame(orig_furcation=2,
                         asymmetry=0.685,
                         min_burden_branch=1,
                         max_burden_branch=5,
                         multifurcation=0.8064516)  

PD43851_stats=data.frame(orig_furcation=2,
                         asymmetry=0.93,
                         min_burden_branch=6,
                         max_burden_branch=6,
                         multifurcation=0)  

PD43850_stats=data.frame(orig_furcation=2,
                         asymmetry=0.60,
                         min_burden_branch=1,
                         max_burden_branch=2,
                         multifurcation=0.4210526)  

# some notes on ABC so you know what's going on 
# there are 3 different algorithms to construct the posterior distribution; rejection, regression-based correction (either local linear regression or neural networks)
# we compute summary stats from the simulation and compare them to the summary stats observed in the data (we do this using the Euclidean distance, d)
# if the distance is less than a threshold (tolerance!) that we set, we accept the set of parameters that got us the simulated summary stats 
# tol is the tolerance rate = percentage of accepted simulations

abc_results=abc(target=PD43850_stats, # vector of OBSERVED summary stats (your data)
                param=all_sim[1:1000,c("mu1","mu2","stage","number")], # these are your parameters with which you simulated the thing 
                sumstat=all_sim[1:1000,c("orig_furcation","asymmetry","min_burden_branch","max_burden_branch","multifurcation")], # these are the summary stats that you get if you run the simulation with those parameters 
                tol=0.01, hcorr=F,method="rejection") # hcorr = should you apply the conditional heterodcedastic model? > need to read up on this but the vignette does not explain why we would apply this as a default

# what does this return?
# gives you the number of accepted simulations and summary statistics for those 
# euclidian distance to observed simulations 

abc_results=abc(target=PD28690_stats, 
                param=all_sim[all_sim$stage<5&all_sim$number%in%c(3,6,9),c("mu1","mu2","death_rate")], 
                sumstat=all_sim[all_sim$stage<5&all_sim$number%in%c(3,6,9),c("orig_furcation","asymmetry","min_burden_branch","max_burden_branch","multifurcation")], 
                tol=0.01, hcorr=F,method="rejection")

abc_results=abc(target=PD43851_stats[,c("orig_furcation","asymmetry","min_burden_branch","max_burden_branch","multifurcation")], 
                param=all_sim[,c("mu1","mu2")], 
                sumstat=all_sim[,c("orig_furcation","asymmetry","min_burden_branch","max_burden_branch","multifurcation")], 
                tol=0.01, hcorr=T,method = "neuralnet")


stage_number_mat=c()
for(stage in 2:6){
  num_range=2:min(30,max(2,2^(stage)-1))
  stage_number_mat=rbind(stage_number_mat,data.frame(stage=stage,number=num_range))
}

cell_alloc_mat=matrix(ncol=6, nrow=30)

for(k in 1:nrow(stage_number_mat)){
  
  n=stage_number_mat[k,"number"]
  s=stage_number_mat[k,"stage"]
  
  cell_alloc_mat[n,s]=sum(abc_results$unadj.values[,c("number")]==n&abc_results$unadj.values[,"stage"]==s)

  }

library(RColorBrewer)
library(gplots)
pdf("PD43850_bottleneck.pdf",width=3.5,height=8)
heatmap.2(cell_alloc_mat,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          dendrogram="none",
          na.color = "gray",
          col=colorRampPalette(brewer.pal(9,"Blues"))(100),
          Rowv = F,
          Colv = F,
          scale='none',
          key=F,
          cexRow=0.8,
          cexCol=1,
          notecex=1.0,
          mar=c(10,6.5))
dev.off()


