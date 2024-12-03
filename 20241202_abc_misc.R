
# simulate a tree until 11 generations
tree=make_tree(6) # create a tree (until 10 cell divisions) 
node_depths = node.depth(tree) # absolute depths 
nodes_to_icm = which(node_depths == 2^(6-4)) # 16 cells to choose from
cells_to_icm = sample(nodes_to_icm, 3) # select 3 cells from 16 to make ICM

# check that you select nodes at the right level 
plot.phylo(tree)
nodelabels(cells_to_icm, cells_to_icm)

# the logic here is that we set the boundary of 100 ICM cells to split into twins 
# the min nr of cells allocated to the ICM is 3 (selected at stage 4); the max nr of cells allocated to the ICM is 30 (at stage 6)
# in either case, at stage 11, we should have at least 96 cells (if 3 selected at stage 6) so we can do the split from the tree we simulated to stage 11 in any case
# maybe an alternative is to assume that split occurred at day 7 latest so we have max 2^8 cells (1 cell div per day + some sort of wiggle in case of faster div or sth)

tree$edge.length = rpois(lambda=mu2,n=nrow(tree$edge)) # add mutations onto each branch of the tree (sample from Poisson distribution with lambda = mutation rate after ZGA)
tree$edge.length[tree$edge[,2]%in%c(find_nodes_gen(1, tree=tree),find_nodes_gen(2, tree=tree))]=rpois(lambda=mu1,n=6) 
# ZGA at the 8 cell stage, so for the first 6 branches (in stage 1 and stage 2), replace with mutations generated with mu1 (basically higher, so we can have more mutations at this stage)
tree_df=as.data.frame(fortify(tree)) # convert to a dt so this is nice to work with

# identify all nodes at a given level (the level you want to take cells to the ICM)
node_depths = node.depth(tree) # absolute depths 
nodes_to_icm = which(node_depths == 2^(6-4)) # 16 cells to choose from
cells_to_icm = sample(nodes_to_icm, 3) # select 3 cells from 16 to make ICM

# I want to retain only ancestors and descendants of the selected nodes 
plot.phylo(tree)
nodelabels(cells_to_icm, cells_to_icm)

nodes_gen=find_nodes_gen(4, tree=tree) # find nodes of the stage of interest (say we take cells to ICM at stage 4 ie 16 cells) 
drop_nodes=sample(nodes_gen,size=2^4-3) # drop nodes at this stage (all cells minus nr of cells allocated to the LCM), say you have 16 cells and take 3 to the ICM
drop_tips=find_children(drop_nodes,tree_df) # drop tips of the nodes that you dropped 

# get the LCM tree only 
tree_subset=drop.tip(tree,tip=drop_tips) # create the subsetted tree
tree_df=as.data.frame(fortify(tree_subset)) # dt of the current tree

# okay now decide at what stage you want twins to split
# let's say that we selected 3 cells at stage 4 to ICM
# we let the ICM cells divide 3 more times each so we have 3 * 2^4 = 24 cells
# once we get 24 cells, we split the embryo into twins 
cells_to_split=find_nodes_gen(4, tree=tree_subset)
print(length(cells_to_split))
drop_tips=find_children(cells_to_split,tree_df)
# drop any cells after this stage
tree_sub=drop.tip(tree,tip=drop_tips)

# okay this is not working yet but what i want
# simulate the tree and selection of cells to the ICM
# let ICM cells divide a few more times
# subset a tree so you only have ancestors of ICM cells and children of those
# select when to split the tree
# select nodes (with specific asymmetry)
# calculate nr and VAF of twin-specific mutations and shared mutations

# then can simluate sequencing or sth like this 






