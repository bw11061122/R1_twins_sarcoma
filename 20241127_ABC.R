####################################################################################################################################
# SCRIPT 7

# Getting started with ABC simulations
# 2024-11-27
# Barbara Walkowiak bw18

# INPUT: 

# OUTPUT:

# Script to run simulations

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

# Script from Tim's paper 
# https://github.com/TimCoorens/PanBody_Phylogenies/blob/main/Analysis/Embryonic_bottleneck_ABC.R

# ABC model to estimate embryonic bottleneck
# Tim Coorens - February 2021 

# Without death rate parameter

options(stringsAsFactors = F)

iter_num=as.numeric(commandArgs(T)[1])

prob_selection=c()
for(stage in 2:6){
  for(number in 2:(2^(stage)-1)){
    probs=c()
    for(k in 1:min(number,2^(stage-1))){
      probs=c(probs,binom.test(x=342,n=500,p = k/number)$p.value)
    }
    prob_selection=rbind(prob_selection,c(stage,number,max(probs),which.max(probs)))
  }
}
colnames(prob_selection)=c("stage","number","pval","picked")
prob_selection=as.data.frame(prob_selection)
prob_selection_final=prob_selection[prob_selection$pval>0.01,]

odd <- function(x) x%%2 != 0
even <- function(x) x%%2 == 0

bi_tree=function(stage){
  vec=1:(2^stage)
  for(n in 1:stage){
    vec=paste0("(",vec[odd(1:length(vec))],
               ",",vec[even(1:length(vec))],")")
  }
  tree=read.tree(text=paste0(vec,";"))
  return(tree)
}

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

sim_tree=function(stage,number,mu1,mu2){
  tree=bi_tree(10)
  tree$edge.length = rpois(lambda=mu2,n=nrow(tree$edge))
  tree$edge.length[tree$edge[,2]%in%c(find_nodes_gen(1, tree=tree),find_nodes_gen(2, tree=tree))]=rpois(lambda=mu1,n=6)
  tree_df=as.data.frame(fortify(tree))
  
  nodes_gen=find_nodes_gen(stage, tree=tree)
  drop_nodes=sample(nodes_gen,size=2^stage-number)
  drop_tips=find_children(drop_nodes,tree_df)
  
  tree_subset=drop.tip(tree,tip=drop_tips)
  tree_df=as.data.frame(fortify(tree_subset))
  
  orig_furcation=sum(tree_df$parent%in%which(tree_df$x==0)&tree_df$x>0)
  multiforcation_score=sum(tree_df$branch.length[tree_df$x<=5&tree_df$node!=Ntip(tree)+1]<1)/sum(tree_df$x<=5&tree_df$node!=Ntip(tree)+1)
  if(is.nan(multiforcation_score))multiforcation_score=0
  d1=length(find_children(nodes=find_nodes_gen(1, tree=tree_subset)[1],tree=tree_df))
  asymmetry=max(d1,Ntip(tree_subset)-d1)/Ntip(tree_subset)
  burden_in_first_branches=tree_df$branch.length[tree_subset$edge[(tree_subset$edge[,1]==Ntip(tree_subset)+1),2]]
  return(c(orig_furcation,asymmetry,max(burden_in_first_branches),min(burden_in_first_branches),multiforcation_score, mu1, mu2, stage, number))
}


stage_number_mat=c()
for(stage in 2:6){
  num_range=2:min(30,max(2,2^(stage)-1))
  stage_number_mat=rbind(stage_number_mat,data.frame(stage=stage,number=num_range))
}

n_iter=5000
sim_results=matrix(0,ncol=9,nrow=n_iter)
for(k in 1:n_iter){
  mu1=sample(seq(0.5,6,by=0.01),1)
  mu2=sample(seq(0.1,1.5,by=0.01),1)
  select=sample(1:nrow(stage_number_mat),1)
  stage=stage_number_mat$stage[select]
  number=stage_number_mat$number[select]
  sim_results[k,]=sim_tree(mu1=mu1,mu2=mu2,stage=stage,number=number)
  
  if (k%%1000==0){
    print(k)
  }
}
colnames(sim_results)=c("orig_furcation","asymmetry","max_burden_branch","min_burden_branch","multifurcation","mu1","mu2","stage","number")
sim_results=as.data.frame(sim_results)
write.table(sim_results,paste0("sim_results_",iter_num,"_round2.txt"))

#--------
# Model with death rate parameter

options(stringsAsFactors = F)
library(abc)
library(ape)
library(ggtree)
library(phytools)

iter_num=as.numeric(commandArgs(T)[1])

odd <- function(x) x%%2 != 0
even <- function(x) x%%2 == 0

bi_tree=function(stage){
  vec=1:(2^stage)
  for(n in 1:stage){
    vec=paste0("(",vec[odd(1:length(vec))],
               ",",vec[even(1:length(vec))],")")
  }
  tree=read.tree(text=paste0(vec,";"))
  return(tree)
}

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

sim_tree_death=function(stage,number,mu1,mu2,death_rate){
  tree=bi_tree(12)
  tree$edge.length = rpois(lambda=mu2,n=nrow(tree$edge))
  tree$edge.length[tree$edge[,2]%in%c(find_nodes_gen(1, tree=tree),find_nodes_gen(2, tree=tree))]=rpois(lambda=mu1,n=6)
  tree_df=as.data.frame(fortify(tree))
  
  nodes_gen=find_nodes_gen(stage, tree=tree)
  drop_nodes=sample(nodes_gen,size=2^stage-number)
  drop_tips=find_children(drop_nodes,tree=tree_df)
  
  tree_subset=drop.tip(tree,tip=drop_tips)
  tree_df=as.data.frame(fortify(tree_subset))
  
  death=tree_subset$edge[,2][which(rbinom(prob=death_rate,size=1,n=nrow(tree_subset$edge))==1)]
  dead_tips=unique(c(death[death<=Ntip(tree_subset)],find_children(death[death>Ntip(tree_subset)],tree_df)))
  while(length(dead_tips)==Ntip(tree_subset)){
    death=tree_subset$edge[,2][which(rbinom(prob=death_rate,size=1,n=nrow(tree_subset$edge))==1)]
    dead_tips=unique(c(death[death<=Ntip(tree_subset)],find_children(death[death>Ntip(tree_subset)],tree_df)))
  }
  
  tree_subset_death=drop.tip(tree_subset,tip=dead_tips)
  tree_df=as.data.frame(fortify(tree_subset_death))
  
  orig_furcation=sum(tree_df$parent%in%which(tree_df$x==0)&tree_df$x>0)
  
  multiforcation_score=sum(tree_df$branch.length[tree_df$x<=5&tree_df$node!=Ntip(tree_subset_death)+1]<1)/sum(tree_df$x<=5&tree_df$node!=Ntip(tree_subset_death)+1)
  if(is.nan(multiforcation_score))multiforcation_score=0
  d1=length(find_children(nodes=find_nodes_gen(1, tree=tree_subset_death)[1],tree=tree_df))
  asymmetry=max(d1,Ntip(tree_subset_death)-d1)/Ntip(tree_subset_death)
  burden_in_first_branches=tree_df$branch.length[tree_subset_death$edge[(tree_subset_death$edge[,1]==Ntip(tree_subset_death)+1),2]]
  return(c(orig_furcation,asymmetry,max(burden_in_first_branches),min(burden_in_first_branches),multiforcation_score, mu1, mu2, stage, number, death_rate))
}


stage_number_mat=c()
for(stage in 2:6){
  num_range=2:min(30,max(2,2^(stage)-1))
  stage_number_mat=rbind(stage_number_mat,data.frame(stage=stage,number=num_range))
}

n_iter=5000
sim_results=matrix(0,ncol=10,nrow=n_iter)
for(k in 1:n_iter){
  death_rate=sample(seq(0,0.10,by=0.001),1)
  mu1=sample(seq(0.5,6,by=0.01),1)
  mu2=sample(seq(0.1,1.5,by=0.01),1)
  select=sample(1:nrow(stage_number_mat),1)
  stage=stage_number_mat$stage[select]
  number=stage_number_mat$number[select]
  sim_results[k,]=sim_tree_death(mu1=mu1,mu2=mu2,stage=stage,number=number,death_rate=death_rate)
  
  if (k%%1000==0){
    print(k)
  }
}
colnames(sim_results)=c("orig_furcation","asymmetry","max_burden_branch","min_burden_branch","multifurcation","mu1","mu2","stage","number","death_rate")
sim_results=as.data.frame(sim_results)
write.table(sim_results,paste0("sim_death_results_",iter_num,"_round2.txt"))

#--------
library(abc)
library(data.table)
all_sim=fread("all_sim_results.txt",data.table=F) #Concatenate and read in results
select=all_sim_stage_number_vec%in%stage_number_mat_vec

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

abc_results=abc(target=PD43850_stats, 
                param=all_sim[select,c("mu1","mu2","stage","number","death_rate")], 
                sumstat=all_sim[select,c("orig_furcation","asymmetry","min_burden_branch","max_burden_branch","multifurcation")], 
                tol=0.01, hcorr=F,method="rejection")

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

for(k in 1:nrow(stage_number_mat)){
  n=stage_number_mat[k,"number"]
  s=stage_number_mat[k,"stage"]
  
  cell_alloc_mat[as.character(n),as.character(s)]=sum(abc_results$unadj.values[,c("number")]==n&
                                                        abc_results$unadj.values[,"stage"]==s)
}
colnames(cell_alloc_mat)=2^(2:6)

library(RColorBrewer)
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


# estimating crypt size 
options(stringsAsFactors = F)

# Load in observed data
physical_dists=readRDS("physical_dists_colon.Rdata")
colon_dists=readRDS("genetic_dists_colon.Rdata")

colon_dists_vec=unlist(colon_dists)
colon_dists_vec=as.numeric(colon_dists_vec[!is.na(colon_dists_vec)])

physical_dists_vec=unlist(physical_dists)
physical_dists_vec=as.numeric(physical_dists_vec[!is.na(physical_dists_vec)])

freq_dist=table(physical_dists_vec)

# Run the simulation
sim_iter=50000
limit=as.numeric(names(freq_dist[freq_dist>10])[sum(freq_dist>10)])

prob_vec_obs_all=prob_vec_obs_all_total=matrix(0,nrow=sim_iter,ncol=limit)
n_iter=sum(physical_dists_vec%in%c(1:limit))
radius_vec=rep(NA,sim_iter)
x_max=max(physical_dists_vec)

freqs=table(physical_dists_vec)
prob_vec=rep(0,x_max)
names(prob_vec)=1:x_max
prob_vec[names(freqs)]=freqs

for(k in 1:sim_iter){
  radius=sample(seq(2,20,by=0.01),1)
  dist_vec=in_circle_vec=rep(NA,n_iter)
  
  for(n in 1:n_iter){
    crypt1_x=sample(seq(-radius,radius,by=1),1)
    y_lim=sqrt(radius^2-crypt1_x^2)
    crypt1_y=sample(seq(-y_lim,y_lim,by=1),1)
    
    dist_vec[n]=dist=sample(1:limit,1,prob=prob_vec[1:limit])
    angle=sample(seq(1,2*pi,by=0.01),1)
    crypt2_x=crypt1_x+cos(angle)*dist
    crypt2_y=crypt1_y+sin(angle)*dist
    
    
    in_circle_vec[n]=sqrt(crypt2_x^2+crypt2_y^2)<=radius
  }
  freqs_sharing_obs=table(dist_vec[in_circle_vec])
  freqs_all_obs=table(dist_vec)
  prob_vec_obs=prob_vec_obs_total=rep(0,limit)
  names(prob_vec_obs)=names(prob_vec_obs_total)=1:limit
  prob_vec_obs[names(freqs_sharing_obs)]=freqs_sharing_obs
  prob_vec_obs_total[names(freqs_all_obs)]=freqs_all_obs
  prob_vec_obs_all[k,]=prob_vec_obs
  prob_vec_obs_all_total[k,]=prob_vec_obs_total
  
  radius_vec[k]=radius
  if (k%%100==0){
    print(k)
  }
}

freqs_sharing2=table(physical_dists_vec[colon_dists_vec>15]) #Select crypts that share more than 15 SNVs
prob_vec_sharing2=rep(0,x_max)
names(prob_vec_sharing2)=1:x_max
prob_vec_sharing2[names(freqs_sharing2)]=freqs_sharing2

prob_vec_obs_all=as.data.frame(prob_vec_obs_all)
prob_vec_obs_all_total=as.data.frame(prob_vec_obs_all_total)
colnames(prob_vec_obs_all)=colnames(prob_vec_obs_all_total)=1:limit

library(abc)
param_df=data.frame(radius=radius_vec)
abc_results=abc(target=prob_vec_sharing2[1:limit], 
                param=param_df, 
                sumstat=prob_vec_obs_all, 
                tol=0.05, hcorr=T,method="neuralnet")

#Generate images
pdf("colon_patch_abc_results.pdf",useDingbats = F)
plot(abc_results, param=param_df)
dev.off()

pdf("colon_sharing_hist.pdf",height=5,width=8)
hist(round(colon_dists_vec),breaks=40,col='steelblue',xlab="Number of shared SNVs",main="",probability = T)
lines(d,lty='dashed')
abline(v=15,lwd=2,col='red')
dev.off()

pdf("crypt_distance_genetic_all_cutoff.pdf")
image(k, col=r,xlab="Physical distance (in crypts)",ylab="Number of shared SNVs",xlim=c(min(k$x),40),ylim=c(min(k$y),80))
abline(h=15,lwd=2,col='red')
dev.off()
