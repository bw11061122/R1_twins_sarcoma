###################################################################################################################################
# SCRIPT 7

# 2024-11-23
# Barbara Walkowiak bw18

# Plotting nice graphs to show on phylogeny 

###################################################################################################################################
# LIBRARIES 

# Load needed libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(glue)
library(comprehenr) # list comprehension in R 
library(stringr)
library(beeswarm)
library(viridis)
library(grid)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

###################################################################################################################################
# INPUT FILES 

# Plotting trees for the phylogeny (semi-manual)


# can I get trees in R?
BiocManager::install("ggtree")
library(ggtree)

# okay try the tutorial thing https://4va.github.io/biodatasci/r-ggtree.html 

library(tidyverse)
library(ggtree)

tree <- read.tree("~/Desktop/1SB/data_tree/tree_newick.nwk")
tree
ggplot(tree) + geom_tree() + theme_tree() 

newick1 = "(A,B);"
newick2 = "((A,B), C);"
newick3 =  "((A,B), (C,D));"
newick4 =  "(((A,B), (C,D)), ((E,F), (G,H)));"

tree1 = read.tree(text = newick1)
tree2 = read.tree(text = newick2)
tree3 = read.tree(text = newick3)
tree4 = read.tree(text = newick4)

ggtree(tree1, layout = 'dendrogram')

ggtree(tree2, layout = 'dendrogram')
ggtree(tree3, layout = 'dendrogram')
ggtree(tree4, layout = 'dendrogram')

# Pie charts to put onto the tree 
create_pie_chart = function(percent, color, name){
  pie_data = data.frame(
    cat = c('shaded', 'empty'),
    value = c(percent, 100-percent)
  )
  ggplot(pie_data, aes(x='', y=value, fill=cat))+
    geom_bar(stat='identity', width=1)+
    coord_polar(theta='y')+
    scale_fill_manual(values=c('lightgrey', color))+
    theme_void()+
    theme(legend.position='none')
  ggsave(glue('Results/20241123_phyloPlot_pieChart_{name}.pdf'), height=2, width=2)
}

# percent = indicates what we think the contributing cell fraction was 
# note that for chr1_388 the VAF HAS TO be inflated but I honestly do not know why at this point
pie_chart_PD62341_1 = create_pie_chart(20, col_normal, 'PD62341_shared')
pie_chart_PD63383_1 = create_pie_chart(100, col_normal, 'PD63383_shared')

pie_chart_PD62341_2 = create_pie_chart(80, col_PD62341, 'PD62341_early')
pie_chart_PD63383_2 = create_pie_chart(0, col_PD63383, 'PD62341_early')

pie_chart_PD62341_3 = create_pie_chart(0, col_PD62341, 'PD63383_early')
pie_chart_PD63383_3.1 = create_pie_chart(40, col_PD63383, 'PD63383_early_1')
pie_chart_PD63383_3.2 = create_pie_chart(20, col_PD63383, 'PD63383_early_2')

create_pie_chart(100, col_PD62341, 'PD62341_100')
create_pie_chart(80, col_PD62341, 'PD62341_80')
create_pie_chart(40, col_PD62341, 'PD62341_40')
create_pie_chart(20, col_PD62341, 'PD62341_20')

create_pie_chart(100, col_normal, 'shared_100')
create_pie_chart(80, col_normal, 'shared_80')
create_pie_chart(50, col_normal, 'shared_40')
create_pie_chart(33, col_normal, 'shared_33')
create_pie_chart(67, col_normal, 'shared_67')
create_pie_chart(40, col_normal, 'shared_40')
create_pie_chart(20, col_normal, 'shared_20')
create_pie_chart(10, col_normal, 'shared_10')
create_pie_chart(53, col_normal, 'shared_53')

create_pie_chart(100, col_PD63383, 'PD63383_100')
create_pie_chart(80, col_PD63383, 'PD63383_80')
create_pie_chart(40, col_PD63383, 'PD63383_40')
create_pie_chart(20, col_PD63383, 'PD63383_20')
create_pie_chart(10, col_PD63383, 'PD63383_10')
create_pie_chart(53, col_PD63383, 'PD63383_53')

create_pie_chart_border = function(percent, color, col_border, name){
  pie_data = data.frame(
    cat = 'shaded',
    value = 100
  )
  ggplot(pie_data, aes(x='', y=value, fill=cat))+
    geom_bar(stat='identity', width=1, color = col_border, size =2.5)+
    coord_polar(theta='y')+
    scale_fill_manual(values=c(color))+
    theme_void()+
    theme(legend.position='none')
  ggsave(glue('Results/20241123_phyloPlot_pieChart_{name}_border.pdf'), height=2, width=2)
}
create_pie_chart_border(100, col_PD62341, col_tumour, 'tumourPD62341_100')
create_pie_chart_border(100, col_PD63383, col_tumour, 'tumourPD63383_100')
