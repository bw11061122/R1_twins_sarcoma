###################################################################################################################################
# SCRIPT 7

# November - December 2024
# Barbara Walkowiak bw18

# Plotting nice graphs to show on phylogeny
# Schematics to show VAF / cell fraction (pie chart / bar chart)

###################################################################################################################################
# LIBRARIES 

# Load needed libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(glue)
library(grid)
library(RColorBrewer)

###################################################################################################################################
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
  ggsave(glue('FiguresMain/Graphics/Graphics_phyloPlot_pieChart_{name}.pdf'), height=2, width=2)
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
create_pie_chart(0, col_PD62341, 'PD62341_0')

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
  ggsave(glue('FiguresMain/Graphics/Graphics_phyloPlot_pieChart_{name}_border.pdf'), height=2, width=2)
}
create_pie_chart_border(100, col_PD62341, col_tumour, 'tumourPD62341_100')
create_pie_chart_border(100, col_PD63383, col_tumour, 'tumourPD63383_100')

###################################################################################################################################
# BAR CHARTS 
# bar charts instead of pie charts

# Function to create stacked bar chart
create_barchart = function(color, percentage, name) {
  
  twin = strsplit(name, "_")[[1]][1]
  
  # dt
  data = data.frame(
    Category = factor(c("Colored", "Remaining"), levels = c("Remaining", "Colored")),
    Value = c(percentage, 100 - percentage)
  )
  
  # barchart
  ggplot(data, aes(x = "Bar", y = Value, fill = Category)) +
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_manual(values = c("Colored" = color, "Remaining" = "grey")) +
    labs(
      title = glue('{twin}'),
      x = NULL,
      y = "Percentage"
    ) +
    theme_classic(base_size=11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )

  ggsave(glue('FiguresMain/Graphics/Graphics_phyloPlot_barchart_{name}.pdf'), height=2, width=2)
}

create_barchart(col_PD62341, 100, 'PD62341_100')  
create_barchart(col_PD62341, 80, 'PD62341_80')  
create_barchart(col_PD62341, 40, 'PD62341_40')  
create_barchart(col_PD62341, 20, 'PD62341_20')  
create_barchart(col_PD62341, 0, 'PD62341_0')  

create_barchart(col_PD63383, 100, 'PD63383_100')  
create_barchart(col_PD63383, 80, 'PD63383_80')  
create_barchart(col_PD63383, 40, 'PD63383_40')  
create_barchart(col_PD63383, 20, 'PD63383_20')  
create_barchart(col_PD63383, 0, 'PD63383_0')  






