library(raster)
library(rgeos )
library(rgdal)
library(ggplot2)
library(ggrepel)
library(tidyverse)

#SET PATH HERE
path_to_wd <- 
  
setwd(path_to_wd)

#Get shapefiles. Working Directory should have a single shapefile with terrains

#set the name of the shapefile here
shape_file <- 
  
shapes <-  readOGR(dsn = paste0(getwd(),'/shps' ), layer = "shape_file" , stringsAsFactors = FALSE )

#Set one or multiple methods here. Default All. Available options:
# SI: 'rabia', 'storie', 'sq'
# Methods: "Rabia", "Storie", "Square Root"
SI <- c('rabia', 'storie', 'sq')
methods <- c("Rabia", "Storie", "Square Root")

color_data <-  c("#d7191c", "#fdae61", "#ffffbf", "#a6d96a", "#1a9641")
labels <- c("  0-10: Permanently not suitable", 
            " 10-25: Marginally not suitable", 
            " 25-50: Marginally suitable",
            " 50-75: Moderately suitable",
            "75-100: Highly suitable")

#main loop. Ensure working directory has outputs and plots subdirectory
for (s in 1:length(SI)){
  
    SI_raster <- raster(paste0(getwd(), '/outputs/', SI[s],'_merged.tif'))
    
    for (sh in 1:length(shapes)){
      
        shape = shapes[sh,]
        suit_rast <- crop(mask(SI_raster, shape), shape)
        valid_cells <- sum ( !is.na(values(suit_rast) ))
        
        table_data <- as.numeric(table(cut(matrix(suit_rast),
                                 breaks = c(-1, 10, 25, 50, 75, 101),
                                 labels = c(0, 10, 25, 50, 75))))
        
        area <-  round( (table_data * (res(SI_raster)[1] * 111000)^2) / 10000, 2)
        percent <- round((table_data / valid_cells) * 100, 2)
        hsize <- 3
        comb_data <- data.frame(percent_data = percent,
                                area_data = area,
                                label_data = labels,
                                colors = color_data)
        
        comb_data <- subset(comb_data, percent_data != 0)
        
        
        comb_data <- comb_data  %>% mutate(x = hsize)
        
        posdata <- comb_data %>% mutate(   csum = rev(cumsum(rev(percent_data))),
                                           pos = percent_data/2 + lead(csum, 1),
                                           pos = ifelse(is.na(pos), percent_data/2, pos)
                                           )
        
        ggplot(comb_data, aes(x = hsize, y = percent_data, fill = label_data)) +
                 geom_col() + 
                 geom_label_repel(data = posdata, aes(y= pos, label = paste0(area_data, 'ha' ,
                                               '\n',"(", percent_data, ')%')),
                           nudge_x =0.5,
                           size =3,
                           fontface = "bold") +
                 coord_polar(theta = "y") +
                 scale_fill_manual(values = color_data) +
                 xlim(c(0.2,hsize+0.5))+
                 theme(panel.background = element_rect(fill = 'transparent'),
                       plot.background = element_rect(fill='transparent', color=NA),
                       panel.grid = element_blank(),
                       axis.title = element_blank(),
                       axis.ticks = element_blank(),
                       axis.text = element_blank()) +
                 guides(fill = guide_legend(title= "Suitability Score"))+
                 ggtitle(paste("Suitability Index Score Terrain:",shape$cc_code ), 
                         subtitle = paste(methods[s], "method"))+
                 theme(legend.position="none",
                       plot.title = element_text(hjust = 0.5, color='black', face='bold', size = 15),
                       plot.subtitle = element_text(hjust = 0.5, color = 'black', face='bold', size = 13))
        
       ggsave(paste0(getwd(),"/plots/", "plot_", shape$cc_code,"_", SI[s],".png"), bg='transparent')
    }
}


write.csv(RSITabData_area, "RSITabData_area3.csv")
write.csv(RSITabData_percent, "RSITabData_percent3.csv")
