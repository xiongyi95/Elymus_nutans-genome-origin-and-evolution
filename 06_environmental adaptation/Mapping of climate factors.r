#!/usr/bin/env Rscript

library(tidyr)
library(ggplot2)
library("sf")
library(rgdal)
library(ggspatial)
library(gridExtra)
library(scales)
library(gstat)
library(sp)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(ncdf4)
library(scatterpie)


setwd("D:/test/haplotype/wc10_present")
tif_file <- "bio15.tif"
raster_bio1 <- raster(tif_file)
plot(raster_bio1)  #?????????????????????????????????
#bio1_color= c("#2892C7","#57A0BA","#78ADAC","#97BD9E","#B5CF8F","#CFDB8A","#E8DE82","#f7cb79","#F7B76D","#f59e5f","#F58653","#F7754D","#FF3333")
#bio1_color= c("#2892C7","#57A0BA","#78ADAC","#97BD9E","#B5CF8F","#CFDB8A","#E8DE82","#f7cb79","#F7B76D","#f59e5f","#F58653","#F7754D","white")
#bio1_color= c("#FF3333","#F7754D","#F58653","#f59e5f","#F7B76D","#f7cb79","#E8DE82","#CFDB8A","#B5CF8F","#97BD9E","#78ADAC","#57A0BA","#2892C7")
#names(bio1_color)=c("(-40)-(-20)","(-35)-(-30)","(-30)-(-25)","(-25)-(-20)","(-20)-(-15)","(-15)-(-10)","(-10)-(-5)","(-5)-0","0-5","5-10","10-15","15-20","20+")
country_shp=sf::st_read("./shapes/world.shp")
rasterbio1=stack("bio15.tif")
#creategroup <- function(gf){
#  colnames(gf)=c("x","y","bio1")
#  gf$level=ifelse(gf$bio1<=-20,"(-20)-",
#                  ifelse(gf$bio1<=-18,"(-20)-(-18)",
#                         ifelse(gf$bio1<=-16,"(-18)-(-16)",
#                                ifelse(gf$bio1<=-14,"(-16)-(-14)",
#                                       ifelse(gf$bio1<=-12,"(-14)-(-12)",
#                                              ifelse(gf$bio1<=-10,"(-12)-(-10)",
#                                                     ifelse(gf$bio1<=-8,"(-10)-(-8)",
#                                                            ifelse(gf$bio1<=-6,"(-8)-(-6)",
#                                                                   ifelse(gf$bio1<=-4,"(-6)-(-4)",
#                                                                          ifelse(gf$bio1<=-2,"(-4)-(-2)",
#                                                                                 ifelse(gf$bio1<=0,"(-2)-0",
#                                                                                        ifelse(gf$bio1<=20,"0-20","20+"
#                                                                                        ))))))))))))
#  return(gf)
#}
bio1_color= c("#2892C7","#57A0BA","#78ADAC","#97BD9E","#B5CF8F","#CFDB8A","#E8DE82","#f7cb79","#F7B76D","#f59e5f","#F58653","#F7754D","white")

names(bio1_color)=c("(-3)-","(-3)-(-2)","(-2)-(-1)","(-1)-1","1-3","3-5","5-7","7-9","9-11","11-13","13-15","15-17","17+")
creategroup <- function(gf){
  colnames(gf)=c("x","y","bio1")
  gf$level=ifelse(gf$bio1<=50,"(-3)-",
                  ifelse(gf$bio1<=65,"(-3)-(-2)",
                         ifelse(gf$bio1<=80,"(-2)-(-1)",
                                ifelse(gf$bio1<=95,"(-1)-1",
                                       ifelse(gf$bio1<=110,"1-3",
                                              ifelse(gf$bio1<=125,"3-5",
                                                     ifelse(gf$bio1<=140,"5-7",
                                                            ifelse(gf$bio1<=155,"7-9",
                                                                   ifelse(gf$bio1<=170,"9-11",
                                                                          ifelse(gf$bio1<=200,"11-13",
                                                                                 ifelse(gf$bio1<=250,"13-15",
                                                                                        ifelse(gf$bio1<=300,"15-17","17+"
                                                                                        ))))))))))))
  return(gf)
}
resdf<-as.data.frame(rasterbio1,xy=TRUE)%>%drop_na()
resdf=creategroup(resdf)


pie <- read.table("bio12_pie.txt",header = T)
pie$long <- as.numeric(pie$long)
pie$lati <- as.numeric(pie$lati)

p_bio1 = ggplot() +
  geom_sf(fill="#DADADA", data=country_shp) +
  geom_raster(aes(x = x, y = y, fill = factor(level)), data = resdf, stat = "identity") +
  geom_sf(fill="transparent", data=country_shp, size=0.5) +
  coord_sf(xlim = c(70,105), ylim = c(27,45)) +
  scale_fill_manual(values = c("#2892C7","#57A0BA","#78ADAC","#97BD9E","#B5CF8F","#CFDB8A","#E8DE82","#f7cb79","#F7B76D","#f59e5f","#F58653","#F7754D","white","#31646c","#dcdfd2"), name = "Level") +
  labs(x="Longitude", y="Latitude") +
  theme_bw() +
  theme(
    text = element_text(family="serif"),
    axis.ticks.length = unit(0.25, "lines"),
    axis.ticks = element_line(colour="black", size = 0.6),
    axis.text.x = element_text(size=12, colour = "black"),
    axis.text.y = element_text(size=12, colour = "black"),
    plot.title = element_text(size = 15L, hjust = 0),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.key.size = unit(0.1, 'inch'),
    legend.title = element_text(size=10.5, color="black", vjust=0.5, hjust=0.5),
    legend.position = "right",
    legend.background = element_rect(colour= "grey", fill= "white", size=0.6),
    legend.text = element_text(size=7.3, color="black", vjust=0.5, hjust=0.5),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.minor = element_line(colour = "white", size=0.1, linetype = 4),
    plot.margin = unit(c(0.1,0.1,0.1,0.1), "mm")
  ) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location="bl", style = north_arrow_nautical(fill = c("grey40","white"), line_col = "grey20")) +
  geom_scatterpie(data=pie, aes(x=long, y=lati, group=pop, r=0.5), cols=names(pie)[4:5])

ggsave(p_bio1, file="bio15.pdf", height=8, width=15)

