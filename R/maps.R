# Ggplot data for plotting maps

## Plotting all predictions####

### Lgeo entire region####
#### Climate only ####
pred_cur <- rast('output/sdms/lgeo/pred_rasters/present_latrodectus_geometricus.tif')

my_palette <- colorRampPalette(c('#4575b4','#ffffbf','#d73027'))(100)
plot(pred_cur, col = my_palette)
plot(as(world, "Spatial"),add=T,col = NA, border = "black")
plot(as((states %>% filter(admin == 'United States of America')),'Spatial'), add=T,col=NA,border = "black")
points(data.frame(lgeo_filt$x,lgeo_filt$y), pch = 20, cex = 0.5, col = 'black')





plotdata <- data.frame(crds(pred_cur$max),as.data.frame(pred_cur$max))
# Without points
ggplot(data = plotdata) +
  geom_tile(aes(x=x,y=y,fill=max))+
  # coord_fixed()+
  scale_fill_gradient2(name = 'Predicted suitability',
                       low = '#4575b4',mid= '#ffffbf', high = '#d73027',
                       midpoint = min(plotdata$max)+((max(plotdata$max))-min(plotdata$max))/2)+
  geom_sf(data = world, fill=NA)+
  geom_sf(data = states, fill = NA)+
  coord_sf(xlim=c(min(plotdata$x),max(plotdata$x)),ylim=c(15,max(plotdata$y)))
theme(panel.grid.major = element_line(color = gray(.9),linetype = "dashed", size = 0),
      panel.background = element_rect(fill = "aliceblue"),
      plot.title = element_text(face = "italic"))
#coord_sf(xlim=c(min(plotdata$x),max(plotdata$x)),ylim=c(min(plotdata$y),max(plotdata$y)))
ggsave('output/sdms/lgeo/prediction_maps/present_latrodectus_geometricus_woPoints.png',
       width = 8, height = 8,dpi = 600)

# With points
ggplot(data = plotdata) +
  geom_tile(aes(x=x,y=y,fill=max))+
  #coord_fixed()+
  scale_fill_gradient2(name = 'Predicted suitability',
                       low = '#4575b4',mid= '#ffffbf', high = '#d73027',
                       midpoint = min(plotdata$max)+((max(plotdata$max))-min(plotdata$max))/2)+
  geom_sf(data = world, fill=NA)+
  geom_sf(data = states, fill = NA)+
  #Adding points
  geom_point(data=lgeo_filt,aes(x=x,y=y),size=0.5,alpha=0.5)+
  coord_sf(xlim=c(min(plotdata$x),max(plotdata$x)),ylim=c(15,max(plotdata$y)))+
  #coord_sf(xlim=c(-min(plotdata$x),max(plotdata$x)),ylim=c(min(plotdata$y),max(plotdata$y)))+
  theme(panel.grid.major = element_line(color = gray(.9),linetype = "dashed", size = 0),
        panel.background = element_rect(fill = "aliceblue"),
        plot.title = element_text(face = "italic"))
ggsave('output/sdms/lgeo/prediction_maps/present_latrodectus_geometricus_withPoints.png',
       width = 8, height = 8,dpi = 600)

## Focused on NE region
# Without points
ggplot(data = plotdata) +
  geom_tile(aes(x=x,y=y,fill=max))+
  # coord_fixed()+
  scale_fill_gradient2(name = 'Predicted suitability',
                       low = '#4575b4',mid= '#ffffbf', high = '#d73027',
                       midpoint = min(plotdata$max)+((max(plotdata$max))-min(plotdata$max))/2)+
  geom_sf(data = world, fill=NA)+
  geom_sf(data = states, fill = NA)+
  coord_sf(xlim=c(-86.155776,-67.793821),ylim=c(36.405098,45.062625))
theme(panel.grid.major = element_line(color = gray(.9),linetype = "dashed", size = 0),
      panel.background = element_rect(fill = "aliceblue"),
      plot.title = element_text(face = "italic"))
#coord_sf(xlim=c(min(plotdata$x),max(plotdata$x)),ylim=c(min(plotdata$y),max(plotdata$y)))
ggsave('output/sdms/lgeo/prediction_maps/present_latrodectus_geometricus_woPoints.png',
       width = 8, height = 8,dpi = 600)

#### Climate and urban ####
pred_cur_urb <- rast('output/sdms/lgeo/pred_rasters/present_latrodectus_geometricus_urb.tif')
plotdata <- data.frame(crds(pred_cur$max),as.data.frame(pred_cur$max))
ggplot(data = plotdata) +
  geom_tile(aes(x=x,y=y,fill=max))+
  # coord_fixed()+
  scale_fill_gradient2(name = 'Predicted suitability',
                       low = '#4575b4',mid= '#ffffbf', high = '#d73027',
                       midpoint = min(plotdata$max)+((max(plotdata$max))-min(plotdata$max))/2)+
  geom_sf(data = world, fill=NA)+
  geom_sf(data = states, fill = NA)+
  coord_sf(xlim=c(min(plotdata$x),max(plotdata$x)),ylim=c(15,max(plotdata$y)))
theme(panel.grid.major = element_line(color = gray(.9),linetype = "dashed", size = 0),
      panel.background = element_rect(fill = "aliceblue"),
      plot.title = element_text(face = "italic"))
#coord_sf(xlim=c(min(plotdata$x),max(plotdata$x)),ylim=c(min(plotdata$y),max(plotdata$y)))
ggsave('output/sdms/lgeo/prediction_maps/present_latrodectus_geometricus_urb.png',
       width = 8, height = 8,dpi = 600)

### Lgeo NE region####
#### Climate only ####

#### Climate and urban ####
### Native species####

### Overlaps####
lapply(list.files('output/sdms/overlap/rasters/',full.names = T, pattern = '.tif'),function(x){
  r <- rast(x)
  plotdata <- data.frame(crds(r),as.data.frame(r)) %>%
    rename(suit = 'max') %>% 
    mutate(suit = case_when(suit < (max(suit)-0.1) ~ NA,
                            .default = suit))
  #suit < (max(suit)-2) ~ 1))
  ggplot(data = plotdata) +
    geom_tile(aes(x=x,y=y,fill=suit))+
    coord_fixed()+
    scale_fill_gradient(name = 'Average suitability',
                        low = '#ffffbf', high = '#d73027')+
    #low = '#4575b4',mid= '#ffffbf', high = '#d73027',
    #midpoint = max(plotdata$max)-0.2)+
    theme(panel.grid.major = element_line(color = gray(.9),linetype = "dashed", size = 0),
          panel.background = element_rect(fill = "aliceblue"),
          plot.title = element_text(face = "italic"))
  #coord_sf(xlim=c(min(plotdata$x),max(plotdata$x)),ylim=c(min(plotdata$y),max(plotdata$y)))
  ggsave(str_replace(x,'rasters','maps') %>% str_replace('.tif','.png'),
         width = 8, height = 8,dpi = 600)
})
