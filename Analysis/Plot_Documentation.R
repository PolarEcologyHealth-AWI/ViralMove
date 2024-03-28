data_folder <- "/bioing/data/PathogenTransport/ViralMove_data"

load(file  = glue::glue("{data_folder}/Data/shp_breed_winter.RData"))
load(glue::glue("{data_folder}/Data/Map/eaafMap.rda"))
load(glue::glue("{data_folder}/Data/simulation_list.rda"))
load(glue::glue("{data_folder}/Data/Map/eaafMap.rda"))
load(glue::glue("{data_folder}/Data/mudflatTab.rda"))

###Distribution of breeding and wintering

spCols1 <- list(c("darkblue","cornflowerblue"), c("chartreuse4","chartreuse1"), c("brown","brown1"),c("darkgoldenrod3", "darkgoldenrod1" ))
names = c("Limosa lapponica", "Calidris canutus", "Calidris ferruginea","Calidris ruficollis")

plot_list <- lapply(1:4, function (x) {

plot <- ggplot() +
  geom_sf(data = eaafMap$grid, color = "grey80", lty = 3) +
  geom_sf(data = eaafMap$map, fill = "grey80") +
  geom_sf(data = tmp_shp$geometry[which(tmp_shp$SCINAME == names[x])], aes(fill = factor(tmp_shp$SEASONAL[which(tmp_shp$SCINAME == names[x])]))) +
  scale_fill_manual(name = ifelse(x ==1, "Season:",""), 
                    values = c("2" = spCols1[[x]][1], "3" = spCols1[[x]][2]),  
                    labels = c("2" = "Breeding", "3" = "Winter"))+
  labs(title = names(spParms)[x])+
  theme_void()+
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", legend.title.align = 0.5)
})


combined_plot <- plot_grid(plotlist = plot_list, ncol = 4)

ggsave("shp_plots.jpg", path = glue::glue("{data_folder}/Figures"), plot = combined_plot, width = 32, height = 13.5, bg = "white", units = "cm", dpi = 300)

####Combination Plot


load(glue::glue("{data_folder}/Data/simulation_list.RData"))

plot_list <- lapply(1:4, function (x) {
  filtered_data_winter <- simulation_list %>%
    filter(species == names(spParms)[x]) %>% group_by(winter_id) %>% summarize(n=sum(n))
  filtered_data_breed <- simulation_list %>%
    filter(species == names(spParms)[x]) %>% group_by(breeding_id) %>% summarize(n=sum(n))
  
  plot <- ggplot() +
    geom_sf(data = eaafMap$grid, color = "grey80", lty = 3) +
    geom_sf(data = eaafMap$map, fill = "grey80") +
    geom_sf(data = mudflatTab$geometry[filtered_data_winter$winter_id] %>% st_centroid(), color = spCols[x], alpha = 0.5, aes(size = filtered_data_winter$n)) +
    geom_sf(data = mudflatTab$geometry[filtered_data_breed$breeding_id] %>% st_centroid(), color = spCols[x], alpha = 0.8, aes(size = filtered_data_breed$n)) +
    scale_size_continuous(name = ifelse(x == 1,"Population Size:",""))+
    labs(title = names(spParms)[x])+
    theme_void()+
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", legend.title.align = 0.5)

})

library(cowplot)

combined_plot <- plot_grid(plotlist = plot_list, ncol = 4)

ggsave("combination_plots.jpg", path = glue::glue("{data_folder}/Figures"), plot = combined_plot, width = 32, height = 13.5, bg = "white", units = "cm", dpi = 300)

