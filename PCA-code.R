rm(list=ls(all=TRUE))
datadir <- "/Users/jarzyna.1/Documents/Dropbox/Manuscripts/2023_MEE_MAJ-JHS_PCA-Method/PCA-Biodiversity-Change/"
setwd(datadir)

#############################################################################
## Load required packages
#############################################################################
require(tidyverse)
require(here)
require(svglite)
require(viridis)
require(lubridate)
require(scales)
require(stars)
require(raster)
require(rnaturalearth)
require(grid)
require(maptools)
require(rgdal)
require(sf)
require(ggplot2)
require(dplyr)

select <- dplyr::select


##################################################
##################################################
##  Case study 1: Change in avian species richness and phylogenetic diversity (Breeding Bird Survey, BBS)
##################################################

##################################################
## Set initial values
##################################################
set.seed(7890)

##################################################
## Set projections
##################################################
### Original projection
crs_wgs84 <- st_crs(4326)

### US National Atlas Equal Area
new_proj = st_crs(2163)

### Create a bounding box to crop
bbox <- st_bbox(c(xmin = -120, ymin = 17, xmax = -78, ymax = 41), crs = new_proj)

##################################################
##  Read data in
##################################################
## Species richness (SR) and phylogneteic diversity (PD) data
## is organized as follows: rows are time steps (here, years, from 1969 through 2013)
## columns are values of biodiversity for each index (SR, PD) and site (BBS route)
data_bbs_mat <- readRDS(file = "data/cs1_data_mat.rds")
## Indexing data -- coordinates and BBS route id -- for each column in data_bbs_mat
index_bbs_df <- readRDS(file = "data/cs1_index_df.rds")

##################################################
##  Perform Principal Component Analysis, PCA
##################################################
### This line removes NA columns (there are none in our data)
na_col_test <- apply(is.na(data_bbs_mat),2,sum) > 0
data_bbs_mat <- data_bbs_mat[,!na_col_test]

##  Perform PCA with scaling using SVD
pca_fit <- prcomp(data_bbs_mat, retx=TRUE, scale=TRUE)

##  Show the values
pc_importance <- summary(pca_fit)
pc_importance


##################################################
##  Scree Plot
##################################################
## Create Scree plot
Scree <- data.frame(cbind(Component=seq(1,dim(pc_importance$importance)[2]),t(pc_importance$importance)))
Scree$EigenVal <- Scree$Standard.deviation^2
Scree_portion <- Scree[1:50,]

## Scree plot - proportion of variance explained
p <- ggplot(Scree_portion, aes(x = Component, y  = Proportion.of.Variance*100)) %>%
	+ geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
	+ scale_x_continuous(name = "Principal Component", limits=c(1,12),breaks=c(seq(1,12,1))) %>%
	+ scale_y_continuous(name = "Proportion of variance (%)", limits=c(0,13),breaks=c(0,4,8,12)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
        panel.border = element_rect(fill=NA, colour = "white", size=1),
        axis.line = element_line(color = 'black', size=1.5),
        plot.title = element_text(size=15, vjust=2, family="sans"),
        axis.text.x = element_text(colour='black',size=22),
        axis.text.y = element_text(colour='black',size=22),
        axis.title.x = element_text(colour='black',size=22),
        axis.title.y = element_text(colour='black',size=22),
        axis.ticks = element_line(color = 'black', size=1.5),
        axis.ticks.length=unit(0.3,"cm"),
        legend.position="none",
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p

p <- ggplot(Scree_portion, aes(x = Component, y  = Proportion.of.Variance*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,45),breaks=seq(1,45,10)) %>%
  + scale_y_continuous(name = "Proportion of variance (%)", limits=c(0,13),breaks=c(0,4,8,12)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
          panel.border = element_rect(fill=NA, colour = "white", size=1),
          axis.line = element_line(color = 'black', size=1.5),
          plot.title = element_text(size=15, vjust=2, family="sans"),
          axis.text.x = element_text(colour='black',size=22),
          axis.text.y = element_text(colour='black',size=22),
          axis.title.x = element_text(colour='black',size=22),
          axis.title.y = element_text(colour='black',size=22),
          axis.ticks = element_line(color = 'black', size=1.5),
          axis.ticks.length=unit(0.3,"cm"),
          legend.position="none",
          legend.text=element_text(size=20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
p

## plot of cumulative variance
p <- ggplot(Scree_portion, aes(x = Component, y  = Cumulative.Proportion*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,45),breaks=seq(1,45,10)) %>%
  + scale_y_continuous(name = "Proportion of variance (%)", limits=c(0,100),breaks=seq(0,100,20)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
          panel.border = element_rect(fill=NA, colour = "white", size=1),
          axis.line = element_line(color = 'black', size=1.5),
          plot.title = element_text(size=15, vjust=2, family="sans"),
          axis.text.x = element_text(colour='black',size=22),
          axis.text.y = element_text(colour='black',size=22),
          axis.title.x = element_text(colour='black',size=22),
          axis.title.y = element_text(colour='black',size=22),
          axis.ticks = element_line(color = 'black', size=1.5),
          axis.ticks.length=unit(0.3,"cm"),
          legend.position="none",
          legend.text=element_text(size=20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
p

p <- ggplot(Scree_portion, aes(x = Component, y  = Cumulative.Proportion*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,12),breaks=seq(1,12,1)) %>%
  + scale_y_continuous(name = "Proportion of variance (%)", limits=c(0,60),breaks=seq(0,60,10)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
          panel.border = element_rect(fill=NA, colour = "white", size=1),
          axis.line = element_line(color = 'black', size=1.5),
          plot.title = element_text(size=15, vjust=2, family="sans"),
          axis.text.x = element_text(colour='black',size=22),
          axis.text.y = element_text(colour='black',size=22),
          axis.title.x = element_text(colour='black',size=22),
          axis.title.y = element_text(colour='black',size=22),
          axis.ticks = element_line(color = 'black', size=1.5),
          axis.ticks.length=unit(0.3,"cm"),
          legend.position="none",
          legend.text=element_text(size=20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
p



###########################################################################
###  Download background data
###########################################################################
index_bbs_sf = st_as_sf(index_bbs_df, coords = c("Xcoord", "Ycoord"), crs = 4326 )

usa_states <- ne_states(country = 'United States of America', returnclass = 'sf')
usa_states <- fortify(usa_states)


##################################################
##  Plot Loadings Maps
##################################################
var_list <- c("pd", "sr")

### Loop through maps
for(i in seq(1,3)){
	### Create temporary loading dataframe
	loading_temp <- data.frame(index_uniq=rownames(pca_fit$rotation), load=pca_fit$rotation[,i])
	
  ### Create plotting dataframe
	plot_df <- index_bbs_sf %>%
		left_join(loading_temp, by = c("id" = "index_uniq"))

	### Find the limits
	plot_lims <- c(-1,1)*max(abs(plot_df$load), na.rm=TRUE)

	for (j in seq(1, 2)){
	var_j <- var_list[[j]]
  plot_j <- plot_df %>% filter(var == var_j)
  
  if (i == 2) {
    plot_j$load <- plot_j$load*(-1)
    ## We flipped scores and loadings for PC2 to ease interpretation. 
    ## Because scores and loadings are multiplied, flipping both signs does not change the meaning. 
    ## We chose positive loading to emphasize the strongest absolute loading values 
    ## and the most contiguous spatial region (western US in PC1, midwest and east in PC2).
  
	### Create plot
    p <- ggplot(data = plot_j) %>%
    	+ geom_sf(data = usa_states, fill = NA, alpha = 1, size=0.2) %>%
    	+ geom_sf(alpha = 0.95, aes(colour = load), size = 3) %>%
      #+ scale_colour_distiller(name = paste0("Loading\n", toupper(var_j), "\nPC ",i), palette = "RdBu" , limits = plot_lims) %>%
      + scale_colour_distiller(name = paste0("Loading\n", toupper(var_j), "\nPC ",i), palette = "RdBu" , limits = c(-0.08,0.08)) %>%
      #	+ coord_sf(xlim = c(-115, -100), ylim = c(29, 50), expand = FALSE) %>%
    	+ coord_sf(xlim = c(-27e5, 27e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
    	+ xlab("Longitude") %>%
    	+ ylab("Latitude") %>%
    	+ theme_bw(8) %>%
    	+ theme(legend.position="bottom") %>%
    	+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    #p
    
  } else {
    ### Create plot
    p <- ggplot(data = plot_j) %>%
      + geom_sf(data = usa_states, fill = NA, alpha = 1, size=0.2) %>%
      + geom_sf(alpha = 0.95, aes(colour = load), size = 3) %>%
      #+ scale_colour_distiller(name = paste0("Loading\n", toupper(var_j), "\nPC ",i), palette = "RdBu" , limits = plot_lims) %>%
      + scale_colour_distiller(name = paste0("Loading\n", toupper(var_j), "\nPC ",i), palette = "RdBu" , limits = c(-0.08,0.08)) %>%
      #	+ coord_sf(xlim = c(-115, -100), ylim = c(29, 50), expand = FALSE) %>%
      + coord_sf(xlim = c(-27e5, 27e5), ylim = c(1.7e5, 32e5), crs = 5070, expand = FALSE)  %>%
      + xlab("Longitude") %>%
      + ylab("Latitude") %>%
      + theme_bw(8) %>%
      + theme(legend.position="bottom") %>%
      +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    #p
  }
	### Save plot
	ggsave(p, file=paste0("graphics/",var_j, "_bbs_pc_",i,"_loading_map.png"), width=6, height=4, dpi=600)
	}
}


##################################################
##  Plot PC Scores
##################################################
## Loop through scores
for(i in seq(1, 3)){
	## Create scores
	plot_df <- data.frame(year = as.numeric(rownames(pca_fit$x)), score = pca_fit$x[,i])
	
	if (i == 2){
	plot_df$score <- plot_df$score*(-1)  
	## We flipped scores and loadings for PC2 to ease interpretation. 
	## Because scores and loadings are multiplied, flipping both signs does not change the meaning. 
	## We chose positive loading to emphasize the strongest absolute loading values 
	## and the most contiguous spatial region (western US in PC1, midwest and east in PC2).
	
	## Create plot
	p <- ggplot(plot_df, aes(x=year, y=score)) %>%
	  + geom_line(size=2) %>%
	  + geom_point(alpha=1,size=4, pch=16) %>%
	  + geom_hline(yintercept = 0, linetype = "dashed") %>%
	  + scale_x_continuous(name = "Time step (year)") %>%
	  + scale_y_continuous(name = "Score") %>%
	  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
	          panel.border = element_rect(fill=NA, colour = "white", size=1),
	          axis.line = element_line(color = 'black', size=1.5),
	          plot.title = element_text(size=15, vjust=2, family="sans"),
	          axis.text.x = element_text(colour='black',size=22),
	          axis.text.y = element_text(colour='black',size=22),
	          axis.title.x = element_text(colour='black',size=22),
	          axis.title.y = element_text(colour='black',size=22),
	          axis.ticks = element_line(color = 'black', size=1.5),
	          axis.ticks.length=unit(0.3,"cm"),
	          legend.position="none",
	          legend.text=element_text(size=20),
	          legend.title=element_blank(),
	          panel.grid.major = element_blank(), 
	          panel.grid.minor = element_blank())
	#p
	} else {
	  ## Create plot
	  p <- ggplot(plot_df, aes(x=year, y=score)) %>%
	    + geom_line(size=2) %>%
	    + geom_point(alpha=1,size=4, pch=16) %>%
	    + geom_hline(yintercept = 0, linetype = "dashed") %>%
	    + scale_x_continuous(name = "Time step (year)") %>%
	    + scale_y_continuous(name = "Score") %>%
	    + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
	            panel.border = element_rect(fill=NA, colour = "white", size=1),
	            axis.line = element_line(color = 'black', size=1.5),
	            plot.title = element_text(size=15, vjust=2, family="sans"),
	            axis.text.x = element_text(colour='black',size=22),
	            axis.text.y = element_text(colour='black',size=22),
	            axis.title.x = element_text(colour='black',size=22),
	            axis.title.y = element_text(colour='black',size=22),
	            axis.ticks = element_line(color = 'black', size=1.5),
	            axis.ticks.length=unit(0.3,"cm"),
	            legend.position="none",
	            legend.text=element_text(size=20),
	            legend.title=element_blank(),
	            panel.grid.major = element_blank(), 
	            panel.grid.minor = element_blank())
	}
	### Save plot
	ggsave(p, file=paste0("graphics/bbs_pc_",i,"_scores.png"), width=8, height=6, dpi=600)

}



##################################################
##################################################
##  Case study 2: Seasonal,variation in avian functional richness and functional evenness (eBird Status and Trends, S&Ts)
##################################################

###########################################################################
## Set initial values
###########################################################################
### Set seed so everyone gets the same random numbers (reproducible example)
set.seed(7890)


###########################################################################
## Set projections
###########################################################################
### Original projection
crs_sinu <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#new_proj = st_crs(5070)
### US National Atlas Equal Area
new_proj = st_crs(2163)

### Create a bounding box to crop
bbox <- st_bbox(c(xmin = -126, ymin = 17, xmax = -108, ymax = 41), crs = new_proj)

##################################################
##  Read data in: NOTE, DATA FOR eBIRD STATUS & TRENDS WILL BE UPLOADED UPON ACCEPTNCE OF THE MANUSCRIPT
##################################################
## Functional richness (FRic) and functional evenness (FEve) data
## is organized as follows: rows are time steps (here, weeks, from 1 through 52)
## columns are values of biodiversity for each index (FRic and FEve) and site (grid cell)
data_ebird_mat <- readRDS(file = "data/cs2_data_mat.rds")
## Indexing data -- coordinates and grid cell id -- for each column in data_bbs_mat
index_ebird_df <- readRDS(file = "data/cs2_index_df.rds")


##################################################
##  Perform Principal Component Analysis, PCA
##################################################
### This line removes NA columns (there are none in our data)
na_col_test <- apply(is.na(data_ebird_mat),2,sum) > 0
data_ebird_mat <- data_ebird_mat[,!na_col_test]

##  Perform PCA with scaling using SVD
pca_fit <- prcomp(data_ebird_mat, retx=TRUE, scale=TRUE)

##  Show the values
pc_importance <- summary(pca_fit)
pc_importance


###########################################################################
###  Download background data
###########################################################################
usa_states <- ne_states(country = 'United States of America', returnclass = 'sf')
usa_states <- fortify(usa_states)
crs_sinu <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
new_proj = st_crs(2163)

iso_list <- c("US-CA", "US-WA", "US-OR", "US-ID", "US-UT", "US-NV", "US-AZ")
region_test <- usa_states$iso_3166_2 %in% iso_list

region_states <- usa_states[region_test,]

### Merge into a single polygons
region_states <- region_states %>%
  st_transform(new_proj)


##################################################
##  Scree Plot
##################################################

### Create Scree plot
Scree <- data.frame(cbind(Component=seq(1,dim(pc_importance$importance)[2]),t(pc_importance$importance)))
Scree$EigenVal <- Scree$Standard.deviation^2
Scree_portion <- Scree[1:52,]

###
p <- ggplot(Scree_portion, aes(x = Component, y  = Proportion.of.Variance*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,52),breaks=c(seq(1,52,10))) %>%
  + scale_y_continuous(name = "Proportion of variance (%)", limits=c(0,20),breaks=c(0,5,10,15,20)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
          panel.border = element_rect(fill=NA, colour = "white", size=1),
          axis.line = element_line(color = 'black', size=1.5),
          plot.title = element_text(size=15, vjust=2, family="sans"),
          axis.text.x = element_text(colour='black',size=22),
          axis.text.y = element_text(colour='black',size=22),
          axis.title.x = element_text(colour='black',size=22),
          axis.title.y = element_text(colour='black',size=22),
          axis.ticks = element_line(color = 'black', size=1.5),
          axis.ticks.length=unit(0.3,"cm"),
          legend.position="none",
          legend.text=element_text(size=20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
p

p <- ggplot(Scree_portion, aes(x = Component, y  = Proportion.of.Variance*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,12),breaks=c(seq(1,12,1))) %>%
  + scale_y_continuous(name = "Proportion of variance (%)", limits=c(0,20),breaks=c(0,5,10,15,20)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
          panel.border = element_rect(fill=NA, colour = "white", size=1),
          axis.line = element_line(color = 'black', size=1.5),
          plot.title = element_text(size=15, vjust=2, family="sans"),
          axis.text.x = element_text(colour='black',size=22),
          axis.text.y = element_text(colour='black',size=22),
          axis.title.x = element_text(colour='black',size=22),
          axis.title.y = element_text(colour='black',size=22),
          axis.ticks = element_line(color = 'black', size=1.5),
          axis.ticks.length=unit(0.3,"cm"),
          legend.position="none",
          legend.text=element_text(size=20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
p


###
p <- ggplot(Scree_portion, aes(x = Component, y  = Cumulative.Proportion*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,52),breaks=seq(1,52,10)) %>%
  + scale_y_continuous(name = "Cumulative variance (%)", limits=c(0,100),breaks=seq(0,100,20)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
          panel.border = element_rect(fill=NA, colour = "white", size=1),
          axis.line = element_line(color = 'black', size=1.5),
          plot.title = element_text(size=15, vjust=2, family="sans"),
          axis.text.x = element_text(colour='black',size=22),
          axis.text.y = element_text(colour='black',size=22),
          axis.title.x = element_text(colour='black',size=22),
          axis.title.y = element_text(colour='black',size=22),
          axis.ticks = element_line(color = 'black', size=1.5),
          axis.ticks.length=unit(0.3,"cm"),
          legend.position="none",
          legend.text=element_text(size=20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
p


p <- ggplot(Scree_portion, aes(x = Component, y  = Cumulative.Proportion*100)) %>%
  + geom_line(size=2) %>%
  + geom_point(alpha=1,size=4, pch=16) %>%
  + scale_x_continuous(name = "Principal Component", limits=c(1,12),breaks=seq(1,12,1)) %>%
  + scale_y_continuous(name = "Cumulative variance (%)", limits=c(0,100),breaks=seq(0,100,10)) %>%
  + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
          panel.border = element_rect(fill=NA, colour = "white", size=1),
          axis.line = element_line(color = 'black', size=1.5),
          plot.title = element_text(size=15, vjust=2, family="sans"),
          axis.text.x = element_text(colour='black',size=22),
          axis.text.y = element_text(colour='black',size=22),
          axis.title.x = element_text(colour='black',size=22),
          axis.title.y = element_text(colour='black',size=22),
          axis.ticks = element_line(color = 'black', size=1.5),
          axis.ticks.length=unit(0.3,"cm"),
          legend.position="none",
          legend.text=element_text(size=20),
          legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
p


##################################################
##  Plot Loading Maps, WARNING--THIS IS A VERY COMPUTATIONALLY INTENSIVE STEP, BEST TO RUN ON A SUPERCOMPUTING PLATFORM
##################################################
var_list <- c( "feve", "fric")

### Loop through maps
for(i in seq(1,3)){
  
  #pc_folder <- file.path(write_path, paste0("pc_", i))
  #dir.create(pc_folder, recursive=TRUE, showWarnings = FALSE)
  
  ### Create temporary loading dataframe
  loading_temp <- data.frame(id=rownames(pca_fit$rotation), load=pca_fit$rotation[,i])
  
  loading_temp <- loading_temp %>%
    separate( col = "id", into = c("a", "b", "c"), sep = "_", remove = FALSE) %>%
    mutate(index_uniq = paste0(b, "_", c)) %>%
    mutate(var = a) %>%
    select(id, var, index_uniq, load)
  
  if (i == 2){
    loading_temp$load <- (loading_temp$load)*(-1)
    ## We flipped scores and loadings for PC2 to ease interpretation. 
    ## Because scores and loadings are multiplied, flipping both signs does not change the meaning. 
    ## We chose to flip the sign so that positive loading emphasizes peaks during the migrations seasons 
    ## to be consistent with PC3 (where positive loadings also are tied to peaks during autumn migration)
    
    ### Create plotting dataframe
    plot_df <- index_ebird_df %>%
      select(id, index_uniq, y_j, x_j, lon, lat) %>%
      full_join(loading_temp, by = c("id" = "id"))
    
    ### Find the limits
    plot_lims <- c(-1,1)*max(abs(plot_df$load), na.rm=TRUE)
    
    for (j in seq(1, length(var_list))){
      var_j <- var_list[[j]]
      ### Create raster
      raster_temp <- plot_df %>% filter(var == var_j) %>% select(lon, lat, load)
      raster_temp <- raster_temp %>% drop_na()
      raster_temp1 <- rasterFromXYZ(raster_temp, crs= crs_sinu)
      
      ### Convert to stars and reproject
      raster_temp <- st_as_stars(raster_temp1) %>%
        st_transform( new_proj)
      
      ### Crop the plot
      #raster_temp <- raster_temp[bbox]
      
      p <- ggplot() +
        geom_stars(data = raster_temp) +
        scale_fill_distiller(name = paste0("Loading\n", toupper(var_j), "\nPC ",i), palette = "RdBu" , limits = plot_lims) +
        geom_sf(data = region_states, fill = NA, alpha = 1, size=0.2) +
        xlab("Longitude") +
        ylab("Latitude") +
        theme_bw(8) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }
    
  } else {
    plot_df <- index_df %>%
      select(id, index_uniq, y_j, x_j, lon, lat) %>%
      full_join(loading_temp, by = c("id" = "id"))
    
    ### Find the limits
    plot_lims <- c(-1,1)*max(abs(plot_df$load), na.rm=TRUE)
    
    for (j in seq(1, length(var_list))){
      var_j <- var_list[[j]]
      ### Create raster
      raster_temp <- plot_df %>% filter(var == var_j) %>% select(lon, lat, load)
      raster_temp <- raster_temp %>% drop_na()
      raster_temp1 <- rasterFromXYZ(raster_temp, crs= crs_sinu)
      
      ### Convert to stars and reproject
      raster_temp <- st_as_stars(raster_temp1) %>%
        st_transform( new_proj)
      
      ### Crop the plot
      #raster_temp <- raster_temp[bbox]
      
      p <- ggplot() +
        geom_stars(data = raster_temp) +
        scale_fill_distiller(name = paste0("Loading\n", toupper(var_j), "\nPC ",i), palette = "RdBu" , limits = plot_lims) +
        geom_sf(data = region_states, fill = NA, alpha = 1, size=0.2) +
        xlab("Longitude") +
        ylab("Latitude") +
        theme_bw(8) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }

     ### Save plot
    ggsave(p, file=paste0("graphics/",var_j, "_ebird_pc_",i,"_loading_map.png"), width=6, height=4, dpi=600)
  }
}  


##################################################
##  Plot PC scores
##################################################

### Loop through scores
for(i in seq(1, 3)){
  
  ### Create loadings
  plot_df <- data.frame(week = seq(1,52), score = pca_fit$x[,i])
  
  if (i == 2){
    ## We flipped scores and loadings for PC2 to ease interpretation. 
    ## Because scores and loadings are multiplied, flipping both signs does not change the meaning. 
    ## We chose to flip the sign so that positive scores emphasizes peaks during the migration seasons 
    ## to be consistent with PC3 (where positive scores also are tied to peaks during autumn migration)
    
    p <- ggplot(plot_df, aes(x=week, y=score*(-1))) %>%
      + geom_line(size=2) %>%
      + geom_point(alpha=1,size=4, pch=16) %>%
      + geom_hline(yintercept = 0, linetype = "dashed") %>%
      + scale_x_continuous(name = "Time step (week)") %>%
      + scale_y_continuous(name = "Score") %>%
      + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
              panel.border = element_rect(fill=NA, colour = "white", size=1),
              axis.line = element_line(color = 'black', size=1.5),
              plot.title = element_text(size=15, vjust=2, family="sans"),
              axis.text.x = element_text(colour='black',size=22),
              axis.text.y = element_text(colour='black',size=22),
              axis.title.x = element_text(colour='black',size=22),
              axis.title.y = element_text(colour='black',size=22),
              axis.ticks = element_line(color = 'black', size=1.5),
              axis.ticks.length=unit(0.3,"cm"),
              legend.position="none",
              legend.text=element_text(size=20),
              legend.title=element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())
  } else {

  ### Create plot
  p <- ggplot(plot_df, aes(x=week, y=score)) %>%
    + geom_line(size=2) %>%
    + geom_point(alpha=1,size=4, pch=16) %>%
    + geom_hline(yintercept = 0, linetype = "dashed") %>%
    + scale_x_continuous(name = "Time step (week)") %>%
    + scale_y_continuous(name = "Score") %>%
    + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
            panel.border = element_rect(fill=NA, colour = "white", size=1),
            axis.line = element_line(color = 'black', size=1.5),
            plot.title = element_text(size=15, vjust=2, family="sans"),
            axis.text.x = element_text(colour='black',size=22),
            axis.text.y = element_text(colour='black',size=22),
            axis.title.x = element_text(colour='black',size=22),
            axis.title.y = element_text(colour='black',size=22),
            axis.ticks = element_line(color = 'black', size=1.5),
            axis.ticks.length=unit(0.3,"cm"),
            legend.position="none",
            legend.text=element_text(size=20),
            legend.title=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
  }
  
  ### Save plot
  ggsave(p, file=paste0("graphics/ebird_pc_",i,"_scores.png"), width=8, height=6, dpi=600)
  
}
