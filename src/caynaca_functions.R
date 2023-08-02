#------------------ UTILITY FUNCTIONS -------------------------------------------


#------------------ WORKFLOW FUNCTIONS -------------------------------------------
# in_sp_dt <- tar_read(rawdata_sp)
# in_env_dt <- tar_read(rawdata_env)

#--- format_spdata -------------
format_spdata <- function(in_sp_dt) {
  skim(in_sp_dt)
  
  spcols <- which(sapply(
    in_sp_dt[, -c('caso', 'total')], 
    is.numeric))
  
  in_sp_dt %>%
    setnames('estado_deFlujo', 'estado_de_flujo')
  
  in_sp_dt[, total_relative := 100*total/max(total), by=sitio] %>% #Compute relative total abundance (compared to max at that site)
    .[, taxo_richness := rowSums(.SD>0, na.rm=T), .SDcols=names(spcols)] #Compute taxonomic richess
  
  return(in_sp_dt)
  
}

#--- format envdata -------------
format_envdata <- function(in_env_dt) {
  skim(in_env_dt)
  
  in_env_dt[, unique(intermitencia)]
  
  
  in_env_dt[caudal %in% c(NA, 'SECO'), caudal := 0] %>% #Convert caudal/discharge to numeric
    .[, caudal := as.numeric(caudal)] %>%
    setnames(c('estado de Flujo', 'Bloque'), 
             c('estado_de_flujo', 'bloque')) %>%
    .[, #Fill NAs by site for altura, pendiente, orden (height, slope, and stream order) 
      c('altura', 'pendiente', 'orden') := .SD[
        !is.na(altura),][1, list(altura, pendiente, orden)], 
      by='sitio'] %>%
    .[estado_de_flujo == 'D', 
      c('AH_med', 'AH_max', 'prof_med', 'veloc_med') := 0] 
  
  in_env_dt[, cross_section_area := AH_med*prof_med] %>%
    .[, relative_cross_section_area := cross_section_area/max(cross_section_area,
                                                              na.rm=T), 
      by=sitio]
  
  #No site with more than 100% for summed portions of bed sediment granularity
  in_env_dt[, which(sum(c(bloque, piedra, grava, arena), na.rm=T) > 100),
            by=caso]
  
  #Did not record sediment when dry, so use average values from sites when flowing
  sed_cols <- c('bloque', 'piedra', 'grava', 'arena')
  in_env_dt[, 
            (sed_cols) :=
              sapply(.SD, simplify=F, function(x) {
                na.aggregate(x, FUN=mean)
              }), 
            by=sitio, 
            .SDcols = sed_cols]
  
  return(in_env_dt)
}

#--- plot_spdata -------------
plot_spdata <- function(in_spenv_dt, in_sp_dt) {
  
  #Boxplots of abundance in each site by species
  spcols <- which(sapply(
    in_sp_dt[, -c('caso', 'total')], 
    is.numeric))
  
  charcols <- which(sapply(in_sp_dt, is.character))
  
  plot_sp_abundance_1 <- 
    melt(in_sp_dt[, c(names(charcols),
                      names(spcols[1:(length(spcols)/2)]))
                  , with=F], 
         id.vars=c(names(charcols))
    ) %>%
    ggplot(aes(x=sitio, y=value)) +    
    geom_boxplot() +
    facet_wrap(~variable, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(size=7, angle=90))
  
  plot_sp_abundance_2 <-
    melt(in_sp_dt[, c(names(charcols),
                      names(spcols[((length(spcols)/2)+1):length(spcols)]))
                  , with=F], 
         id.vars=c(names(charcols))
    )%>%
    ggplot(aes(x=sitio, y=value)) +   
    geom_boxplot() +
    facet_wrap(~variable, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(size=7, angle=90))
  
  #Line plot of total abundance over time
  #by site, colored by site's intermittency status
  plot_abundance_vs_time <- ggplot(in_spenv_dt, aes(x=fecha, y=total, 
                                                    group=sitio)) + 
    geom_line(aes(color=intermitencia)) + 
    scale_y_sqrt() + 
    theme_bw()
  
  #Line plot of relative abundance over time (compared to max N at that site)
  #by site, colored by site's intermittency status
  plot_relative_abundance_vs_time <- ggplot(in_spenv_dt, 
                                            aes(x=fecha, y=total_relative)) + 
    geom_line(aes(color=intermitencia, group=sitio)) + 
    geom_smooth(aes(x=as.numeric(factor(fecha)), color=intermitencia), 
                linewidth=2, span=0.5, se=F
    ) +
    coord_cartesian(ylim=c(0,100), clip="off") + 
    theme_bw() 
  
  #Line plot of richness over time 
  #by site, colored by site's intermittency status
  plot_richness_vs_time <- ggplot(in_spenv_dt, 
                                  aes(x=fecha, y=taxo_richness)
  )+ 
    geom_line(aes(color=intermitencia, group=sitio)) + 
    geom_smooth(aes(x=as.numeric(factor(fecha)), color=intermitencia), 
                linewidth=2, span=0.5, se=F
    ) +
    coord_cartesian(ylim=c(0,40), expand=c(0,0), clip="off") + 
    theme_bw() 
  
  #Scatterplot of relative abundance against relative channel cross section area
  #across all sites and dates, colored by site's intermittency status
  plot_relative_abundance_vs_relative_area <- ggplot(
    in_spenv_dt, 
    aes(x=relative_cross_section_area,
        y=total_relative)) + 
    geom_point(aes(color=intermitencia)) +
    geom_quantile(aes(color=intermitencia), quantile=c(0.2,0.5,0.8)
                  #,linewidth=2, span=1, method='lm', se=F
    ) +
    coord_cartesian(ylim=c(0,100)) + 
    theme_bw() 
  
  #Scatterplot of abundance against channel cross section area
  #across all sites and dates, colored by site's intermittency status
  plot_abundance_vs_area <- ggplot(in_spenv_dt, 
                                   aes(x=cross_section_area,
                                       y=total)) + 
    
    geom_point(aes(color=intermitencia)) +
    geom_quantile(aes(color=intermitencia), quantile=c(0.2,0.5,0.8),
                  #,linewidth=2, span=1, method='lm', se=F
    ) +
    theme_bw() 
  
  
  return(list(
    plot_sp_abundance_1 = plot_sp_abundance_1,
    plot_sp_abundance_2 = plot_sp_abundance_2,
    plot_abundance_vs_time = plot_abundance_vs_time,
    plot_relative_abundance_vs_time = plot_relative_abundance_vs_time,
    plot_richness_vs_time = plot_richness_vs_time,
    plot_relative_abundance_vs_relative_area = plot_relative_abundance_vs_relative_area,
    plot_abundance_vs_area = plot_abundance_vs_area 
  ))
}

#--- plot_envdata -------------
plot_envdata <- function(in_env_dt) {
  
  env_melt <- melt(in_env_dt, 
                   id.vars=c('cuenca', 'presencia', 'sitio',
                             'fecha', 'intermitencia',
                             'estado_de_flujo'))
  
  env_boxplots_bysites <- ggplot(env_melt,
                         aes(x=sitio, y=value)) +
    geom_jitter(aes(color=estado_de_flujo)) +
    geom_boxplot(alpha=1/2) +
    facet_wrap(~variable, scales='free') +
    theme_bw() +
    theme(axis.text.x = element_text(size=7, angle=30))
  
  env_distribs_plot <- ggplot(env_melt,
                        aes(x=value)) +
    geom_histogram() +
    facet_wrap(~variable, scales='free') +
    theme_bw() 
  
  return(list(
    env_boxplots_bysites = env_boxplots_bysites,
    env_distribs_plot = env_distribs_plot
    )  
  )
}


#--- fix_net_dangles --------------
fix_net_dangles <- function(in_net, in_tolerance, IDcol='OBJECTID_1') {
  
  if (inherits(in_net, 'character')) {
    #Import segments for basin and project them using the projection of the points
    net <- sf::st_read(in_net, quiet=T)
  } else if (inherits(in_net, 'sf')) {
    net <- in_net
  }
  
  #Remove 0-length segments
  net_sub <- net[!(as.numeric(st_length(net)) == 0),]
  
  
  #Check disconnected reaches 
  net_int <- sf::st_intersection(net_sub, net_sub)
  
  p_net <- ggplot() +
    geom_sf(data = net)
  
  #Establish those to reconnect
  IDcoldupli <- paste0(IDcol, '.1')
  segments_to_correct <- net_sub[
    !(net_sub[[IDcol]] %in% 
        unique(
          net_int[net_int[[IDcol]] !=net_int[[IDcoldupli]],][[IDcol]])
    ),] 
  
  segments_to_snap_to <- net_sub[
    (net_sub[[IDcol]] %in% 
       unique(
         net_int[net_int[[IDcol]] !=net_int[[IDcoldupli]],][[IDcol]])
    ),] 
  
  #Correct network geometry by snapping the start and end points of all lines
  #to correct to nearest line with a tolerance
  #-- Get start and endpoints + join attributes back to the points
  segments_to_correct_extremities <- st_line_sample(segments_to_correct, 
                                                    sample=c(0,1))  %>%
    st_cast("POINT") %>% #Convert from MULTIPOINT (including both start and end points) to POINT
    cbind(
      segments_to_correct[rep(seq_len(nrow(segments_to_correct)), each=2),], #join attributes
      data.frame(pos = rep(c('start', 'end'), length(.)/2)) #add start and end attribute
    ) %>%
    .[,-(which(names(.) == 'geometry.1'))] #Remove LINESTRING attribute
  
  #Check segments to correct and their extremities
  # p_net +
  #   geom_sf(data=net_int) +
  #   geom_sf(data=segments_to_correct, color='red') +
  #   geom_sf(data = segments_to_correct_extremities, color='red')
  
  #Snap points
  #first computing a line between site and snapping place on nearest segment
  vect_extremities <- vect(segments_to_correct_extremities)
  
  sitesnap_l <- terra::nearest(vect_extremities,
                               vect(segments_to_snap_to), 
                               centroids = F, lines = T)
  sitesnap_l[[IDcol]] <- segments_to_correct_extremities[[IDcol]]
  sitesnap_l[['pos']] <- segments_to_correct_extremities[['pos']]
  sitesnap_l[['snap_dist']] <- perim(sitesnap_l) #Compute snapping distance
  
  #convert the line to a point (the line's end point)
  #Only keep snapped points that were moved within the snapping tolerance
  snapped_extremities <- sitesnap_l[sitesnap_l$snap_dist < in_tolerance &
                                      sitesnap_l$snap_dist > 0,] %>% 
    terra::as.points(.) %>%
    .[duplicated(paste0(.[[IDcol]][,1], .[['pos']][,1])),] 
  
  #Edit start and end vertices in network
  net_sub_edit <- copy(net_sub)
  
  for (i in 1:nrow(net_sub[net_sub[[IDcol]] %in% 
                           snapped_extremities[[IDcol]][,1], ])
  ) {
    line <- net_sub[net_sub[[IDcol]] %in% 
                      snapped_extremities[[IDcol]][,1], ][i,]
    #print(line[[IDcol]])
    
    new_node <- snapped_extremities[
      snapped_extremities[[IDcol]][,1] == line[[IDcol]], ]
    
    if (new_node$pos == 'start') {
      net_sub_edit[net_sub_edit[[IDcol]] %in% 
                     snapped_extremities[[IDcol]][,1], ][i,]$geometry[[1]][1,] <-
        as.vector(geom(new_node)[,c('x', 'y')])
    } else if (new_node$pos == 'end') {
      net_sub_edit[net_sub_edit[[IDcol]] %in% 
                     snapped_extremities[[IDcol]][,1], ][i,]$geometry[[1]][nrow(line$geometry[[1]]),]  <-
        as.vector(geom(new_node)[,c('x', 'y')])
    }
  }
  
  #Check intersections
  net_sub_edit_int <- sf::st_intersection(net_sub_edit, net_sub_edit)
  p_net +
    geom_sf(data=net_sub_edit_int) +
    geom_sf(data=segments_to_correct, color='red')
  #st_write(net_sub_edit, file.path(resdir, 'check.shp'), append=F)
  
  return(net_sub_edit)
}

#--- Snap sites to nearest segment ---------------------------------------------
# in_sites_point_path = tar_read(sites_shp_path)
# in_sitesSQL=""
# in_segments = tar_read(net_shp_path)
# out_path = sites_snapped_out
# proj = 'original'
# pt_IDcol = 'C___DIGO_N'
# nearby_ID = T
# seg_IDcol = 'OBJECTID_1'
# overwrite = F

#Custom snap method. Lighter and faster than other tests options
#sf::st_snap doesn't work
#maptools::snapPointsToLines requires SpatialPoint and SpatialLines - too heavy/slow for this dataset
#The option used here relies on terra::nearest, which is faster and returns only 
# one line for each point compared to sf::st_nearest_points
snap_sites <- function(in_sites_point_path, 
                       in_sitesSQL="", 
                       in_segments, 
                       out_path,
                       proj = 'custom', 
                       pt_IDcol = NULL,
                       nearby_ID = F,
                       seg_IDcol = NULL,
                       overwrite = F) {
  
  sitesp <- terra::vect(in_sites_point_path, query = in_sitesSQL)
  
  #Project sites 
  # Global datasets tend to be in geographic coordinates. The unit of these 
  # coordinates are decimal degrees, whose west-east length decreases
  # with increasing latitude. Therefore, for identifying the nearest line,
  # which is based on distance calculation, the point dataset needs to be 
  # projected. However, no single projection is valid for the entire planet. 
  # Consequently, for each basin, the sites are projected using a custom 
  # projection which minimizes distortions for distance calculations within the
  # sites' bounding box.
  if (proj == 'custom') {
    if (nrow(sitesp) > 1 & 
        ((xmax(sitesp) != xmin(sitesp)) | (ymax(sitesp) != ymin(sitesp)))
    ){
      sitesp_proj <- terra::project(sitesp, 
                                    dist_proj(sitesp))
    } else {
      #if only one site, project to UTM
      sitesp_proj <- terra::project(
        sitesp,
        paste0('+proj=utm +zone=', 
               floor((xmin(sitesp) + 180) / 6) + 1,
               ' +datum=WGS84 +units=m +no_defs +ellps=WGS84')
      )
      
    }
  } else if (proj == 'original') {
    sitesp_proj <- sitesp
  } else if (inherits(proj, 'numeric')) {
    sitesp_proj <- terra::project(sitesp, crs(paste0('EPSG:', proj)))
  }
  
  #Import segments for basin
  segs <- terra::vect(in_segments)
  
  
  #Project them using the projection of the points
  if (crs(segs) != crs(sitesp_proj)) {
    segs <- terra::project(segs, crs(sitesp_proj))
  }
  
  #Make sure extent is right (bug)
  actual_segs_ext <- t(
    as.data.table(geom(segs))[
      , list(min(x, na.rm=T), max(x, na.rm=T),
             min(y, na.rm=T), max(y, na.rm=T))])
  segs <- crop(segs, actual_segs_ext)
  
  #Check that segment has unique IDs
  if (sum(duplicated(segs[[seg_IDcol]])) > 0) {
    seg_IDcol <- 'new_segID'
    segs[[seg_IDcol]] <- seq_along(segs)
  } 
  
  #Snap points (fastest way in R, it seems):
  #first computing a line between site and snapping place on nearest segment
  sitesnap_l <- terra::nearest(sitesp_proj, segs, centroids = F, lines = T)
  sitesnap_l[[pt_IDcol]] <- sitesp_proj[[pt_IDcol]]
  
  #convert the line to a point (the line's end point)
  sitesnap_p <- terra::as.points(sitesnap_l) %>%
    .[duplicated(.[[pt_IDcol]]),]
  
  #Join ID of nearest line to that point
  if (nearby_ID == T) {
    sitesnap_p[[seg_IDcol]] <- terra::nearby(
      sitesnap_p, segs, k=1)[, 'k1'] %>%
      as.data.frame(segs)[., seg_IDcol] 
  }
  
  #Reproject points to WGS84
  if (proj == 'custom') {
    sitesnap_p <- terra::project(sitesnap_p, "+proj=longlat +datum=WGS84")
  }
  
  terra::writeVector(sitesnap_p,
                     out_path,
                     overwrite=overwrite)
  
  return(out_path)
}

#--- Compute geodesic distance -------------------------------------------------
compute_eucdist <- function(in_sites_path, IDcol) {
  sites <- terra::vect(in_sites_path)
  distmat <- terra::distance(x=sites, y=sites)
  rownames(distmat) <- colnames(distmat) <- sites[[IDcol]][[1]]
  
  return(distmat)
}


#--- Compute network distance --------------------------------------------------
# in_sites_snapped <- tar_read(sites_snapped_path)
# in_net <- tar_read(net_formatted)
# IDcol <- sites_IDcol

compute_netdist <- function(in_sites_snapped, 
                            in_net,
                            IDcol) {
  #Import sites
  if (inherits(in_sites_snapped, 'character')) {
    #Import segments for basin and project them using the projection of the points
    sites <- sf::st_read(in_sites_snapped, quiet=T)   
  } else if (inherits(in_sites_snapped, 'sf')) {
    sites <- in_sites_snapped
  }
  
  #Slightly round coordinates to make sure that edges are connected
  in_net <- sf::st_cast(in_net, "LINESTRING")
  st_geometry(in_net) = st_geometry(in_net) %>%
    lapply(function(x) round(x, 5)) %>%
    st_sfc(crs = st_crs(in_net))
  
  st_geometry(sites) = st_geometry(sites) %>%
    lapply(function(x) round(x, 5)) %>%
    st_sfc(crs = st_crs(sites))
  
  #Format sfnetwork to compute distance matrix
  sfnet <-in_net %>%
    sfnetworks::as_sfnetwork(directed=F) %>% #Convert segments to sfnetwork
    tidygraph::convert(to_spatial_simple) %>% #Remove pseudonodes
    tidygraph::convert(to_spatial_smooth)  %>% #Remove loops
    tidygraph::convert(to_spatial_subdivision) %>% #in sfnetwork, edges that aren’t connected at terminal nodes are considered disconnected. so deal with that
    activate("edges") %>%
    mutate(length = edge_length())#Re-compute edge length
  
  #Compute distance matrix (keeping only segments downstream of sites)
  dist_mat <- sfnetworks::st_network_blend(sfnet, sites) %>% #Append sites to sfnetwork
    sfnetworks::st_network_cost(from = sites, #Compute network distances
                                to = sites,
                                weights = 'length')
  
  rownames(dist_mat) <- sites[[IDcol]]
  colnames(dist_mat) <-sites[[IDcol]]
  
  #Check for dangles in sfnetwork topology
  # open_ended_nodes <- sfnet %>%
  #   activate("nodes") %>%
  #   mutate(degree = centrality_degree()) %>%
  #   st_as_sf() %>%
  #   mutate(row = row_number()) %>%
  #   filter(degree == 1)
  # 
  # disconnected_edges <- sfnet %>%
  #   activate("edges") %>%
  #   st_as_sf() %>%
  #   filter(from %in% open_ended_nodes$row | to %in% open_ended_nodes$row)
  # 
  # mapview(sfnet %>% activate("edges") %>% st_as_sf(), layer.name="rail network") +
  #   mapview(disconnected_edges, color="red", layer.name = "nodes with only 1 edge") +
  #   mapview(open_ended_nodes, color="red", col.regions="red", layer.name = "edges of 1 edge nodes") +
  #   mapview(sites, color='black', col.regions="black")
  
  return(dist_mat)
}


#--- Compute environmental distance--------------------------------------------
# in_env_dt <- tar_read(env_dt)
# IDcol = sites_IDcol

compute_envdist <- function(in_env_dt, IDcol) {
  #Isolate actual env attributes
  envcols  <- c('AH_max', 'pH', 'conductivity_esp', 'oxigen_sat', 'TDS',
                'AH_med', 'prof_med', 'veloc_med', 'altura', 'pendiente', 
                'orden', 'caudal', 'bloque', 'piedra', 'grava', 'arena')
  
  #Transform data to normal distribution
  in_env_dt[, `:=`(AH_max_sqrt = sqrt(AH_max),
                   AH_med_sqrt = sqrt(AH_med),
                   prof_med_sqrt = sqrt(prof_med),
                   veloc_med_sqrt = sqrt(veloc_med),
                   caudal_log10 = log10(caudal+0.01))]
  envcols_edit <- plyr::mapvalues(envcols, 
                                  c('AH_max', 'AH_med', 'prof_med', 
                                    'veloc_med', 'caudal'),
                                  c('AH_max_sqrt', 'AH_med_sqrt', 'prof_med_sqrt', 
                                    'veloc_med_sqrt', 'caudal_log10')
  )
  
  #z-standardise data
  env_dt_norm <- in_env_dt[, sapply(.SD, function(x) {
    scale(x, center=T, scale=T)
  }),
  .SDcols = envcols_edit] %>%
    as.data.table
  
  #Check correlation plot
  env_corr_plot <- ggcorrplot::ggcorrplot(
    round(
      cor(env_dt_norm, use='pairwise.complete.obs'), 
      2)
    )
  
  #Compute Gower's distance
  #With each category of variable weighted by 1
  # Hydraulic variables (5 including 2 nearly indentical for AH):
  #   'AH_max_sqrt', 'AH_med_sqrt' 
  #   'prof_med_sqrt', 'veloc_med_sqrt', 'caudal_log10'
  # Physiographic variables (3):
  #   'altura', 'pendiente', 'orden'
  # Physico-chem variables (4 including 2 nearly identical for conduct and TDS):
  #   'pH', 'conductivity_esp', 'oxigen_sat', 'TDS'
  # Sediments (4):
  #   'bloque', 'piedra', 'grava', 'arena'
  env_dist_weighted <- cluster::daisy(env_dt_norm, 
                                      metric = 'gower',
                                      weights = c(
                                        1/8, #AH_max_sqrt
                                        1/3, #pH
                                        1/6, #conductivity_esp
                                        1/3, #oxigen_sat
                                        1/6, #TDS
                                        1/8, #AH_med_sqrt
                                        1/4, #prof_med_sqrt
                                        1/4, #veloc_med_sqrt
                                        1/3, #altura
                                        1/3, #pendiente
                                        1/3, #orden
                                        1/4, #caudal_log10
                                        1/4, #bloque
                                        1/4, #piedra
                                        1/4, #grava
                                        1/4 #arena
                                      )
  )
  
  env_dist_unweighted <- cluster::daisy(env_dt_norm, 
                                      metric = 'gower',
                                      weights = c(
                                        1/2, #AH_max_sqrt
                                        1, #pH
                                        1/2, #conductivity_esp
                                        1, #oxigen_sat
                                        1/2, #TDS
                                        1/2, #AH_med_sqrt
                                        1, #prof_med_sqrt
                                        1, #veloc_med_sqrt
                                        1, #altura
                                        1, #pendiente
                                        1, #orden
                                        1, #caudal_log10
                                        1, #bloque
                                        1, #piedra
                                        1, #grava
                                        1 #arena
                                      )
  )
  
  return(list(
    env_dist_weighted,
    env_dist_unweighted
  ))
}

#----- Compute alpha diversity for each site and year -------------------------


#----- Compute spatial beta diversity for each year ---------------------------
#Pairwise Bray–Curtis dissimilarity distance

#----- Compute temporal beta diversity for each site --------------------------

#----- Compute Gamma diversity for each year ----------------------------------
#Total number of species
#Total number of genera
#Total number of family

#----- Compute nMDS trajectories for each site over time ----------------------
#nMDS across all sites and dates
#then plot
