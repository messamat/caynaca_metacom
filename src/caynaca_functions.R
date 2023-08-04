#------------------ UTILITY FUNCTIONS -------------------------------------------
#---- box.cox.chord ---------------------------------------------------
#' Compute the box.cox.chord transformation on quantitative community composition 
#' data for any exponent. Usual exponents are larger than or equal to 0.
#'
#' Arguments --
#' @param mat : matrix or data.frame of quantitative non-negative community 
#'    composition data (frequencies, biomasses, energy measures, etc.)
#' @param bc.exp : Box-Cox exponent to the data before chord transformation. 
#'    Usual exponent values are {1, 0.5, 0.25, 0}, where 
#'    bc.exp=1: no transformation; 
#'    bc.exp=0.5: square-root transformation; 
#'    bc.exp=0.25: fourth-root (or double square-root) transformation; 
#'    bc.exp=0: log(y+1) transformation (default value). 
#'    Default value: bc.exp=0 (log(y+1) transformation).
#'
#' Value --
#' A Box-Cox+chord transformed matrix of the same size as the original data matrix.
#'
#' Author:: Pierre Legendre (Legendre and Brocard 2018)
#' License: GPL (>=2)

box.cox.chord <- 
  function(mat, 
           bc.exp=0) 
  { 
    # Internal function
    vec.norm <- function(vec)  sqrt(sum(vec^2))
    #
    chck <- apply(mat, 1, sum)
    if(any(chck == 0)) stop("Rows",which(chck==0)," of the data matrix sum to 0")
    #
    # Apply the user-selected Box-Cox exponent (bc.exp) to the frequency data
    if(bc.exp==0) {
      tmp <- log(mat+1) 
    } else { 
      tmp <- mat^bc.exp 
    }
    row.norms <- apply(tmp, 1, vec.norm)
    #
    # Apply the chord transformation to matrix "tmp" before returning it
    res <- sweep(tmp, 1, row.norms, "/")
  }

#---- BCD ---------------------------------------------------
#' Box-Cox transformation: find the best exponent to reach multivariate normality.
#'
#' Box-Cox-Dagnelie (BCD) method – Transform the data using different exponents. 
#' Default: exponents in the [0,1] interval by steps of 0.1. Test the multivariate 
#' normality of the data after each transformation using function dagnelie.test() 
#' from package ade4.
#' Note: the Dagnelie test requires that n > (rank+1) where 'n' is the number of 
#' obsevations and 'rank' is the rank of the covariance matrix.
#'
#' Arguments --
#' @param mat  Multivariate data table (object class: matrix or data.frame).
#' @param bc.exp  vector of exponents from the Box-Cox series for transformation, 
#'    for example bc.exp = c(0,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
#'    Positive and negative exponents are allowed by the function, although 
#'    it is recommended to use exponent values in the range [0,1].
#' @param chord  Chord-transform the data and recompute the Dagnelie test of 
#'    normality. Default: chord=TRUE; if chord=FALSE, do not transform data and 
#'    do not recompute the test.
#'
#' Value --
#' A table showing the Box-Cox exponent in the first column of each row. In 
#' columns 2 and 3, one finds the Shapiro-Wilk W statistic (BC_W) of the Dagnelie 
#' test of multivariate normality and the associated p-value (BC_p-val) after the 
#' exponent has been applied to the original data. Columns 4 and 5 show the same  
#' statistics(BC.chord_W and BC.chord_p-val) after the chord transformation has    
#' been applied to the Box-Cox transformed data.
#' 
#' References --
#'  Dagnelie, P. 1975. L'analyse statistique a plusieurs variables. 
#'  Les Presses agronomiques de Gembloux, Gembloux, Belgium.
#'
#'  Legendre, P. and L. Legendre. 2012. Numerical ecology, 3rd English
#'  edition. Elsevier Science BV, Amsterdam, The Netherlands.
#'
#' Author  Pierre Legendre
#' License GPL (>=2)

BCD <- 
  function(mat, 
           bc.exp=c(0,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
           chord=TRUE)
  {
    # Internal function
    vec.norm <- function(vec)  sqrt(sum(vec^2))
    #
    require(ade4)
    epsilon <- sqrt(.Machine$double.eps)
    mat <- as.matrix(mat)
    n <- nrow(mat)
    p <- ncol(mat)
    n.exp <- length(bc.exp)
    #
    if(chord) {
      res <- matrix(NA,n.exp,5)
      colnames(res) <- c("BC.exp","BC_W","BC_p-val","BC.chord_W","BC.chord_p-val")
    } else {
      res <- matrix(NA,n.exp,3)
      colnames(res) <- c("BC.exp", "BC_W", "BC_p-val")	
    }
    res[,1] <- bc.exp
    #
    if(any(mat < 0)) stop("Negative values not allowed in community data", 
                          call. = FALSE)
    
    chck1 <- apply(mat, 1, sum)
    if(any(chck1 == 0)) stop("One or several rows of 'mat' sum to 0", 
                             call. = FALSE)
    
    chck2 <- apply(mat, 2, var)
    keep.spec <- which(chck2 > epsilon)
    if(length(keep.spec) < p) {
      cat(length(keep.spec),"species have variances > 0 and were kept\n")
      cat("Species",which(chck2 <= epsilon)," were excluded\n")
      mat2 <- mat[,keep.spec] 
    } else { mat2 <- mat }
    #
    for(k in 1:n.exp) {
      if(bc.exp[k]==0) {
        # If BC exponent = 0, compute log(x+1)
        # Add 1 to the data before log transformation
        tmp <- log(mat2+1)                
        # Add 1 to the data before applying a negative exponent
      } else if(bc.exp[k]<0) { tmp <- (mat2+1)^bc.exp[k]
      # No transformation when bc.exp=1
      } else if(bc.exp[k]==1) { tmp <- mat2
      # Apply the exponent to the data
      } else { tmp <- mat2^bc.exp[k] }
      #
      tmp2 <- dagnelie.test(tmp)
      if((max(tmp2$D)-min(tmp2$D)) < epsilon)
        stop("All D values are equal, Dagnelie's test cannot be computed. ",
             "Check the data.", call. = FALSE)
      res[k,2] <- tmp2$Shapiro.Wilk$statistic
      res[k,3] <- tmp2$Shapiro.Wilk$p.value
      if(chord) {
        # Apply the chord transformation to matrix "tmp"
        row.norms <- apply(tmp, 1, vec.norm)
        mat3 <- sweep(tmp, 1, row.norms, "/")
        tmp2 <- dagnelie.test(mat3)
        res[k,4] <- tmp2$Shapiro.Wilk$statistic
        res[k,5] <- tmp2$Shapiro.Wilk$p.value
      }
    }
    res
  }

#------------------ WORKFLOW FUNCTIONS -------------------------------------------
# in_sp_dt <- tar_read(rawdata_sp)
# in_env_dt <- tar_read(rawdata_env)

#--- format_spdata -------------
format_spdata <- function(in_sp_dt) {
  skim(in_sp_dt)
  
  spcols <- names(
    which(sapply(
      in_sp_dt[, -c('caso', 'total')], 
      is.numeric))
  )
  
  in_sp_dt %>%
    setnames('estado_deFlujo', 'estado_de_flujo')
  
  in_sp_dt[, total_relative := 100*total/max(total), by=sitio] %>% #Compute relative total abundance (compared to max at that site)
    .[, taxo_richness := rowSums(.SD>0, na.rm=T), 
      .SDcols=spcols] #Compute taxonomic richess
  
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
  
  estado_de_flujo_ts <- in_env_dt[, list(
    estado_de_flujo = .SD[, .N, by=estado_de_flujo]$estado_de_flujo,
    rel_freq = .SD[, .N, by=estado_de_flujo]$N/.N), 
    by=fecha] %>%
    merge(expand.grid(c(paste0('fecha', seq(1,6))),
                      c('F', 'IP', 'D')),
          by.x=c('fecha', 'estado_de_flujo'),
          by.y=c('Var1', 'Var2'),
          all.y=T
    ) %>%
    .[is.na(rel_freq), rel_freq := 0]
          
  estado_de_flujo_tsplot <- ggplot(data=estado_de_flujo_ts) +
    geom_area(aes(x=fecha, y=rel_freq, fill=estado_de_flujo, group=estado_de_flujo),
              position='stack') +
    coord_cartesian(expand=c(0,0)) +
    theme_classic()
  
  return(list(
    env_boxplots_bysites = env_boxplots_bysites,
    env_distribs_plot = env_distribs_plot,
    estado_de_flujo_ts = estado_de_flujo_ts,
    estado_de_flujo_tsplot = estado_de_flujo_tsplot
  ))
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

#----- Compute spatial beta diversity for each date ---------------------------
#in_sp_dt <- tar_read(sp_dt)
compute_spatial_beta <- function(in_sp_dt) {
  #----- Prepare data ----------------------------------------------------------
  #Shroeder and Jenkins for choice of indices
  #Add some based on Anderson et al. 2011
  spcols <- names(
    which(sapply(
      in_sp_dt[, -c('caso', 'total', 'total_relative', 'taxo_richness')], 
      is.numeric))
  )
  
  # ggplot(data=melt(in_sp_dt[, spcols, with=F])) +
  #   geom_histogram(aes(x=value)) +
  #   scale_x_sqrt()
  
  #----- Jaccard Index for presence-absence by fecha ---------------------------
  presabs_dt <- in_sp_dt[, lapply(.SD, function(x) {as.numeric(x > 0)}),
                         .SDcols= spcols] %>%
    cbind(in_sp_dt[,c('caso', 'fecha'), with=F])
  
  jaccard_mats <- lapply(unique(presabs_dt$fecha), function(t) {
    subdt <- presabs_dt[fecha == t,]
    jaccard_dist <- vegan::vegdist(subdt[fecha == t, spcols, with=F],
                                   method = 'jaccard', diag=T)
    jaccard_mat <- as.matrix(jaccard_dist)
    jaccard_mat[upper.tri(jaccard_mat)] <- NA
    colnames(jaccard_mat) <- subdt$caso
    rownames(jaccard_mat) <- subdt$caso
    
    return(jaccard_mat)
  })
  names(jaccard_mats) <- unique(presabs_dt$fecha)
  
  #----- pres-abs-based total beta diversity, turnover and nestedness by fecha -----
  beta_jaccard <- lapply(unique(presabs_dt$fecha), function(t) {
    subdt <- presabs_dt[fecha == t,]
    sub_dist_jac <- betapart.core(subdt[, spcols, with=F]) %>%
      beta.multi(index.family="jac")
    return(sub_dist_jac)
  }) %>%
    rbindlist %>%
    .[, fecha := unique(presabs_dt$fecha)]
  
  #----- Bray–Curtis dissimilarity distance (untransformed) by fecha -----------
  bray_mats <- lapply(unique(in_sp_dt$fecha), function(t) {
    subdt <- in_sp_dt[fecha == t,]
    bray_dist <- vegan::vegdist(subdt[fecha == t, spcols, with=F],
                                method = 'bray', diag=T)
    bray_mat <- as.matrix(bray_dist)
    bray_mat[upper.tri(bray_mat)] <- NA
    colnames(bray_mat) <- subdt$caso
    rownames(bray_mat) <- subdt$caso
    
    return(bray_mat)
  })
  names(bray_mats) <- unique(in_sp_dt$fecha)
  
  #----- Abundance-based total beta diversity, turnover and nestedness -----
  beta_bray <- lapply(unique(in_sp_dt$fecha), function(t) {
    #fourth-root transform
    
    subdt <- in_sp_dt[fecha == t,]
    sub_dist_bray <- betapart.core.abund(subdt[, spcols, with=F]) %>%
      beta.multi.abund(index.family="bray")
    return(sub_dist_bray)
  }) %>%
    rbindlist %>%
    .[, fecha := unique(in_sp_dt$fecha)]

  #----- nMDS with Bray-Curtis -------------------------------------------------
  nmds_bray <- metaMDS(as.matrix(in_sp_dt[total !=0, spcols, with=F]),
                       distance = "bray", trymax=200, autotransform=T)
  
  nmds_bray_dt <- cbind(in_sp_dt[total !=0, -spcols, with=F],
                        nmds_bray$points
  )
  
  (nmds_bray_dt)
  
  #----- Chord dissimilarity matrix across all sites and dates -----------------
  #Determine Best Box-Cox+chord transformation coefficient
  best_bc_chord_exponent <- BCD(in_sp_dt[total!=0, spcols, with=F]) %>%
    as.data.table %>%
    .[which.max(BC.chord_W), BC.exp]
  
  #Compute transformation
  sp_chord <- box.cox.chord(mat = in_sp_dt[total!=0, spcols, with=F], 
                            bc.exp=best_bc_chord_exponent) 
  
  #Check output
  # ggplot(data=melt(sp_chord)) +
  #   geom_histogram(aes(x=value)) 
  
  #Compute box-chord distance
  sp_chord_dist <- dist(sp_chord) %>% 
    as.matrix
  rownames(sp_chord_dist) <- in_sp_dt[total!=0, caso]
  colnames(sp_chord_dist) <- in_sp_dt[total!=0, caso]
  
  #----- Chord dissimilarity matrix across all sites and dates -----------------
  chord_abundance_pca <- prcomp(x= sp_chord)
  screeplot(chord_abundance_pca)
  summary(chord_abundance_pca) #First 2 PCs capture only 22% of variance, first 4 capture 36% -- weak PCA
  
  pca_chord_dt <- cbind(
    in_sp_dt[total!=0, -spcols, with=F],
    scores(chord_abundance_pca)[, c('PC1', 'PC2', 'PC3')]
  )
  
  #----- Extra stuff ---------------------------------------
  # %>%
  #   merge(in_sp_dt[, c('caso', 'cuenca', 'sitio', 'fecha', 'intermitencia',
  #                      'estado_de_flujo', 'total', 'taxo_richness'),
  #                  with = F],
  #         by.x = 'caso.x', by.y='caso') %>%
  #   merge(in_sp_dt[, c('caso', 'cuenca', 'sitio', 'fecha', 'intermitencia',
  #                      'estado_de_flujo', 'total', 'taxo_richness'),
  #                  with = F],
  #         by.x = 'caso.y', by.y='caso')
  
  return(list(
    jaccard_mats = jaccard_mats,
    beta_jaccard = beta_jaccard,
    bray_mats = bray_mats,
    beta_bray = beta_bray,
    nmds_bray_dt =  nmds_bray_dt,
    sp_chord_dist = sp_chord_dist,
    chord_abundance_pca = chord_abundance_pca,
    pca_chord_dt = pca_chord_dt
  ))
} 
#----- Plot spatial beta div --------------------------------------------------
#in_spatial_beta <- tar_read(spatial_beta)

plot_spatial_beta <- function(in_spatial_beta) {
  
  #Plot beta diversity and its components over time based on pres-abs Jaccard
  beta_jaccard_melt <- in_spatial_beta$beta_jaccard %>%
    melt(id.vars='fecha')
  
  plot_beta_jaccard_fecha <- ggplot() +
    geom_area(data = beta_jaccard_melt[variable != 'beta.JAC'],
              aes(x=fecha, y=value, 
                  position = 'stack',
                  fill=variable, group=variable)) +
    geom_line(data = beta_jaccard_melt[variable == 'beta.JAC'],
              aes(x=fecha, y=value, group=variable, color=variable),
              size=2) + 
    coord_cartesian(expand=c(0,0)) +
    theme_classic() + 
    theme(legend.title = element_blank())
  
  #Plot beta diversity and its components over time based on abundance Bray-Curtis
  beta_bray_melt <- in_spatial_beta$beta_bray %>%
    melt(id.vars='fecha')
  
  plot_beta_bray_fecha <- ggplot() +
    geom_area(data = beta_bray_melt[variable != 'beta.BRAY'],
              aes(x=fecha, y=value, 
                  position = 'stack',
                  fill=variable, group=variable)) +
    geom_line(data = beta_bray_melt[variable == 'beta.BRAY'],
              aes(x=fecha, y=value, group=variable, color=variable),
              size=2) + 
    coord_cartesian(expand=c(0,0)) +
    theme_classic() + 
    theme(legend.title = element_blank())
  
  #Plot trajectories
  #MDS for Bray-Curtis
  ggplot(data=in_spatial_beta$nmds_bray_dt[intermitencia=='isolpool',],
         aes(x=MDS1, y=MDS2, color=estado_de_flujo)) +
    geom_point() +
    geom_line(aes(group=sitio))
  
  plot_nMDS_time <- ggplot(data=in_spatial_beta$nmds_bray_dt,
                           aes(x=MDS1, y=MDS2, 
                               shape=estado_de_flujo, color=intermitencia)) +
    geom_point(size=4, alpha=0.85) +
    facet_wrap(~fecha) +
    theme_classic()
  
  
}



#----- Compute nMDS trajectories for each site over time ----------------------
#nMDS across all sites and dates
#then plot


#----- Compute temporal beta diversity for each site --------------------------

#----- Compute Gamma diversity for each year ----------------------------------
#Total number of species
#Total number of genera
#Total number of family

