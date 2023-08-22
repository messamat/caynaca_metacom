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

#---- direct network ----------------------------------------------------------
# in_net <- tar_read(net_formatted)
# idcol <- 'OBJECTID_1'
# outlet_id <- 245

direct_network <- function(in_net,
                           idcol,
                           outlet_id
) {
  #------------------ Split lines at intersections -----------------------------
  st_precision(in_net) <- 0.05 #Reduce precision to make up for imperfect geometry alignments
  
  #Remove artefacts in network
  in_net <- in_net[in_net$length > 0.1 &
                     in_net[[idcol]] != 73,] 
  
  #Get outlet
  outlet_p <-  st_cast(in_net[in_net[[idcol]] == outlet_id,], "POINT") %>%
    .[nrow(.),]
  
  sfnet <- as_sfnetwork(in_net) %>%
    activate(edges) %>%
    arrange(edge_length()) %>%
    tidygraph::convert(to_spatial_simple) %>% #Remove loops
    tidygraph::convert(to_spatial_smooth) %>% #Remove pseudo nodes (doesn't work well)
    tidygraph::convert(to_spatial_subdivision) #Split at intersections
  
  net<- activate(sfnet, "edges") %>% #Grab edges
    st_as_sf() 
  net$newID <- seq_len(nrow(net)) #Create new IDs because of merging and resplitting
  
  #Get newID for outlet
  outlet_newID <- sf::st_intersection(net,
                                      outlet_p)[['newID']]
  
  #------------------ Identify dangle points -----------------------------------
  search_dangles <- function(in_network, idcol) {
    if (nrow(in_network) > 1) {
      #Get sfnetwork
      sfnet <- as_sfnetwork(in_network) %>%
        activate("edges") %>%
        arrange(edge_length()) %>%
        tidygraph::convert(to_spatial_subdivision) 
      
      edges <- activate(sfnet, "edges") %>%
        st_as_sf()
      nodes <- activate(sfnet, "nodes") %>%
        st_as_sf()                                
      
      #Intersect nodes and edges
      nodes_edges_inters <- sf::st_intersection(edges,
                                                nodes) %>%
        as.data.table
      
      #Identify nodes that intersect with only one edge (dangle points)
      dangle_points <- nodes_edges_inters[
        !(duplicated(nodes_edges_inters$.tidygraph_node_index) |
            duplicated(nodes_edges_inters$.tidygraph_node_index, fromLast = T)),
      ]
      
      #Plot network with dangle points
      netp <- ggplot() +
        geom_sf(data = edges, linewidth=1.2, color='blue') +
        geom_sf(data=  st_as_sf(dangle_points))
      
      print(netp)
      
      #Return newID for dangle points
      return(dangle_points[[idcol]])
    } else {
      return(NULL)
    }
  }
  
  #Set up loop that will iteratively identify dangle points, remove the associated
  #lowest-order edges from network, then re-identify dangle points, removing the 
  #associated lowest-order edges from network, and so on, iteratively, until
  #only the outlet edge remains
  
  dangles <- 'go!'
  order <- 1
  net_list <- list() #List into which each subsequent set of edges will be written
  
  while (length(dangles) > 0) {
    dangles <- search_dangles(in_network=net, idcol='newID') #identify dangle points
    
    dangle_segs_boolean <- (net[['newID']] %in% dangles &
                              net[['newID']] != outlet_newID)
    
    net[dangle_segs_boolean, 'stream_order'] <- order #assign the associated edges a stream order
    net_list[[order]] <-  net[dangle_segs_boolean,] #Write these edges to the list
    net <- net[!(dangle_segs_boolean),] #Remove these edges from the network
    order <- order + 1 
    
    net  <- as_sfnetwork(net) %>% 
      activate(edges) %>%
      arrange(edge_length()) %>%
      tidygraph::convert(to_spatial_smooth) %>% #Re-dissolve edges
      tidygraph::convert(to_spatial_simple) %>%
      st_as_sf()
    net$newID <- seq_len(nrow(net)) #Re-assign new IDs
    
    outlet_newID <- sf::st_intersection(net, #Re-dentify ID for outlet
                                        outlet_p)[['newID']]
    print(outlet_newID)
    
    st_precision(net) <- 0.05 #Make sure the correct precision is set up so that even edges that don't perfectly match can be linked
    
    write_sf(net[, c(idcol, 'newID')], #Write out intermediate network layer
             file.path(resdir, paste0('check', order, '.shp')),
             overwrite = T
    )
  }
  #Add the outlet edge
  net_list[[order]] <- net
  net_list[[order]]$stream_order <- order
  
  #Create a single sf
  out_net <- do.call(rbind, lapply(net_list, st_sf))
  
  return(out_net)
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
    .[, alpha_div := rowSums(.SD>0, na.rm=T), 
      .SDcols=spcols] #Compute taxonomic richess
  
  in_sp_dt[, gamma_div := melt(.SD)[value > 0,][!duplicated(variable), .N],
           .SDcols = spcols,
           by=fecha]
  
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
  plot_richness_vs_time <- ggplot(in_sp_dt, 
                                  aes(x=fecha)
  )+ 
    geom_line(data=in_sp_dt[!duplicated(fecha),],
              aes(x=fecha, y=gamma_div, group=cuenca, 
                  lty=cuenca), color='black',
              linewidth=1.5) +
    geom_line(aes(y=alpha_div, color=intermitencia, group=sitio)) + 
    geom_smooth(aes(x=as.numeric(factor(fecha)), y=alpha_div,
                    color=intermitencia), 
                linewidth=2, span=0.5, se=F
    ) +
    scale_y_continuous(name='Diversity') +
    scale_linetype(name='', labels=c('Gamma diversity')) +
    scale_color_manual(
      name = stringr::str_wrap('Alpha diversity by long-term flow regime', 30),
      values = c('#d73027', '#fdb863', '#4575b4'),
      labels = c('Intermittent: dry',
                 'Intermittent: disconnected pools',
                 'Perennial')) +
    guides(linetype = guide_legend(order = 1), 
           colour = guide_legend(order = 2)) +
    coord_cartesian(ylim=c(0,65), expand=c(0,0), clip="off") + 
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
  
  return(list(
    xy = geom(sites)[, c('x','y')],
    distmat = distmat
  ))
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
    sf::st_sfc(crs = st_crs(in_net))
  
  st_geometry(sites) = st_geometry(sites) %>%
    lapply(function(x) round(x, 5)) %>%
    sf::st_sfc(crs = st_crs(sites))
  
  #Format sfnetwork to compute distance matrix
  sfnet <-in_net %>%
    sfnetworks::as_sfnetwork(directed=F) %>% #Convert segments to sfnetwork
    tidygraph::convert(to_spatial_simple) %>% #Remove pseudonodes
    tidygraph::convert(to_spatial_smooth)  %>% #Remove loops
    tidygraph::convert(to_spatial_subdivision) %>% #in sfnetwork, edges that aren’t connected at terminal nodes are considered disconnected. so deal with that
    activate("edges") %>%
    dplyr::mutate(length = edge_length())#Re-compute edge length
  
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
  
  env_corr_plot <- ggcorrplot::ggcorrplot(
    round(
      cor(in_env_dt[, envcols_edit, with=F],
          use='pairwise.complete.obs'), 
      2)
  )
  
  
  #z-standardise data
  env_dt_norm <- in_env_dt[, sapply(.SD, function(x) {
    scale(x, center=T, scale=T)
  }),
  .SDcols = envcols_edit] %>%
    as.data.table %>%
    cbind(fecha=in_env_dt$fecha)
  
  #Check correlation plot
  
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
  
  full_weights <- c(
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
  
  env_dist_weighted <- cluster::daisy(env_dt_norm[, envcols_edit, with=F], 
                                      metric = 'gower',
                                      weights = full_weights
  )
  
  
  partial_weights <- c(
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
  
  env_dist_unweighted <- cluster::daisy(env_dt_norm[, envcols_edit, with=F], 
                                        metric = 'gower',
                                        weights = partial_weights
  )
  
  env_dist_weighted_byfecha <- lapply(unique(env_dt_norm$fecha), function(in_fecha) {
    subdt <- env_dt_norm[fecha == in_fecha, -'fecha', with=F]
    env_dist_weighted <- cluster::daisy(subdt, 
                                        metric = 'gower',
                                        weights = full_weights
    )
    return(env_dist_weighted)
  })
  names(env_dist_weighted_byfecha) <- unique(env_dt_norm$fecha)
  
  
  return(list(
    gow_dist_weighted = env_dist_weighted,
    gow_dist_unweighted = env_dist_unweighted,
    gow_dist_weighted_byfecha = env_dist_weighted_byfecha
  ))
}

#--- Compute spatial beta diversity for each date ---------------------------
#in_sp_dt <- tar_read(sp_dt)
compute_spatial_beta <- function(in_sp_dt) {
  #----- Prepare data ----------------------------------------------------------
  #Shroeder and Jenkins for choice of indices
  #Add some based on Anderson et al. 2011
  spcols <- names(
    which(sapply(
      in_sp_dt[, -c('caso', 'total', 'total_relative', 'alpha_div', 'gamma_div')], 
      is.numeric))
  )
  
  # ggplot(data=melt(in_sp_dt[, spcols, with=F])) +
  #   geom_histogram(aes(x=value)) +
  #   scale_x_sqrt()
  
  #----- Jaccard Index for presence-absence by fecha ---------------------------
  presabs_dt <- in_sp_dt[, lapply(.SD, function(x) {as.numeric(x > 0)}),
                         .SDcols= spcols] %>%
    cbind(in_sp_dt[,c('caso', 'fecha', 'presencia'), with=F])
  
  jaccard_mats <- lapply(unique(presabs_dt$fecha), function(t) {
    subdt <- presabs_dt[fecha == t,]
    jaccard_dist <- vegan::vegdist(subdt[fecha == t, spcols, with=F],
                                   method = 'jaccard', diag=T, upper=T)
    jaccard_mat <- as.matrix(jaccard_dist)
    #jaccard_mat[upper.tri(jaccard_mat)] <- NA
    colnames(jaccard_mat) <- subdt$caso
    rownames(jaccard_mat) <- subdt$caso
    
    return(jaccard_mat)
  })
  names(jaccard_mats) <- unique(presabs_dt$fecha)
  
  #----- pres-abs-based total beta diversity, turnover and nestedness by fecha (Baselga) -----
  # beta_jaccard <- lapply(unique(presabs_dt$fecha), function(t) {
  #   subdt <- presabs_dt[fecha == t,]
  #   sub_dist_jac <- betapart.core(subdt[, spcols, with=F]) %>%
  #     beta.multi(index.family="jac")
  #   return(sub_dist_jac)
  # }) %>%
  #   rbindlist %>%
  #   .[, fecha := unique(presabs_dt$fecha)]
  
  #----- pres-abs-based total beta diversity, replacement and richness difference (Legendre) -----
  beta_jaccard <- lapply(unique(presabs_dt$fecha), function(t) {
    subdt <- presabs_dt[fecha == t,]
    sub_dist_jac <- beta.div.comp(subdt[presencia=='si', 
                                        spcols, with=F],
                                  coef = 'J', quant = F) 
    
    return(sub_dist_jac)
  })
  
  beta_jaccard_dt <- lapply(beta_jaccard, 
                            function(out) as.data.table(t(out$part))) %>%
    rbindlist %>%
    .[, fecha := unique(presabs_dt$fecha)]
  
  
  #----- Bray–Curtis dissimilarity distance (untransformed) by fecha -----------
  bray_mats <- lapply(unique(in_sp_dt$fecha), function(t) {
    subdt <- in_sp_dt[fecha == t, ]%>%
      .[, (spcols) := sapply(.SD, simplify=F,
                             function(x) x^(1/5)),
        .SDcols = spcols]
    bray_dist <- vegan::vegdist(subdt[fecha == t, spcols, with=F],
                                method = 'bray', diag=T, upper=T)
    bray_mat <- as.matrix(bray_dist)
    #bray_mat[upper.tri(bray_mat)] <- NA
    colnames(bray_mat) <- subdt$caso
    rownames(bray_mat) <- subdt$caso
    
    return(bray_mat)
  })
  names(bray_mats) <- unique(in_sp_dt$fecha)
  
  # #----- Abundance-based total beta diversity and decomposition (Baselga) -----
  # beta_bray <- lapply(unique(in_sp_dt$fecha), function(t) {
  #   subdt <- in_sp_dt[fecha == t,]
  #   sub_dist_bray <- betapart.core.abund(subdt[, spcols, with=F]) %>%
  #     beta.multi.abund(index.family="bray")
  #   return(sub_dist_bray)
  # }) %>%
  #   rbindlist %>%
  #   .[, fecha := unique(in_sp_dt$fecha)]
  
  #----- Abundance-based total beta diversity and decomposition (Legendre - NOT transformed) -----
  beta_ruzicka <- lapply(unique(in_sp_dt$fecha), function(t) {
    subdt <- in_sp_dt[fecha == t,]
    sub_dist_ruzicka  <- beta.div.comp(subdt[presencia=='si',
                                             (spcols), with=F],
                                       coef = 'J', quant = T
    )
    return(sub_dist_ruzicka)
  })
  
  
  beta_ruzicka_dt <- lapply(beta_ruzicka, 
                            function(out) as.data.table(t(out$part))) %>%
    rbindlist %>%
    .[, fecha := unique(in_sp_dt$fecha)]
  
  #----- Abundance-based total beta diversity and decomposition (Legendre - 5th rt transformed) -----
  beta_ruzicka_trans <- lapply(unique(in_sp_dt$fecha), function(t) {
    subdt <- in_sp_dt[fecha == t,] %>%
      .[, (spcols) := sapply(.SD, simplify=F,
                             function(x) x^(1/5)),
        .SDcols = spcols]
    sub_dist_ruzicka  <- beta.div.comp(subdt[presencia=='si',
                                             (spcols), with=F],
                                       coef = 'J', quant = T
    )
    return(sub_dist_ruzicka)
  })
  
  
  beta_ruzicka_trans_dt <- lapply(beta_ruzicka_trans, 
                                  function(out) as.data.table(t(out$part))) %>%
    rbindlist %>%
    .[, fecha := unique(in_sp_dt$fecha)]
  
  #----- nMDS with Bray-Curtis -------------------------------------------------
  nmds_bray <- in_sp_dt[total !=0, spcols, with=F] %>%
    vegan::vegdist(method = 'bray') %>%
    metaMDS( distance = "bray", trymax=500, maxit=999,
             autotransform=T)
  
  print(nmds_bray)
  
  nmds_bray_dt <- cbind(in_sp_dt[total !=0, -spcols, with=F],
                        nmds_bray$points
  )
  
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
  #                      'estado_de_flujo', 'total', 'alpha_div'),
  #                  with = F],
  #         by.x = 'caso.x', by.y='caso') %>%
  #   merge(in_sp_dt[, c('caso', 'cuenca', 'sitio', 'fecha', 'intermitencia',
  #                      'estado_de_flujo', 'total', 'alpha_div'),
  #                  with = F],
  #         by.x = 'caso.y', by.y='caso')
  
  #--- Function outputs ------------
  return(list(
    presabs_dt = presabs_dt,
    jaccard_mats = jaccard_mats,
    beta_jaccard = beta_jaccard,
    beta_jaccard_dt = beta_jaccard_dt,
    bray_mats = bray_mats,
    #beta_bray = beta_bray,
    beta_ruzicka = beta_ruzicka,
    beta_ruzicka_dt = beta_ruzicka_dt,
    beta_ruzicka_trans = beta_ruzicka_trans,
    beta_ruzicka_trans_dt = beta_ruzicka_trans_dt,
    nmds_bray = nmds_bray,
    nmds_bray_dt =  nmds_bray_dt,
    sp_chord_dist = sp_chord_dist,
    chord_abundance_pca = chord_abundance_pca,
    pca_chord_dt = pca_chord_dt
  ))
} 
#--- Plot spatial beta div --------------------------------------------------
#in_spatial_beta <- tar_read(spatial_beta)

plot_spatial_beta <- function(in_spatial_beta) {
  
  #Plot beta diversity and its components over time based on pres-abs Jaccard
  beta_jaccard_melt <- in_spatial_beta$beta_jaccard_dt %>%
    melt(id.vars='fecha') %>%
    .[variable %in% c('BDtotal', 'Repl', 'RichDif'),]
  
  plot_beta_jaccard_fecha <- ggplot() +
    geom_area(data = beta_jaccard_melt[variable != 'BDtotal'],
              aes(x=fecha, y=value, 
                  position = 'stack',
                  fill=variable, group=variable)) +
    geom_line(data = beta_jaccard_melt[variable == 'BDtotal'],
              aes(x=fecha, y=value, group=variable, color=variable),
              size=1) + 
    scale_color_manual(values ='black',
                       labels = 'Total beta diversity') +
    scale_fill_discrete(labels=c('Replacement', 'Richness difference')) +
    scale_x_discrete(name='Date') +
    scale_y_continuous(limits=c(0,1), 
                       name='Value') +
    coord_cartesian(expand=c(0,0)) +
    theme_classic() + 
    theme(legend.title = element_blank())
  
  # #Plot beta diversity and its components over time based on abundance Bray-Curtis
  # beta_bray_melt <- in_spatial_beta$beta_bray %>%
  #   melt(id.vars='fecha')
  # 
  # plot_beta_bray_fecha <- ggplot() +
  #   geom_area(data = beta_bray_melt[variable != 'BDtotal'],
  #             aes(x=fecha, y=100*value, 
  #                 position = 'stack',
  #                 fill=variable, group=variable)) +
  #   geom_line(data = beta_bray_melt[variable == 'BDtotal'],
  #             aes(x=fecha, y=100*value, group=variable, color=variable),
  #             size=1) + 
  #   scale_color_manual(values ='black',
  #                      labels = 'Total beta diversity') +
  #   scale_fill_discrete(labels=c('Replacement',
  #                                'Richness difference')) +
  #   scale_x_discrete(name='Date') +
  #   scale_y_continuous(limits=c(0,100), 
  #                      name='Value') +
  #   coord_cartesian(expand=c(0,0)) +
  #   theme_classic() + 
  #   theme(legend.title = element_blank())
  
  #Plot beta diversity and its components over time based on abundance Ruzicka
  beta_ruzicka_melt <- in_spatial_beta$beta_ruzicka_dt %>%
    melt(id.vars='fecha') %>%
    .[variable %in% c('BDtotal', 'Repl', 'RichDif'),]
  
  plot_beta_ruzicka_fecha <- ggplot() +
    geom_area(data = beta_ruzicka_melt[variable != 'BDtotal'],
              aes(x=fecha, y=value,
                  position = 'stack',
                  fill=variable, group=variable)) +
    geom_line(data = beta_ruzicka_melt[variable == 'BDtotal'],
              aes(x=fecha, y=value, group=variable, color=variable),
              size=1) +
    scale_color_manual(values ='black',
                       labels = 'Total beta diversity') +
    scale_fill_discrete(labels=c('Replacement',
                                 'Richness difference')) +
    scale_x_discrete(name='Date') +
    scale_y_continuous(limits=c(0,1),
                       name='Value') +
    coord_cartesian(expand=c(0,0)) +
    theme_classic() +
    theme(legend.title = element_blank())
  
  #Plot beta diversity and its components over time based on abundance Ruzicka (after transforming data)
  beta_ruzicka_trans_melt <- in_spatial_beta$beta_ruzicka_trans_dt %>%
    melt(id.vars='fecha') %>%
    .[variable %in% c('BDtotal', 'Repl', 'RichDif'),]
  
  plot_beta_ruzicka_trans_fecha <- ggplot() +
    geom_area(data = beta_ruzicka_trans_melt[variable != 'BDtotal'],
              aes(x=fecha, y=value,
                  position = 'stack',
                  fill=variable, group=variable)) +
    geom_line(data = beta_ruzicka_trans_melt[variable == 'BDtotal'],
              aes(x=fecha, y=value, group=variable, color=variable),
              size=1) +
    scale_color_manual(values ='black',
                       labels = 'Total beta diversity') +
    scale_fill_discrete(labels=c('Replacement',
                                 'Richness difference')) +
    scale_x_discrete(name='Date') +
    scale_y_continuous(limits=c(0,1),
                       name='Value') +
    coord_cartesian(expand=c(0,0)) +
    theme_classic() +
    theme(legend.title = element_blank())
  
  #Plot trajectories
  #nMDS for Bray-Curtis
  # ggplot(data=in_spatial_beta$nmds_bray_dt[intermitencia=='isolpool',],
  #        aes(x=MDS1, y=MDS2, color=estado_de_flujo)) +
  #   geom_point() +
  #   geom_line(aes(group=sitio))
  
  nMDS_dt_forplot <- in_spatial_beta$nmds_bray_dt %>%
    merge(
      setnames(
        expand.grid(c(paste0('fecha', seq(1,6))),
                    unique(.$sitio)),
        c('fecha', 'sitio')
      ),
      by=c('fecha', 'sitio'),
      all.y=T
    ) %>%
    .[, `:=`(MDS1_lag = lag(MDS1),
             MDS2_lag = lag(MDS2),
             estado_de_flujo_lag = lag(estado_de_flujo)),
      by=sitio]
  
  nMDS_dt_forplot
  
  plot_nmds_time <- ggplot() +
    geom_point(data= nMDS_dt_forplot,
               aes(x = MDS1_lag, y = MDS2_lag, 
                   shape=estado_de_flujo_lag),
               color = 'grey',
               size=4, alpha=0.5) +
    geom_point(data = nMDS_dt_forplot[is.na(MDS1),],
               aes(x = MDS1_lag, y = MDS2_lag, 
                   shape=estado_de_flujo_lag),
               color = 'black',
               size=4, alpha=0.4) +
    geom_segment(data= nMDS_dt_forplot,
                 aes(x=MDS1_lag, xend=MDS1,
                     y=MDS2_lag, yend=MDS2),
                 color = 'grey', alpha=0.5) +
    geom_point(data= nMDS_dt_forplot,
               aes(x = MDS1, y = MDS2, 
                   shape=estado_de_flujo, color=intermitencia),
               size=4, alpha=0.85) +
    scale_shape_discrete(name='Flow state at time step',
                         labels = c('Flowing',
                                    'Disconnected pools',
                                    '')) +
    scale_color_manual(name='Long-term flow regime',
                       values = c('#d73027', '#fee090', '#4575b4', '#ffffff'),
                       labels = c('Intermittent: dry',
                                  'Intermittent: disconnected pools',
                                  'Perennial',
                                  '')) +
    scale_x_continuous(name = 'MDS1') +
    scale_y_continuous(name = 'MDS2') + 
    coord_fixed() +
    facet_wrap(~fecha) +
    theme_classic()
  
  return(list(
    plot_beta_jaccard_fecha = plot_beta_jaccard_fecha,
    plot_beta_ruzicka_fecha = plot_beta_ruzicka_fecha,
    plot_beta_ruzicka_trans_fecha = plot_beta_ruzicka_trans_fecha,
    plot_nmds_time = plot_nmds_time
  ))
}


#--- Compute mantel test -----------------------------------------------------
# in_spenv_dt = tar_read(spenv_dt)
# in_spatial_beta = tar_read(spatial_beta)
# in_env_dist = tar_read(env_dist)
# in_euc_dist = tar_read(euc_dist)
# in_net_dist = tar_read(net_dist)
# rep = 999

compute_mantel <- function(
    in_spenv_dt,
    in_spatial_beta,
    in_euc_dist,
    in_env_dist,
    in_net_dist,
    rep = 999) {
  #We used traditional Mantel tests for spatial distances 
  #(topographic and network distance) and Mantel tests corrected by spatial 
  #autocorrelation through Moran spectral randomization (MSR; Crabot et al. 2019) 
  #for environmental distances (water chemistry and flow regime difference)
  #(999 runs for each test) to identify the relative effect of environmental 
  #and spatial filters on community composition in each year. The MSR correction 
  #was able to remove the spurious spatial dependence from our environmental 
  #distances, providing a correlation value that reflected the net environmental 
  #importance. Each Mantel test included the Bray–Curtis dissimilarity based on
  #macroinvertebrate abundances as response variable and one environmental 
  #(water chemistry, hydrological) or spatial (topographic, network) distance 
  #as predictor. [from Arguelles et al. 2020]
  
  #------------------ Prepare data ---------------------------------------------
  #Get latitute/northing and longitude/easting
  xy <- as.matrix(in_euc_dist$xy)
  
  #Extract pres-abs data
  presabs_dt_spcols_fecha <- in_spatial_beta$presabs_dt[
    , -c('caso', 'presencia')]
  
  spcols <-  names(presabs_dt_spcols_fecha[, -'fecha', with=F])
  
  #Extract abundance data
  spcols_abund <- in_spenv_dt[, names(presabs_dt_spcols_fecha),
                              with=F]
  spcols_abund_trans <- spcols_abund %>%
    .[, (spcols) := sapply(.SD, simplify=F,
                           function(x) x^(1/5)),
      .SDcols = spcols]
  
  #Check multivariate empirical variogram for Jaccard distance
  vario_presabs <- adespatial::variogmultiv(presabs_dt_spcols_fecha[fecha=='fecha1', 
                                                                    -'fecha', with=F],
                                            xy,
                                            nclass=5)
  
  plot(
    vario_presabs$d,
    vario_presabs$var,
    ty='b',
    pch=20,
    xlab='Distance',
    ylab=("C(distance")
  )
  
  #Check multivariate empirical variogram for Bray-Curtis distance
  vario_abund <- spcols_abund[fecha=='fecha1', 
                              -'fecha', with=F] %>%
    .[, sapply(.SD, function(x) x^1/5)] %>%
    variogmultiv(xy,
                 nclass=5)
  
  plot(
    vario_abund$d,
    vario_abund$var,
    ty='b',
    pch=20,
    xlab='Distance',
    ylab=("C(distance")
  )
  
  #------------------ Test Crabot et al. 2019 procedure -----------------------
  run_msr_mantel <- function(in_xy,
                             in_resp = NULL,
                             in_pred_dist,
                             in_resp_dist,
                             sqrt_pred_dist,
                             simple_test,
                             rep = 999,
                             nullify_NAs = TRUE,
                             verbose = F
  ) {
    
    if (nullify_NAs) {in_resp_dist[is.na(in_resp_dist)] <- 0} #Double 0s
    
    if (sqrt_pred_dist) {in_pred_dist <- sqrt(in_pred_dist)}
    
    if (simple_test) {
      mantel_out<- ade4::mantel.randtest(
        as.dist(in_pred_dist), 
        as.dist(in_resp_dist)
      )
      
      r <-  mantel_out$obs - mantel_out$expvar["Expectation"]
      names(r) <- ""
      signif <- mantel_out$pvalue
      test <- 'ade4::mantel.randtest - simple'
      
    } else { #Run partial mantel test
      #If the pred distance matrix is non-euclidean
      if (!(ade4::is.euclid(as.dist(in_pred_dist)))) {
        #adespatial::listw.explore() #Too explore different SWM models
        #Perform a standard partial mantel test
        mantel_out <- vegan::mantel.partial(xdis = as.dist(in_pred_dist), 
                                            ydis = as.dist(in_resp_dist),
                                            zdis = dist(in_xy),
                                            permutations = rep)
        r <- mantel_out$statistic
        signif <- mantel_out$signif
        test <- 'vegan::mantel.partial - standard partial'
        
        #If the pred distance matrix is euclidean
      } else {
        
        # Optimize spatial weight matrix before Moran Spectral Randomization
        sw_candidates <- adespatial::listw.candidates(in_xy, style='C') 
        
        sw_selected <- adespatial::listw.select(
          x = in_resp,
          candidates = sw_candidates,
          MEM.autocor = 'positive',
          method = 'FWD',
          p.adjust = TRUE,
          verbose = verbose
        )
        
        #If there was no significant positive spatial structure
        if (is.null(sw_selected$best.id)) {
          
          #Perform a standard partial mantel test
          mantel_out <- vegan::mantel.partial(xdis = as.dist(in_pred_dist), 
                                              ydis = as.dist(in_resp_dist),
                                              zdis = dist(in_xy),
                                              permutations = rep)
          r <- mantel_out$statistic
          signif <- mantel_out$signif
          test <- 'vegan::mantel.partial - standard partial'
          
        } else {
          #If a significant positive spatial structure was detected, perform 
          #Mantel test based on spatially constrained randomizations using
          #Moran spectral randomization with the best spatial weighting matrix
          R2_best <- round(sw_selected$candidates$R2Adj.select[sw_selected$best.id], 3)
          #Compute SW based on best model
          lw1 <- sw_candidates[[names(sw_selected$best.id)]]
          
          #Run standar Mantel
          standard_mantel <- ade4::mantel.randtest(m1=as.dist(in_pred_dist), 
                                                   m2=as.dist(in_resp_dist)
          )
          
          #Then 
          mantel_out <- adespatial::msr(standard_mantel, 
                                        lw1, 
                                        nrepet=rep)
          
          r <-  mantel_out$obs - mantel_out$expvar["Expectation"]
          names(r) <- ""
          signif <- mantel_out$pvalue
          test <- 'adespatial::msr - msr partial'
        }
      }
    }
    
    return(list(
      mantel_out = mantel_out,
      r = r,
      signif = signif,
      test = test)
    )
  }
  
  #----------------- Compute partial mantel test -------------------------------
  #Jaccard distance (presence-absence data) - Environmental distance
  print("Jaccard distance (presence-absence data) - Environmental distance")
  jaccard_mantel_envdist_list <- lapply(
    unique(presabs_dt_spcols_fecha$fecha),
    function(in_fecha) { 
      print(in_fecha)
      return(
        run_msr_mantel(
          in_xy = xy,
          in_resp = presabs_dt_spcols_fecha[fecha==in_fecha, -'fecha', with=F],
          in_pred_dist = as.matrix(in_env_dist$gow_dist_weighted_byfecha[[in_fecha]]),
          in_resp_dist =  in_spatial_beta$jaccard_mats[[in_fecha]],
          simple_test = FALSE,
          sqrt_pred_dist = TRUE,
          nullify_NAs = TRUE,
          rep = 999
        )
      )}
  )
  
  #Jaccard distance (presence-absence data) - Network distance
  print("Jaccard distance (presence-absence data) - Network distance")
  jaccard_mantel_netdist_list <- lapply(
    unique(presabs_dt_spcols_fecha$fecha),
    function(in_fecha) { 
      print(in_fecha)
      return(
        run_msr_mantel(
          in_xy = xy,
          in_resp = presabs_dt_spcols_fecha[fecha==in_fecha, -'fecha', with=F],
          in_pred_dist = as.matrix(in_net_dist),
          in_resp_dist =  in_spatial_beta$jaccard_mats[[in_fecha]],
          simple_test = FALSE,
          sqrt_pred_dist = FALSE,
          nullify_NAs = TRUE,
          rep = 999
        )
      )}
  )
  
  #Bray-Curtis distance (abundance data) - Environmental distance
  print("Bray-Curtis distance (abundance data) - Environmental distance")
  bray_mantel_envdist_list <- lapply(
    unique(spcols_abund_trans$fecha),
    function(in_fecha) { 
      print(in_fecha)
      return(
        run_msr_mantel(
          in_xy = xy,
          in_resp = spcols_abund_trans[fecha=='fecha1', -'fecha', with=F],
          in_pred_dist = as.matrix(in_env_dist$gow_dist_weighted_byfecha[[in_fecha]]),
          in_resp_dist =  in_spatial_beta$bray_mats[[in_fecha]],
          simple_test = FALSE,
          sqrt_pred_dist = TRUE,
          nullify_NAs = TRUE,
          rep = 999
        )
      )}
  )
  
  #Bray-Curtis distance (abundance data) - Network distance
  print("Bray-Curtis distance (abundance data) - Network distance")
  bray_mantel_netdist_list <- lapply(
    unique(spcols_abund_trans$fecha),
    function(in_fecha) { 
      print(in_fecha)
      return(
        run_msr_mantel(
          in_xy = xy,
          in_resp = spcols_abund_trans[fecha=='fecha1', -'fecha', with=F],
          in_pred_dist = as.matrix(in_net_dist),
          in_resp_dist =  in_spatial_beta$bray_mats[[in_fecha]],
          simple_test = FALSE,
          sqrt_pred_dist = FALSE,
          nullify_NAs = TRUE,
          rep = 999
        )
      )}
  )
  
  #------------------ Run simple Mantel test against euclidean distance --------
  #Jaccard distance (presence-absence data) - Euclidean distance
  print("Jaccard distance (presence-absence data) - Euclidean distance")
  jaccard_mantel_eucdist_list <- lapply(
    unique(presabs_dt_spcols_fecha$fecha),
    function(in_fecha) { 
      print(in_fecha)
      
      out_list <- run_msr_mantel(
        in_xy = xy,
        in_pred_dist = as.matrix(in_euc_dist$distmat),
        in_resp_dist =  in_spatial_beta$jaccard_mats[[in_fecha]],
        simple_test = TRUE,
        sqrt_pred_dist = FALSE,
        nullify_NAs = TRUE,
        rep = 999
      )
      return(out_list)
    }
  )
  
  #Bray-Curtis distance (abundance data) - Euclidean distance
  print("Bray-Curtis distance (abundance data) - Euclidean distance")
  bray_mantel_eucdist_list <- lapply(
    unique(presabs_dt_spcols_fecha$fecha),
    function(in_fecha) { 
      print(in_fecha)
      
      out_list <- run_msr_mantel(
        in_xy = xy,
        in_pred_dist = as.matrix(in_euc_dist$distmat),
        in_resp_dist =  in_spatial_beta$bray_mats[[in_fecha]],
        simple_test = TRUE,
        sqrt_pred_dist = FALSE,
        nullify_NAs = TRUE,
        rep = 999
      )
      return(out_list)
    }
  )
  
  mantel_netdist_eucdist <- mantel_out<- ade4::mantel.randtest(
    as.dist(in_euc_dist$distmat), 
    as.dist(in_net_dist)
  )
  
  return(list(
    jaccard_mantel_envdist_list = jaccard_mantel_envdist_list,
    jaccard_mantel_netdist_list = jaccard_mantel_netdist_list,
    jaccard_mantel_eucdist_list = jaccard_mantel_eucdist_list,
    bray_mantel_envdist_list = bray_mantel_envdist_list,
    bray_mantel_netdist_list = bray_mantel_netdist_list,
    bray_mantel_eucdist_list = bray_mantel_eucdist_list,
    mantel_netdist_eucdist = mantel_netdist_eucdist
  ))
}

#--- Plot Mantel tests results ----------------------------------------------
#in_mantel_test_list <- tar_read(mantel_test_list)

plot_mantel_tests <- function(in_mantel_test_list) {
  
  #Format all mantel test results in a single dt
  mantel_dt <- 
    lapply(list('jaccard', 'bray'), function(dissimilarity) {
      lapply(
        list('envdist', 'netdist', 'eucdist'), function(mat) {
          lapply(
            in_mantel_test_list[[paste0(
              dissimilarity, '_mantel_', mat, '_list')]], 
            function(x) x[2:4]) %>%
            rbindlist %>%
            .[, `:=`(variable = mat,
                     dissim = dissimilarity,
                     fecha = paste0('fecha', seq(1,6)))]
        }) %>%
        rbindlist
    }) %>%
    rbindlist
  
  mantel_plot <- ggplot(mantel_dt, aes(x=fecha, y=r)) +
    geom_path(aes(color=variable, linetype = dissim,
                  group = interaction(dissim, variable))) +
    geom_point(aes(color=variable, shape = dissim), size=3) +
    geom_point(data = mantel_dt[signif > 0.05,], 
               aes(shape = dissim), color='grey', size=3) +
    scale_color_discrete(name = '', 
                         labels=c('Environment', 'Euclidean distance',
                                  'Network distance')) +
    scale_shape_discrete(name='Dissimilarity measure',
                         labels=c('Bray-Curtis (abundance)', 
                                  'Jaccard (presence-absence)')) +
    scale_linetype_discrete(name='Dissimilarity measure',
                            labels=c('Bray-Curtis (abundance)', 
                                     'Jaccard (presence-absence)')) +
    theme_classic() +
    theme(axis.title.x = element_blank())
  
  return(list(
    mantel_plot = mantel_plot,
    mantel_dt = mantel_dt
  ))
  
  
  
  
  
  
  
  
}

#--- Download basemap data -------
#out_path <- file.path(resdir, 'basemap_data')
#in_net <- tar_read(net_formatted)

download_basemap <- function(out_path
                             # , in_net
) {
  if (!dir.exists(out_path)) {
    dir.create(out_path)
  }
  
  #------- Download administrative boundaries----------------------------------
  admin <- geodata::world(resolution=3, path=out_path)
  
  bolivia_boundaries <- admin[admin$NAME_0 == 'Bolivia']
  
  #------- Download low-res DEM for all of Bolivia -----------------------------
  elev_bolivia <- geodata::elevation_30s(country='Bolivia',
                                         path=file.path(out_path, 'strm'))
  
  
  #------- Download high-res DEM for watershed ---------------------------------
  #net_bbox <- ext(terra::project(x=vect(in_net), crs(elv_bolivia)))
  # elev_tiles <- geodata::elevation_3s(lon=net_bbox[1:2], lat=net_bbox[3:4],
  #                                     path = file.path(out_path, 'strm'))
  
  tile <- '23_16'
  dem_path = file.path(out_path, paste0('srtm_', tile, '.zip'))
  if (!file.exists(dem_path)) {
    download.file(
      paste0('https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_',
             tile, '.zip'), 
      dem_path)
  }
  
  elev_net <- unzip(dem_path, exdir = out_path) %>%
    grep('.tif', ., value = T) %>%
    terra::rast(.)
  
  # # tile_id_list <- c('23_14', '23_15', '23_16', '23_17',
  # #          '24_14', '24_15', '24_16', '24_17',
  # #          '25_16', '25_15', '25_16', '25_17')
  # 
  # elev_tiles <- lapply(tile_id_list, function(tile) {
  #   print(tile)
  #   dem_path = file.path(out_path, paste0('srtm_', tile, '.zip'))
  #   if (!file.exists(dem_path)) {
  #     download.file(
  #       paste0('https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_',
  #              tile, '.zip'), 
  #       dem_path)
  #   }
  #   
  #   elev <- unzip(dem_path, exdir = out_path) %>%
  #     grep('.tif', ., value = T) %>%
  #     terra::rast(.)
  #   
  #   return(elev)
  # })
  # 
  # #------- Mosaick and crop DEM  -----------------------------------------------
  # #dem_out_path <- file.path(out_path, 'dem_bolivia')
  # 
  # elev_crop <- sprc(elev_tiles) %>%
  #   merge %>%  
  #   crop(ext(bolivia_boundaries) + 0.1) %>%
  #   mask(bolivia_boundaries)
  # 
  # names(elev_crop) <- 'elevation'
  # 
  # out_elev_rast <- file.path(out_path, 'elev_crop.tif')
  # writeRaster(elev_crop, out_elev_rast, overwrite = T)
  # 
  # #lc <- geodata::landcover('trees')
  
  #------- Function outputs ----------------------------------------------------
  return(list(
    admin = terra::serialize(admin, NULL),
    elev_bolivia = terra::serialize(elev_bolivia, NULL),
    elev_net = terra::serialize(elev_net, NULL)
  ))
}


#--- Create hillshade ---------------------------------------------------------
create_hillshade <- function(in_dem, z_exponent, write=F, out_path) {
  #From Dr. Dominic Royé - https://dominicroye.github.io/en/2022/hillshade-effects/
  #terra::rast(in_basemaps$elev_path) %>%
  elev <- in_dem %>%
    terra::project("epsg:32720") %>%
    .^(z_exponent)
  
  # estimate the slope
  sl <- terra::terrain(elev, "slope", unit = "radians")
  
  # estimate the aspect or orientation
  asp <- terra::terrain(elev, "aspect", unit = "radians")
  
  # pass multiple directions to shade()
  hillmulti <- purrr::map(c(270, 15, 60, 330), function(dir){ 
    shade(sl, asp, 
          angle = 45, 
          direction = dir,
          normalize= TRUE)}
  ) %>%
    rast %>%
    sum
  
  if (write) {
    terra::writeRaster(hillmulti, out_path, overwrite = T)
    return(
      out_rast
    )
  } else {
    return(serialize(hillmulti, NULL))
  }
  
}


#--- Map sites --------------------------------------------------------------
in_spenv_dt = tar_read(spenv_dt)
in_net = tar_read(net_directed)
in_sites_path = tar_read(sites_path)
in_basemaps <- tar_read(basemaps)
in_hillshade_bolivia <- tar_read(hillshade_bolivia)
in_hillshade_net <- tar_read(hillshade_net)

map_caynaca <- function(in_spenv_dt, 
                        in_net,
                        in_sites_path,
                        in_basemaps,
                        in_hillshade_bolivia,
                        in_hillshade_net) {
  #------------ Make watershed map ---------------------------------------------
  netbbox <- ext(vect(in_net))
  
  elev_net <- unserialize(in_basemaps$elev_net) %>%
    terra::project("epsg:32720") %>%
    crop(netbbox + 1000) 
  
  #Load and crop hillshade
  hilldf_net <- unserialize(in_hillshade_net) %>%
    crop(netbbox + 1000) 
  
  #Format Hillshade map - https://dieghernan.github.io/202210_tidyterra-hillshade/
  # normalize names
  names(hilldf_net) <- "shades"
  # Make palette
  pal_greys <- hcl.colors(1000, "Grays")
  # Use a vector of colors
  index <- hilldf_net %>%
    tidyterra::mutate(index_col = rescale(shades, to = c(1, length(pal_greys)))) %>%
    tidyterra::mutate(index_col = round(index_col)) %>%
    tidyterra::pull(index_col)
  # Get cols
  vector_cols <- pal_greys[index]
  
  axis_ext <- vect(in_net) %>%
    project("EPSG:4326") %>%
    ext() %>%
    as.vector()
  
  br_y <- seq(axis_ext[3], axis_ext[4], length.out = 1000) %>%
    pretty(n = 3) %>%
    round(3) %>%
    unique()
  
  br_x <- seq(axis_ext[1], axis_ext[2], length.out = 1000) %>%
    pretty(n = 3) %>%
    round(3) %>%
    unique()
  
  hill_plot <- ggplot() +
    geom_spatraster(
      data = hilldf_net, fill = vector_cols, maxcell = Inf,
      alpha = 0.6
    ) +
    scale_x_continuous(breaks=br_x) +
    scale_y_continuous(breaks=br_y)
  
  #Format full map
  r_limits <- minmax(elev_net$srtm_23_16) %>% as.vector() 
  r_limits <-  c(floor(r_limits[1] / 500), 
                 ceiling(r_limits[2] / 500)) * 500 %>%
    pmax(0)
  
  #Test palette
  # elevt_test <- ggplot() +
  #   geom_spatraster(data = elev_net, aes(fill=srtm_23_16))
  # plot_pal_test <- function(pal) {
  #   elevt_test +
  #     scale_fill_hypso_tint_c(
  #       limits = r_limits,
  #       palette = pal
  #     ) +
  #     ggtitle(pal) +
  #     theme_minimal()
  # }
  # 
  # plot_pal_test("etopo1_hypso")
  # plot_pal_test("dem_poster")
  # plot_pal_test("spain")
  # plot_pal_test("pakistan")
  # plot_pal_test("utah_1")
  # plot_pal_test("wiki-2.0_hypso")
  
  base_plot <- hill_plot +
    # Avoid resampling with maxcell
    geom_spatraster(data =  elev_net, maxcell = Inf) +
    scale_fill_hypso_tint_c(
      limits = r_limits,
      palette = "etopo1_hypso",
      alpha = 0.3,
      labels = label_comma(),
      # For the legend I use custom breaks
      breaks = c(
        seq(0, 500, 100),
        seq(750, 1500, 250),
        2000
      )
    )
  
  net_plot <- base_plot +
    geom_sf(data = in_net, aes(linewidth=stream_order), color='blue') +
    scale_linewidth_continuous(range=c(0.5,1.7)) +
    geom_spatvector(data = vect(in_sites_path), 
                    color='black',
                    size = 1.75) +
    #theme_minimal(base_family = "notoserif") +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          panel.background = element_rect(color='white'))
  
  
  #------------ Make inset map -------------------------------------------------
  
  
  admin <- unserialize(in_basemaps$admin)%>%
    terra::project("epsg:32720")
  elev <- rast(in_basemaps$elev_path) %>%
    terra::project("epsg:32720")
  
  netbbox <- st_bbox(in_net)
  
  #Convert the hillshade to XYZ
  hilldf <- rast(in_hillshade) %>%
    as.data.frame(xy = TRUE)
  
  ggplot() +
    geom_raster(data = hilldf,
                aes(x, y, fill = sum),
                show.legend = FALSE) +
    scale_fill_distiller(palette = "Greys") +
    geom_spatraster(data = elev,
                    aes(fill=BOL_elv_msk),
                    alpha=0.2) +
    geom_sf(data = in_net) +
    # scale_fill_hypso_tint_c(breaks = c(180, 250, 500, 1000,
    #                                    1500,  2000, 2500,
    #                                    3000, 3500, 4000)) +
    # guides(fill = guide_colorsteps(barwidth = 20,
    #                                barheight = .5,
    #                                title.position = "right")) +
    labs(fill = "m") +
    xlim(c(netbbox[1]-5000, netbbox[3]+5000)) +
    ylim(c(netbbox[2]-5000, netbbox[4]+5000)) +
    theme_void() +
    theme(legend.position = "bottom")
  
}
