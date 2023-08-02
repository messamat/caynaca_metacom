source("src/caynaca_packages.R")
source("src/caynaca_functions.R")

rootdir <- here::here()
datadir <- file.path(rootdir, 'data')
resdir <- file.path(rootdir, 'results')

sites_snapped_out <- file.path(resdir, 'sites_snapped.shp')
sites_IDcol <- 'C___DIGO_N'
overwrite <- T

#Inspect data
list(
  #Download h90 basins
  tar_target(
    rawdata_path,
    file.path(datadir, 'alltaxa.xls')
  )
  ,
  
  tar_target(
    net_shp_path,
    file.path(datadir, 'MatrixNaborA', 'FNetRiver.shp')
  )
  ,
  
  tar_target(
    sites_shp_path,
    file.path(datadir, 'MatrixNaborA', 'FSites.shp')
  )
  ,
  
  tar_target(
    sp_rawdt,
    readxl::read_xls(rawdata_path, sheet="alltaxa") %>%
      setDT
  ),
  
  tar_target(
    env_rawdt,
    readxl::read_xls(rawdata_path, sheet="environ") %>%
      setDT
  )
  ,
  
  tar_target(
    sp_dt,
    format_spdata(in_sp_dt = sp_rawdt)
  )
  ,
  
  tar_target(
    env_dt,
    format_envdata(in_env_dt = env_rawdt)
  )
  ,
  
  tar_target(
    spenv_dt,
    merge(sp_dt,
          env_dt[, c('caso', 
                     names(.SD)[!(names(.SD) %in% names(sp_dt))]),
                 with=F], 
          by='caso')
  )
  ,
  
  tar_target(
    sp_plots,
    plot_spdata(in_spenv_dt = spenv_dt, 
                in_sp_dt = sp_dt)
  ),
  
  tar_target(
    env_plots,
    plot_envdata(in_env_dt = env_dt)
  )
  ,
  
  tar_target(
    net_formatted,
    fix_net_dangles(in_net = net_shp_path,
                    in_tolerance = 7) %>%
      dplyr::mutate(length = as.numeric(st_length(.))) %>%
      .[st_is_valid(.),]
  ),
  
  #Snap sites
  tar_target(
    sites_snapped_path,
    snap_sites(in_sites_point_path = sites_shp_path,
               in_sitesSQL="",
               in_segments = net_formatted,
               out_path = sites_snapped_out,
               proj = 'original',
               pt_IDcol = sites_IDcol,
               nearby_ID = T,
               seg_IDcol = 'OBJECTID_1',
               overwrite = overwrite)
  )
  ,
  
  #Compute euclidean geographic distance
  tar_target(
    euc_dist,
    compute_eucdist(in_sites_path=sites_snapped_path,
                    IDcol = sites_IDcol)
  )
  ,
  
  #Compute network distance
  tar_target(
    net_dist,
    compute_netdist(in_sites_snapped = sites_snapped_path, 
                    in_net = net_formatted,
                    IDcol = sites_IDcol)
  )
  ,

  #Compute environmental distance matrix
  tar_target(
    env_dist,
    compute_envdist(
      in_env_dt = env_dt,
      IDcol = sites_IDcol
    )
  )

)








#LAter on: Compute topographic distance (terrain is accidented)







#traditional Mantel tests for spatial distances (topographic and network distance)
#Bray–Curtis dissimilarity based on macroinvertebrate abundances as response variable

#Mantel tests corrected by spatial autocorrelation through Moran spectral randomization 
#(MSR; Crabot, Clappe, Dray, & Datry, 2019) for environmental distances (999 runs for each test) 
#to identify the relative effect of environmental and spatial filters on community composition in each year.
#Bray–Curtis dissimilarity based on macroinvertebrate abundances as response variable

#Before summing Mantel r-values, we tested the collinearity between distances (Canedo-Arguelles et al. 2020).
#Each pair of distances was only weakly correlated (Pearson r: −0.20 to 0.25),

#More advanced methods

#Make a map of the sites, river network, watershed - with site color showing flow permanence

#Make a plot of Mantel's R over time


#https://public.igb-berlin.de/index.php/s/agciopgzXjWswF4/download?path=%2Fr.stream.order%2Forder_vect_tiles20d&files=order_vect_segment_h10v10.gpkg