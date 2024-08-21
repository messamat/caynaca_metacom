suppressMessages(source("src/caynaca_packages.R"))
source("src/caynaca_functions.R")

rootdir <- here::here()
datadir <- file.path(rootdir, 'data')
resdir <- file.path(rootdir, 'results')

sites_snapped_out <- file.path(resdir, 'sites_snapped.shp')
sites_IDcol <- 'C___DIGO_N'
overwrite <- T

#Inspect data
list(
  tar_target(
    rawdata_path,
    file.path(datadir, 'alltaxa_mod_mar2024_edit.xls'),
    format = "file"
  )
  ,
  
  tar_target(
    traitsdata_path,
    file.path(datadir, 'bio_traits.csv'),
    format = "file"
  )
  ,
  
  tar_target(
   familydata_path,
   file.path(datadir, 'relative_abundance_bio_traits.xlsx'),
   format = "file"
  )
  ,
  
  tar_target(
    net_path,
    file.path(datadir, 'MatrixNaborA', 'FNetRiver.shp'),
    format = "file"
  )
  ,
  
  tar_target(
    sites_path,
    file.path(datadir, 'MatrixNaborA', 'FSites.shp'),
    format = "file"
  )
  ,
  
  tar_target(
    sp_rawdt,
    readxl::read_xls(rawdata_path, sheet="alltaxa") %>%
      setDT
  ),
  
  tar_target(
    traits_rawdt,
    fread(traitsdata_path) %>%
      setDT
  )
  ,
  
  tar_target(
    family_rawdt,
    readxl::read_xlsx(familydata_path, sheet="trsposed") %>%
      setDT %>%
      setnames("Policentropodidae", "Polycentropodidae") %>%
      setnames(function(x) {str_trim(tolower(x))})
  )
  ,
  
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
                     names(env_dt)[!(names(env_dt) %in% names(sp_dt))]),
                 with=F], 
          by='caso') %>%
      .[, estado_de_flujo := factor(
        estado_de_flujo, 
        levels=c('F', 'IP', 'D'),
        labels = c('Flowing', 'Isolated pools', 'Dry bed'))
      ]
  )
  ,
  
  tar_target(
   traits_dt,
   format_traitsdata(in_traits_rawdt = traits_rawdt)
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
  
  # tar_target(
  #   func_plots,
  #   plot_funcdata(
  #     
  #   )
  # )
  
  tar_target(
    net_formatted,
    fix_net_dangles(in_net = net_path,
                    in_tolerance = 7) %>%
      dplyr::mutate(length = as.numeric(st_length(.))) %>%
      .[st_is_valid(.),]
  ),
  
  #Snap sites
  tar_target(
    sites_snapped_path,
    snap_sites(in_sites_point_path = sites_path,
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
  ),
  
  tar_target(
    spatial_taxo_beta,
    compute_spatial_taxo_beta(
      in_sp_dt = sp_dt
    )
  )
  ,
  
  tar_target(
    spatial_beta_plots,
    plot_spatial_beta(
      in_spatial_beta = spatial_taxo_beta
    )
  )
  ,

  tar_target(
    mantel_dt,
    compute_mantel(in_spenv_dt = spenv_dt,
                   in_spatial_beta = spatial_taxo_beta,
                   in_env_dist = env_dist,
                   in_euc_dist = euc_dist,
                   in_net_dist = net_dist)
  ),
  
  tar_target(
    mantel_plots,
    plot_mantel_tests(in_mantel_dt = mantel_dt)
  ),
  
  tar_target(
    basemaps,
    download_basemap(out_path = file.path(resdir, 'basemap_data'))
  ),
  
  tar_target(
    net_directed,
    direct_network(in_net = net_formatted,
                   idcol = 'OBJECTID_1',
                   outlet_id = 245
                   )
  )
  ,
  
  tar_target(
    hillshade_bolivia,
    create_hillshade(in_dem = unserialize(basemaps$elev_bolivia),
                     z_exponent = 1.3)
  ),
  
  tar_target(
    hillshade_net,
    create_hillshade(in_dem = unserialize(basemaps$elev_net),
                     z_exponent = 1.3)
  )
  ,

  tar_target(
    mapped_caynaca,
    map_caynaca(
      in_spenv_dt = spenv_dt,
      in_net = net_directed,
      in_sites_path = sites_path,
      in_basemaps = basemaps,
      in_hillshade_bolivia = hillshade_bolivia,
      in_hillshade_net = hillshade_net,
      out_plot = file.path(resdir, 'figures_manuscript', 'study_map.png')
    )
  ),
  
  tar_target(
    mapped_data,
    map_data(in_spenv_dt = spenv_dt, 
             in_net = net_directed,
             in_sites_path = sites_path,
             sites_IDcol = sites_IDcol)
  )
)


#Later on: Compute topographic distance (terrain is accidented)

