program offlinediag

  use netcdf
  use input_mod, only:     read_data,                open_input_files,        &
                      get_dimensions,                       axes_init,        &
              read_surf_geopotential

  use get_grid_fields_mod, only:                                              &
                    grid_fields_init,                 get_grid_fields,        &
                       surf_gradient,             zonal_spectrum_grid

  use vars

  use isentropic_mod, only:                                                   &
                     isentropic_init,             sigma_to_isentropes,        &
                isentropic_variables,      compute_isentropic_density

  use local_utilities_mod, only:                                              &
                           sigma_avg,        interpolate_half_to_full,        &
                           mrdnl_avg,                   theta_sfc_pdf,        &
               residual_circulations,              buoyancy_frequency,        &
                      divide_by_psfc,            update_stat_unstable,        &
                      pot_temp_deriv,                  vrtcl_gradient,        &
                      mrdnl_gradient,                       vrtcl_int,        &
                             pdf_avg,                    sigma_avg_3d,        &
                   divide_by_psfc_3d        

  use constants_and_switches_mod, only:                                       &
            reference_sea_level_pres,                              cp,        &
                               rdgas,                           kappa,        &
                                  pi,                           omega,        &
                                grav,                          radius,        &
                               rvgas

  use output_mod, only:                                                       &
                  create_output_file,                 write_variables,        &
                   close_output_file


  implicit none

! ------------------------------ local variables ------------------------------

  real :: time                 ! current time
  integer :: time_counter=1    ! used to calculate frequencies based on daily data
  
  ! local dimensions of gridpoint arrays
  integer ::                                                                  &
       num_lat,            &   ! total number of latitudes
       num_lon,            &   ! total number of longitudes
       num_lev,            &   ! total number of levels
       num_times,          &   ! total number of times in file
       curr_num_times          ! number of times in current segment

  integer i,j,k

! ------------------------------ executable code -----------------------------

  call read_namelist

  call open_input_files(      DataIn,                     Data_Source,        &
                            moisture,                       bucket,           &
                             virtual,                       isentrope,        & 
                            UVarName,                        VVarName,        &
                          VORVarName,                      DIVVarName,        &
                         TEMPVarName,                       TSVarName,        &
                           PSVarName,                                         &
                         SHumVarName,                                         &
                         CondVarName,                     ConvVarName,        &        
                         DiffVarName,                                         &
                   PrecipCondVarName,               PrecipConvVarName,        &
                    SfcFluxLHVarName,                SfcFluxSHVarName,        &
                    SfcFluxSWVarName,               SfcFluxLWDVarName,        &
                   SfcFluxLWUVarName,                SfcFluxQFVarName,        &
                    ToaFluxSWVarName,               ToaFluxLWUVarName,        &
                     DiabCondVarName,                 DiabConvVarName,        &
                     DiabDiffVarName,                  DiabRadVarName,        &
                       DiabSWVarName,                                         &  
                  BucketDepthVarName,          BucketDepthConvVarName,        &
              BucketDepthCondVarName,            BucketDepthLHVarName,        &       
                                               BucketDiffusionVarName,        &                  
                       DragMOVarName,                   DragLHVarName,        &
                       DragSHVarName,                                         &
                       InputFileName)


  call get_dimensions(       num_lat,                         num_lon,        &
                             num_lev,                       num_times )

  ! initialize fields based on dimensions
  call allocate_and_initialize(                                               &
                             num_lat,                         num_lon,        &
                             num_lev )

! get axes
  call axes_init(            deg_lat,                         deg_lon,        &
                             sin_lat,                                         &
                          cos_lat_xy,                      cos_lat_yz )

! read_surf_geopotential must be called before grid_fields_init
  if(.not.smooth_surface)  call read_surf_geopotential(                       &
                           SGVarName,                 SGInputFileName,        &
                           geopot_sfc )

! initialize get_grid_fields module
  call grid_fields_init(      DataIn,                     is_gaussian,        &
                            moisture,                                         &
                    moist_isentropes,                     data_source,        &
                         num_fourier,                         num_lat,        &
                             num_lon,                         num_lev,        &
                          cos_lat_xy,                      geopot_sfc,        &
                               sigma,                              pk,        &
                                  bk )


! initialize isentropic module
  call isentropic_init(   PotTempMin,                      PotTempMax,        &
                      pot_temp_coord,                d_pot_temp_coord )

  
  curr_num_times = num_times / num_segments
  write(*,*) "Starting time ", TimeIn
  write(*,*) 'Analyzing ', curr_num_times, ' instants.'
  
  time = TimeIn
  do while (time < TimeIn + curr_num_times)

     time = time + 1.0

     write(6,*) time

     call read_data(            time,                          u_grid,        &
                              v_grid,                        vor_grid,        &
                            div_grid,                       temp_grid,        &
                             ts_grid,                                         &
                             ps_grid,                       shum_grid,        &
                   dt_shum_cond_grid,               dt_shum_conv_grid,        & 
                   dt_shum_diff_grid,                                         &
                    precip_cond_grid,                precip_conv_grid,        &
                    sfc_flux_lh_grid,                sfc_flux_sh_grid,        &
                    sfc_flux_sw_grid,               sfc_flux_lwd_grid,        &
                   sfc_flux_lwu_grid,             sfc_flux_ocean_grid,        &
                    toa_flux_sw_grid,               toa_flux_lwu_grid,        &
                   dt_temp_cond_grid,               dt_temp_conv_grid,        &
                   dt_temp_diff_grid,                dt_temp_rad_grid,        &
                     dt_temp_sw_grid,                                         &
                   bucket_depth_grid,          bucket_depth_conv_grid,        &
              bucket_depth_cond_grid,            bucket_depth_lh_grid,        &
                                                bucket_diffusion_grid,        &
                  drag_coeff_mo_grid,              drag_coeff_lh_grid,        & 
                  drag_coeff_sh_grid)

     ! calculate virtual temperature consistent with the approximation 
     ! used in the fms code; for use in z, mont_x_grid, N^2, svol
     if (virtual .and. moisture) then
      virtual_temp_grid = temp_grid*(1.0+shum_grid*(rvgas/rdgas-1.0))
     else
      virtual_temp_grid = temp_grid
     endif

     call get_grid_fields(    u_grid,                          v_grid,        &
                            vor_grid,                        div_grid,        &
                           temp_grid,                       shum_grid,        &
                             ts_grid,                                         &                            
                             ps_grid,                       vcos_grid,        &
                              w_grid,                   pot_temp_grid,        &
                     pot_temp_e_grid,             pot_temp_e_sat_grid,        &
                              p_half,                          p_full,        &
                           mont_grid,                    dx_mont_grid,        &
                          sfctn_grid,                           sigma,        &
                       pot_temp_spec,              eape2eke_conv_spec,        &
                           shum_spec,                   sat_shum_spec,        &
                           rhum_spec,                 pot_temp_e_spec,        &
                                                    virtual_temp_grid,        &
                         geopot_grid,                     sat_mr_grid,        &
                       sat_shum_grid,                       rhum_grid)



     ! average (around latitude circle) sigma-coordinate fields (inc. surface)
     call zonal_average_sigma


     if(isentrope) then


     call isentropic_variables(                                               &
                             sin_lat,                   pot_temp_grid,        &
                              p_full,           p_half(:,:,num_lev+1),        &
                            vor_grid,                       vcos_grid,        &
                               k_sfc,                          k_asfc,        &
                     p_isent_grid_nw,           abs_vor_isent_grid_nw,        &
                  vcos_isent_grid_nw)                                         


     pot_vor_isent_grid             = abs_vor_isent_grid_nw
     mrdnl_pot_vor_flux_isent_grid  = vcos_isent_grid_nw * pot_vor_isent_grid
   

!     call theta_sfc_pdf(k_sfc,                   sfc_pot_temp_pdf)

     ! find isentropic density implied by linear interpolation 
     ! routine sigma_to_isentropes
     call compute_isentropic_density(p_isent_grid_nw, dens_isent_grid_nw)


     ! note sfctn_grid is an integral quantity and sfctn_isentr thus has
     ! no density factor (optional argument integral is .true.)

     call sigma_to_isentropes(                                                &
                        pot_temp_grid,                          p_half,       &
                      p_isent_grid_nw,                          k_asfc,       &
                               u_grid,                    u_isent_grid)

     call sigma_to_isentropes(                                                &
                        pot_temp_grid,                          p_half,       &
                      p_isent_grid_nw,                          k_asfc,       &
                      u_grid * u_grid,                u_var_isent_grid)
  
     call sigma_to_isentropes(                                                &
                        pot_temp_grid,                          p_half,       &
                      p_isent_grid_nw,                          k_asfc,       &
                               v_grid,                    v_isent_grid)

     call sigma_to_isentropes(                                                &
                        pot_temp_grid,                          p_half,       &
                      p_isent_grid_nw,                          k_asfc,       &
                            vcos_grid,                 vcos_isent_grid)

     call sigma_to_isentropes(                                                &
                        pot_temp_grid,                          p_half,       &
                      p_isent_grid_nw,                          k_asfc,       &
                      v_grid * v_grid,                v_var_isent_grid)


     call sigma_to_isentropes(                                                &
                        pot_temp_grid,                          p_half,       &
                      p_isent_grid_nw,                          k_asfc,       &
                           sfctn_grid,                sfctn_isent_grid,       &
                              .true. )


     call sigma_to_isentropes(                                                &
                        pot_temp_grid,                          p_half,       &
                      p_isent_grid_nw,                          k_asfc,       &
                            mont_grid,                 mont_isent_grid)


     call sigma_to_isentropes(                                                &
                        pot_temp_grid,                          p_half,       &
                      p_isent_grid_nw,                          k_asfc,       &
                mont_grid * mont_grid,              mont_var_isent_grid)

     call sigma_to_isentropes(                                                &
                       pot_temp_grid,                          p_half,        &
                     p_isent_grid_nw,                           k_sfc,        &
                        dx_mont_grid,                dx_mont_isent_grid)

     if(moisture) then

        call sigma_to_isentropes(                                             &
                       pot_temp_grid,                          p_half,        &
                     p_isent_grid_nw,                           k_sfc,        &
                           shum_grid,                  shum_isent_grid)

        call sigma_to_isentropes(                                             &
                       pot_temp_grid,                          p_half,        &
                     p_isent_grid_nw,                           k_sfc,        &
                  u_grid * shum_grid,       zonal_shum_flux_isent_grid)

        call sigma_to_isentropes(                                             &
                       pot_temp_grid,                          p_half,        &
                     p_isent_grid_nw,                           k_sfc,        &
               vcos_grid * shum_grid,       mrdnl_shum_flux_isent_grid)

        call sigma_to_isentropes(                                             &
                       pot_temp_grid,                          p_half,        &
                     p_isent_grid_nw,                           k_sfc,        &
                           rhum_grid,                  rhum_isent_grid)

        call sigma_to_isentropes(                                             &
                       pot_temp_grid,                          p_half,        &
                     p_isent_grid_nw,                           k_sfc,        &
                       sat_shum_grid,              sat_shum_isent_grid)

        call sigma_to_isentropes(                                             &
                       pot_temp_grid,                          p_half,        &
                     p_isent_grid_nw,                           k_sfc,        &
              u_grid * sat_shum_grid,   zonal_sat_shum_flux_isent_grid)
     
        call sigma_to_isentropes(                                             &
                       pot_temp_grid,                          p_half,        &
                     p_isent_grid_nw,                           k_sfc,        &
           vcos_grid * sat_shum_grid,   mrdnl_sat_shum_flux_isent_grid)

     endif

    !S (1891)
    do i = 1,size(pot_temp_grid,1)
       do k = 1, size(pot_temp_grid,3)
          v_geostr_isent_grid(i, :, k) = dx_mont_isent_grid(i, :, k)  * 1.0 / (2. * omega * sin_lat)
       enddo
    enddo

     ! average (around latitude circle) theta-coordinate fields (inc. surface)
     call zonal_average_theta

     endif


  enddo ! end loop over times

  call time_average_fields(time-TimeIn)

  if(eddy_fluxes) then
     call eddy_fields
  endif

  ! compute mean potential temperature gradients
  call pot_temp_deriv(                                                        &
                             deg_lat,                           sigma,        &
                        pot_temp_avg,                          ps_avg,        &
                   d_dy_pot_temp_avg,             d_dsig_pot_temp_avg,        &
                   d_dp_pot_temp_avg )

  ! get Brunt-Vaeisaelae frequency
  buoyancy_freq_avg = buoyancy_frequency(                                     &
                          p_full,                    virtual_temp_avg,        &
                          pot_temp_avg,              d_dp_pot_temp_avg )

  ! compute mean equivalent potential temperature gradients
  if (moisture)                                                               &
    call pot_temp_deriv(                                                      &
                                  deg_lat,                        sigma,      &
                           pot_temp_e_avg,                       ps_avg,      &
                      d_dy_pot_temp_e_avg,        d_dsig_pot_temp_e_avg,      &
                      d_dp_pot_temp_e_avg )


  ! generate TEM and Modified TEM circulations
  call residual_circulations(                                                 &
                           sfctn_avg,    mrdnl_eddy_pot_temp_flux_avg,        &
        vrtcl_eddy_pot_temp_flux_avg,               d_dy_pot_temp_avg,        &
                 d_dsig_pot_temp_avg,                          ps_avg,        &
                        TEM_res_circ,                   MTEM_res_circ )

  call create_output_file(   num_lat,                         num_lev,        &
                             num_lon                                   )

  call write_variables(      num_lat,                     num_fourier,        &
                             num_lev,                     time-TimeIn,        &
                             num_lon                                   )

  call close_output_file()
  
contains
  
! ----------------------------------------------------------------------------

  subroutine zonal_average_sigma

!                            --- local variables ---

    real, dimension(size(u_grid,1),size(u_grid,2)) :: u_barotr, v_barotr
    real, dimension(size(u_grid,1),size(u_grid,2),size(u_grid,3)) :: ucos_grid ! u * cosine(latitude)

!                            --- executable code ---

! Pressure and density

    ps_avg          = ps_avg           + ps_grid 
    ps_var_avg      = ps_var_avg       + ps_grid*ps_grid
    p_mean           = p_mean            + p_full
    p_var            = p_var             + p_full*p_full

    ts_avg           = ts_avg            + ts_grid

    specific_vol_avg = specific_vol_avg  + sigma_avg_3d(rdgas*virtual_temp_grid/p_full,  &
         ps_grid)

! Zonal Winds

    u_avg            = u_avg             + sigma_avg_3d(u_grid, ps_grid)
    u_barotr           = vrtcl_int(p_half, u_grid)
    u_barotr_avg       = u_barotr_avg        + u_barotr
    u_barotr_var_avg   = u_barotr_var_avg    + u_barotr**2 / ps_grid
    u_var_avg        = u_var_avg         + sigma_avg_3d(u_grid**2, ps_grid)
    mrdnl_u_flux_avg = mrdnl_u_flux_avg  + sigma_avg_3d(u_grid*vcos_grid, ps_grid)
    vrtcl_u_flux_avg = vrtcl_u_flux_avg  + sigma_avg_3d(u_grid*w_grid, ps_grid)
    do k = 1, num_lev
       ucos_grid(:,:,k) = u_grid(:,:,k) * cos_lat_xy
    end do

! Meridional Winds

    v_avg               = v_avg             + sigma_avg_3d(v_grid, ps_grid)
    v_barotr           = vrtcl_int(p_half, v_grid)
    v_barotr_avg       = v_barotr_avg      + v_barotr
    v_barotr_var_avg   = v_barotr_var_avg  + v_barotr**2 / ps_grid
    vcos_avg           = vcos_avg          + sigma_avg_3d(vcos_grid, ps_grid)
    vort_avg           = vort_avg          + sigma_avg_3d(vor_grid, ps_grid)
    vort_var_avg       = vort_var_avg      + sigma_avg_3d(vor_grid**2, ps_grid)
    v_var_avg          = v_var_avg         + sigma_avg_3d(v_grid**2, ps_grid)


    sfctn_avg        = sfctn_avg         +                                    &
        interpolate_half_to_full(sfctn_grid)
    sfctn_var_avg    = sfctn_var_avg     + sfctn_avg**2

! Vertical Winds

    w_avg            = w_avg             + sigma_avg_3d(w_grid, ps_grid)
    w_var_avg        = w_var_avg         + sigma_avg_3d(w_grid**2, ps_grid)

! Temp Quantities

    temp_avg              = temp_avg            + sigma_avg_3d(temp_grid, ps_grid)
    temp_var_avg          = temp_var_avg        + sigma_avg_3d(temp_grid**2, ps_grid)
    zonal_temp_flux_avg   = zonal_temp_flux_avg +                               &
         sigma_avg_3d(u_grid * temp_grid, ps_grid)
    mrdnl_temp_flux_avg   = mrdnl_temp_flux_avg +                               &
         sigma_avg_3d(vcos_grid * temp_grid, ps_grid)
    vrtcl_temp_flux_avg   = vrtcl_temp_flux_avg +                               &
         sigma_avg_3d(w_grid * temp_grid, ps_grid)
    virtual_temp_avg = virtual_temp_avg  + sigma_avg_3d(virtual_temp_grid, ps_grid)

! Potential Temp Quantities

    pot_temp_avg              = pot_temp_avg            + sigma_avg_3d(pot_temp_grid, ps_grid)
    pot_temp_var_avg          = pot_temp_var_avg        + sigma_avg_3d(pot_temp_grid**2, ps_grid)
    zonal_pot_temp_flux_avg   = zonal_pot_temp_flux_avg +                               &
         sigma_avg_3d(u_grid * pot_temp_grid, ps_grid)
    mrdnl_pot_temp_flux_avg   = mrdnl_pot_temp_flux_avg +                               &
         sigma_avg_3d(vcos_grid * pot_temp_grid, ps_grid)
    vrtcl_pot_temp_flux_avg   = vrtcl_pot_temp_flux_avg +                               &
         sigma_avg_3d(w_grid * pot_temp_grid, ps_grid)

! Geopotential Quantities

    z_avg                = z_avg          + sigma_avg_3d(geopot_grid / grav, ps_grid)
    z_var_avg            = z_var_avg      + sigma_avg_3d(geopot_grid**2 / grav**2, ps_grid)
    zonal_z_flux_avg     = zonal_z_flux_avg   +                                 &
         sigma_avg_3d(u_grid * geopot_grid / grav, ps_grid)
    mrdnl_z_flux_avg     = zonal_z_flux_avg   +                                 &
         sigma_avg_3d(vcos_grid * geopot_grid / grav, ps_grid)
    vrtcl_z_flux_avg     = vrtcl_z_flux_avg   +                                 &
         sigma_avg_3d(w_grid * geopot_grid / grav, ps_grid)

! Moisture Quantities
    if(moisture) then
       shum_avg         = shum_avg                + sigma_avg_3d(shum_grid, ps_grid)
       shum_var_avg     = shum_var_avg            + sigma_avg_3d(shum_grid**2, ps_grid)
       zonal_shum_flux_avg = zonal_shum_flux_avg  + sigma_avg_3d(shum_grid*u_grid, ps_grid)
       mrdnl_shum_flux_avg = mrdnl_shum_flux_avg  + sigma_avg_3d(shum_grid*vcos_grid, ps_grid)
       vrtcl_shum_flux_avg = vrtcl_shum_flux_avg  + sigma_avg_3d(shum_grid*w_grid, ps_grid)  

       sat_shum_avg     = sat_shum_avg       + sigma_avg_3d(sat_shum_grid, ps_grid) 
       sat_shum_var_avg = sat_shum_var_avg   + sigma_avg_3d(sat_shum_grid**2, ps_grid)

       rhum_avg         = rhum_avg           + sigma_avg_3d(rhum_grid, ps_grid)
!       rhum_pdf         = rhum_pdf           + pdf_avg_3d(rhum_grid, ps_grid, rhum_bin)

       pot_temp_e_avg      = pot_temp_e_avg        + sigma_avg_3d(pot_temp_e_grid, ps_grid)
       pot_temp_e_sat_avg  = pot_temp_e_sat_avg    + sigma_avg_3d(pot_temp_e_sat_grid, ps_grid)
       pot_temp_e_var_avg  = pot_temp_e_var_avg    + sigma_avg_3d(pot_temp_e_grid**2, ps_grid)
       zonal_pot_temp_e_flux_avg   = zonal_pot_temp_e_flux_avg +                     &
           sigma_avg_3d(u_grid * pot_temp_e_grid, ps_grid)
       mrdnl_pot_temp_e_flux_avg   = mrdnl_pot_temp_e_flux_avg +                     &
           sigma_avg_3d(vcos_grid * pot_temp_e_grid, ps_grid)
       vrtcl_pot_temp_e_flux_avg   = vrtcl_pot_temp_e_flux_avg +                     &
           sigma_avg_3d(w_grid * pot_temp_e_grid, ps_grid)

       precip_cond_avg  = precip_cond_avg + precip_cond_grid
       precip_conv_avg  = precip_conv_avg + precip_conv_grid
       dt_shum_cond_avg   = dt_shum_cond_avg + sigma_avg_3d(dt_shum_cond_grid, ps_grid)
       dt_shum_conv_avg   = dt_shum_conv_avg + sigma_avg_3d(dt_shum_conv_grid, ps_grid)
       dt_shum_diff_avg   = dt_shum_diff_avg + sigma_avg_3d(dt_shum_diff_grid, ps_grid)
       dt_shum_cond_counter = 0
       where(dt_shum_cond_grid.ne.0)
       dt_shum_cond_counter = 1
       endwhere
       dt_shum_conv_counter = 0
       where(dt_shum_conv_grid.ne.0)
       dt_shum_conv_counter = 1
       endwhere
       shum_cond_prob_avg   = shum_cond_prob_avg + sigma_avg_3d(dt_shum_cond_counter, ps_grid)
       shum_conv_prob_avg   = shum_conv_prob_avg + sigma_avg_3d(dt_shum_conv_counter, ps_grid)


       ! calculate precipitation intensity based on daily data
       ! (input precipitation data should be 4xday and time averaged rather than instantaneous)

       !!!!!!MODEL TIME SAVING STEP DEPENDENT: FACTORS 0.25;4

       precip_cond_daily = precip_cond_daily + 0.25*precip_cond_grid
       precip_conv_daily = precip_conv_daily + 0.25*precip_conv_grid
       precip_tot_daily = precip_tot_daily + 0.25*precip_cond_grid + 0.25*precip_conv_grid
       
       if (time_counter == 4) then  ! reached end of day
          days_total = days_total + one_array

          precip_mask = (precip_cond_daily+precip_conv_daily) .gt. precip_daily_threshold
          precip_daily_above_threshold_avg = precip_daily_above_threshold_avg + (precip_cond_daily+precip_conv_daily)*precip_mask                           
          days_above_threshold = days_above_threshold + one_array*precip_mask     

          ! Computes a pdf of excedance frequency
          precip_cond_pdf  = precip_cond_pdf    + pdf_avg(precip_cond_daily, precip_bin)
          precip_conv_pdf  = precip_conv_pdf    + pdf_avg(precip_conv_daily, precip_bin)
          precip_tot_pdf  = precip_tot_pdf    + pdf_avg(precip_tot_daily, precip_bin)

          precip_cond_daily = 0
          precip_conv_daily = 0
          precip_tot_daily  = 0
          time_counter      = 1       
       else       
          time_counter = time_counter + 1       
       endif


! Surface and TOA Fluxes

       sfc_flux_lh_avg    = sfc_flux_lh_avg     + sfc_flux_lh_grid
       sfc_flux_sh_avg    = sfc_flux_sh_avg     + sfc_flux_sh_grid
       sfc_flux_sw_avg    = sfc_flux_sw_avg     + sfc_flux_sw_grid
       sfc_flux_lwd_avg   = sfc_flux_lwd_avg    + sfc_flux_lwd_grid
       sfc_flux_lwu_avg   = sfc_flux_lwu_avg    + sfc_flux_lwu_grid     
       sfc_flux_ocean_avg = sfc_flux_ocean_avg  + sfc_flux_ocean_grid
       toa_flux_sw_avg    = toa_flux_sw_avg     + toa_flux_sw_grid
       toa_flux_lwu_avg   = toa_flux_lwu_avg    + toa_flux_lwu_grid

       drag_coeff_mo_avg  = drag_coeff_mo_avg   + drag_coeff_mo_grid
       drag_coeff_lh_avg  = drag_coeff_lh_avg   + drag_coeff_lh_grid     
       drag_coeff_sh_avg  = drag_coeff_sh_avg   + drag_coeff_sh_grid


! bucket (simplified hydrology) 
       if(bucket) then
          bucket_depth_avg        = bucket_depth_avg   + bucket_depth_grid
          bucket_depth_conv_avg   = bucket_depth_conv_avg   + bucket_depth_conv_grid
          bucket_depth_cond_avg   = bucket_depth_cond_avg   + bucket_depth_cond_grid
          bucket_depth_LH_avg   = bucket_depth_LH_avg   + bucket_depth_LH_grid
          bucket_diffusion_avg    = bucket_diffusion_avg   + bucket_diffusion_grid
       endif

! Heating Rates

       dt_temp_cond_avg = dt_temp_cond_avg + sigma_avg_3d(dt_temp_cond_grid, ps_grid)     
       dt_temp_conv_avg = dt_temp_conv_avg + sigma_avg_3d(dt_temp_conv_grid, ps_grid)        
       dt_temp_diff_avg = dt_temp_diff_avg + sigma_avg_3d(dt_temp_diff_grid, ps_grid)        
       dt_temp_rad_avg  = dt_temp_rad_avg + sigma_avg_3d(dt_temp_rad_grid, ps_grid)     
       dt_temp_sw_avg   = dt_temp_sw_avg + sigma_avg_3d(dt_temp_sw_grid, ps_grid)     

    end if
    
   if(isentrope) then

    mont_sfc  = geopot_sfc + cp * temp_grid(:, :, num_lev)
    call surf_gradient(pot_temp_grid(:,:,num_lev), dx_pot_temp_sfc,                &
         dy_pot_temp_sfc)    

    ! geostrophic surface flux of potential temperature (times f)
     mont_dx_pot_temp_sfc_avg =  mont_dx_pot_temp_sfc_avg  +     &
         mont_sfc * dx_pot_temp_sfc

     do j = 1, num_lat
        mrdnl_geostr_eddy_pot_temp_flux_sfc_avg(:,j) =  - mont_dx_pot_temp_sfc_avg(:,j)  * 1.0 / (2. * omega * sin_lat(j))
     enddo

    pot_temp_spec_avg = pot_temp_spec_avg +  pot_temp_spec
    eape2eke_conv_spec_avg = eape2eke_conv_spec_avg + eape2eke_conv_spec
    shum_spec_avg = shum_spec_avg +  shum_spec
    sat_shum_spec_avg = sat_shum_spec_avg + sat_shum_spec
    rhum_spec_avg = rhum_spec_avg + rhum_spec
    pot_temp_e_spec_avg =  pot_temp_e_spec_avg + pot_temp_e_spec

    endif

 end subroutine zonal_average_sigma
  
! ----------------------------------------------------------------------------

  subroutine zonal_average_theta

     p_isent_avg_nw     = p_isent_avg_nw      + p_isent_grid_nw
     p_var_isent_avg_nw = p_var_isent_avg_nw  + p_isent_grid_nw**2

     exner_isent_avg_nw    = exner_isent_avg_nw  +                                           &
          cp *  (p_isent_grid_nw/reference_sea_level_pres)**kappa  
     
     p_exner_isent_avg_nw    = p_exner_isent_avg_nw +                                         &
          (p_isent_grid_nw/reference_sea_level_pres)**(kappa+1.0) 

     dens_isent_avg_nw     = dens_isent_avg_nw   +  dens_isent_grid_nw

     u_isent_avg           = u_isent_avg       + u_isent_grid
     u_var_isent_avg       = u_var_isent_avg   + u_var_isent_grid

     v_isent_avg           = v_isent_avg       + v_isent_grid
     vcos_isent_avg        = vcos_isent_avg    + vcos_isent_grid
     v_var_isent_avg       = v_var_isent_avg   + v_var_isent_grid

     v_geostr_isent_avg    = v_geostr_isent_avg +  v_geostr_isent_grid

     sfctn_isent_avg       = sfctn_isent_avg   + sfctn_isent_grid

     pot_vor_isent_avg     = pot_vor_isent_avg + pot_vor_isent_grid
     mrdnl_pot_vor_flux_isent_avg =  mrdnl_pot_vor_flux_isent_avg + mrdnl_pot_vor_flux_isent_grid 

     mont_isent_avg  = mont_isent_avg + mont_isent_grid
     mont_var_isent_avg  = mont_var_isent_avg + mont_var_isent_grid

     if(moisture)   then
     ! note linearly interpolated quantities already include density weighting
        shum_isent_avg      = shum_isent_avg   + shum_isent_grid
        zonal_shum_flux_isent_avg  = zonal_shum_flux_isent_avg  + zonal_shum_flux_isent_grid
        mrdnl_shum_flux_isent_avg  = mrdnl_shum_flux_isent_avg  + mrdnl_shum_flux_isent_grid

        rhum_isent_avg      = rhum_isent_avg   + rhum_isent_grid

        sat_shum_isent_avg = sat_shum_isent_avg + sat_shum_isent_grid
        zonal_sat_shum_flux_isent_avg = zonal_sat_shum_flux_isent_avg  + zonal_sat_shum_flux_isent_grid
        mrdnl_sat_shum_flux_isent_avg = mrdnl_sat_shum_flux_isent_avg  + mrdnl_sat_shum_flux_isent_grid
     endif

  end subroutine zonal_average_theta

! ----------------------------------------------------------------------------

  subroutine time_average_fields(dof)

!----------------------------- input arguments --------------------------------

    real, intent(in) :: dof ! degrees of freedom

!----------------------------- local variables  -------------------------------
    
    integer j

!----------------------------- executable code  -------------------------------
    
     ! divide averages by number of degrees of freedom in estimates



    ! divide sigma-averaged fields by mean sfc pressure and dof's
    p_mean               = p_mean/ dof
    p_var                = p_var / dof
    ts_avg               = ts_avg / dof
    ps_avg               = ps_avg / dof  
    ps_var_avg           = ps_var_avg / dof
    specific_vol_avg     = divide_by_psfc_3d(specific_vol_avg/dof,     ps_avg)

    u_avg                = divide_by_psfc_3d(u_avg/dof,              ps_avg)
    u_var_avg            = divide_by_psfc_3d(u_var_avg/dof,          ps_avg)
    mrdnl_u_flux_avg     = divide_by_psfc_3d(mrdnl_u_flux_avg/dof,   ps_avg)
    vrtcl_u_flux_avg     = divide_by_psfc_3d(vrtcl_u_flux_avg/dof,   ps_avg)
    u_barotr_avg         = u_barotr_avg / ( dof * ps_avg )
    v_barotr_avg         = v_barotr_avg / ( dof * ps_avg )
    u_barotr_var_avg     = u_barotr_var_avg / (dof * ps_avg )
    v_barotr_var_avg     = v_barotr_var_avg / (dof * ps_avg )

    v_avg                = divide_by_psfc_3d(v_avg/dof,              ps_avg)
    v_var_avg            = divide_by_psfc_3d(v_var_avg/dof,          ps_avg)
    vcos_avg             = divide_by_psfc_3d(vcos_avg/dof,           ps_avg) 
    sfctn_avg            = sfctn_avg            / dof
    sfctn_var_avg        = sfctn_var_avg       / dof
    vort_avg             = divide_by_psfc_3d(vort_avg/dof,           ps_avg)
    vort_var_avg         = divide_by_psfc_3d(vort_var_avg/dof,       ps_avg)
    w_avg                = divide_by_psfc_3d(w_avg/dof,              ps_avg) 
    w_var_avg            = divide_by_psfc_3d(w_var_avg/dof,          ps_avg)

    pot_temp_avg            = divide_by_psfc_3d(pot_temp_avg/dof,          ps_avg)
    pot_temp_var_avg        = divide_by_psfc_3d(pot_temp_var_avg/dof,      ps_avg) 
    zonal_pot_temp_flux_avg = divide_by_psfc_3d(zonal_pot_temp_flux_avg/dof, ps_avg)
    mrdnl_pot_temp_flux_avg = divide_by_psfc_3d(mrdnl_pot_temp_flux_avg/dof, ps_avg)
    vrtcl_pot_temp_flux_avg = divide_by_psfc_3d(vrtcl_pot_temp_flux_avg/dof, ps_avg)
    pot_temp_spec_avg                      = pot_temp_spec_avg / dof
    pot_temp_spec_avg(1,:,:)               = 0.25 * pot_temp_spec_avg(1,:,:)
    pot_temp_spec_avg(2:num_fourier+1,:,:) = 0.25 * pot_temp_spec_avg(2:num_fourier+1,:,:)


    temp_avg             = divide_by_psfc_3d(temp_avg/dof,           ps_avg)
    temp_var_avg         = divide_by_psfc_3d(temp_var_avg/dof,       ps_avg) 
    zonal_temp_flux_avg  = divide_by_psfc_3d(zonal_temp_flux_avg/dof,  ps_avg)
    mrdnl_temp_flux_avg  = divide_by_psfc_3d(mrdnl_temp_flux_avg/dof,  ps_avg)
    vrtcl_temp_flux_avg  = divide_by_psfc_3d(vrtcl_temp_flux_avg/dof,  ps_avg)
    virtual_temp_avg     = divide_by_psfc_3d(virtual_temp_avg/dof,   ps_avg)

    z_avg                = divide_by_psfc_3d(z_avg/dof,              ps_avg)
    z_var_avg            = divide_by_psfc_3d(z_var_avg/dof,          ps_avg)
    mrdnl_z_flux_avg     = divide_by_psfc_3d(zonal_z_flux_avg/dof,     ps_avg)
    mrdnl_z_flux_avg     = divide_by_psfc_3d(mrdnl_z_flux_avg/dof,     ps_avg)
    vrtcl_z_flux_avg     = divide_by_psfc_3d(vrtcl_z_flux_avg/dof,     ps_avg)

    eape2eke_conv_spec_avg                      = eape2eke_conv_spec_avg  / dof
    eape2eke_conv_spec_avg(1,:,:)               = 0.25 * eape2eke_conv_spec_avg(1,:,:)
    eape2eke_conv_spec_avg(2:num_fourier+1,:,:) = 0.25 * eape2eke_conv_spec_avg(2:num_fourier+1,:,:)

    if(moisture) then
       shum_avg                 = divide_by_psfc_3d(shum_avg/dof,               ps_avg)
       shum_var_avg             = divide_by_psfc_3d(shum_var_avg/dof,           ps_avg)
       zonal_shum_flux_avg      = divide_by_psfc_3d(zonal_shum_flux_avg/dof,    ps_avg)
       mrdnl_shum_flux_avg      = divide_by_psfc_3d(mrdnl_shum_flux_avg/dof,    ps_avg)
       vrtcl_shum_flux_avg      = divide_by_psfc_3d(vrtcl_shum_flux_avg/dof,    ps_avg)
      shum_spec_avg                      = shum_spec_avg  / dof
      shum_spec_avg(1,:,:)               = 0.25 * shum_spec_avg(1,:,:)
      shum_spec_avg(2:num_fourier+1,:,:) = 0.25 * shum_spec_avg(2:num_fourier+1,:,:)

       sat_shum_avg             = divide_by_psfc_3d(sat_shum_avg/dof,           ps_avg)
       sat_shum_var_avg         = divide_by_psfc_3d(sat_shum_var_avg/dof,       ps_avg)
      sat_shum_spec_avg                      = sat_shum_spec_avg  / dof
      sat_shum_spec_avg(1,:,:)               = 0.25 * sat_shum_spec_avg(1,:,:)
      sat_shum_spec_avg(2:num_fourier+1,:,:) = 0.25 * sat_shum_spec_avg(2:num_fourier+1,:,:)

       pot_temp_e_avg              = divide_by_psfc_3d(pot_temp_e_avg/dof,            ps_avg)
       pot_temp_e_var_avg          = divide_by_psfc_3d(pot_temp_e_var_avg/dof,        ps_avg) 
       zonal_pot_temp_e_flux_avg   = divide_by_psfc_3d(zonal_pot_temp_e_flux_avg/dof, ps_avg)
       mrdnl_pot_temp_e_flux_avg   = divide_by_psfc_3d(mrdnl_pot_temp_e_flux_avg/dof, ps_avg)
       vrtcl_pot_temp_e_flux_avg   = divide_by_psfc_3d(vrtcl_pot_temp_e_flux_avg/dof, ps_avg)
       pot_temp_e_spec_avg                      = pot_temp_e_spec_avg / dof
       pot_temp_e_spec_avg(1,:,:)               = 0.25 * pot_temp_e_spec_avg(1,:,:)
       pot_temp_e_spec_avg(2:num_fourier+1,:,:) = 0.25 * pot_temp_e_spec_avg(2:num_fourier+1,:,:)

       pot_temp_e_sat_avg          = divide_by_psfc_3d(pot_temp_e_sat_avg/dof,        ps_avg)

       rhum_avg                           = divide_by_psfc_3d(rhum_avg/dof,               ps_avg)
       rhum_spec_avg                      = rhum_spec_avg  / dof
       rhum_spec_avg(1,:,:)               = 0.25 * rhum_spec_avg(1,:,:)
       rhum_spec_avg(2:num_fourier+1,:,:) = 0.25 * rhum_spec_avg(2:num_fourier+1,:,:)
!       do j=1,num_bin
!          rhum_pdf(:,:,:,j)    = divide_by_psfc_3d(rhum_pdf(:,:,:,j)/dof,        ps_avg)
!       end do

       dt_shum_cond_avg  = divide_by_psfc_3d(dt_shum_cond_avg/dof,   ps_avg)
       dt_shum_conv_avg  = divide_by_psfc_3d(dt_shum_conv_avg/dof,   ps_avg)
       dt_shum_diff_avg  = divide_by_psfc_3d(dt_shum_diff_avg/dof,   ps_avg)
       shum_cond_prob_avg = divide_by_psfc_3d(shum_cond_prob_avg/dof,   ps_avg)
       shum_conv_prob_avg = divide_by_psfc_3d(shum_conv_prob_avg/dof,   ps_avg)
       precip_cond_avg  = precip_cond_avg/dof
       precip_conv_avg  = precip_conv_avg/dof
       precip_daily_above_threshold_avg  = precip_daily_above_threshold_avg/days_above_threshold
       precip_daily_above_threshold_prob_avg = days_above_threshold/days_total

       precip_cond_pdf = 4*precip_cond_pdf/dof
       precip_conv_pdf = 4*precip_conv_pdf/dof
       precip_tot_pdf = 4*precip_tot_pdf/dof

       sfc_flux_lh_avg      = sfc_flux_lh_avg / dof   
       sfc_flux_sh_avg      = sfc_flux_sh_avg / dof   
       sfc_flux_sw_avg      = sfc_flux_sw_avg / dof   
       sfc_flux_lwd_avg     = sfc_flux_lwd_avg / dof   
       sfc_flux_lwu_avg     = sfc_flux_lwu_avg / dof   
       sfc_flux_ocean_avg   = sfc_flux_ocean_avg / dof 
       toa_flux_sw_avg      = toa_flux_sw_avg / dof   
       toa_flux_lwu_avg     = toa_flux_lwu_avg / dof   

       drag_coeff_mo_avg    = drag_coeff_mo_avg / dof
       drag_coeff_lh_avg    = drag_coeff_lh_avg / dof
       drag_coeff_sh_avg    = drag_coeff_sh_avg / dof

       dt_temp_cond_avg     = divide_by_psfc_3d(dt_temp_cond_avg/dof, ps_avg)   
       dt_temp_conv_avg     = divide_by_psfc_3d(dt_temp_conv_avg/dof, ps_avg)   
       dt_temp_diff_avg     = divide_by_psfc_3d(dt_temp_diff_avg/dof, ps_avg)   
       dt_temp_rad_avg      = divide_by_psfc_3d(dt_temp_rad_avg/dof, ps_avg)   
       dt_temp_sw_avg       = divide_by_psfc_3d(dt_temp_sw_avg/dof, ps_avg)   

    endif

    ! added fridoo sept 2012 bucket hydrology
    if(bucket) then
        bucket_depth_avg    = bucket_depth_avg / dof
        bucket_depth_conv_avg    = bucket_depth_conv_avg / dof
        bucket_depth_cond_avg    = bucket_depth_cond_avg / dof
        bucket_depth_LH_avg    = bucket_depth_LH_avg / dof
        bucket_diffusion_avg     = bucket_diffusion_avg / dof
    endif
    ! end fridoo

    if(isentrope) then
       p_isent_avg_nw           = p_isent_avg_nw          / dof
       p_var_isent_avg_nw       = p_var_isent_avg_nw      / dof
       exner_isent_avg_nw       = exner_isent_avg_nw          / dof
       p_exner_isent_avg_nw     = p_exner_isent_avg_nw    / dof
       dens_isent_avg_nw        = dens_isent_avg_nw           / dof
       u_isent_avg              = u_isent_avg              / dof
       u_var_isent_avg          = u_var_isent_avg          / dof
       v_isent_avg              = v_isent_avg              / dof
       vcos_isent_avg           = vcos_isent_avg           / dof
       v_var_isent_avg          = v_var_isent_avg          / dof
       v_geostr_isent_avg       = v_geostr_isent_avg       / dof 
       sfctn_isent_avg          = sfctn_isent_avg          / dof
       pot_vor_isent_avg       = pot_vor_isent_avg         / dof
       pot_vor_var_isent_avg   = pot_vor_var_isent_avg     / dof
       mrdnl_pot_vor_flux_isent_avg = mrdnl_pot_vor_flux_isent_avg / dof
       mont_isent_avg           = mont_isent_avg           / dof
       mont_var_isent_avg       = mont_var_isent_avg       / dof

       mrdnl_geostr_eddy_pot_temp_flux_sfc_avg  = mrdnl_geostr_eddy_pot_temp_flux_sfc_avg / dof
!       sfc_pot_temp_pdf = sfc_pot_temp_pdf / dof

       if(moisture) then
          shum_isent_avg            = shum_isent_avg              / dof
          zonal_shum_flux_isent_avg  = zonal_shum_flux_isent_avg  / dof
          mrdnl_shum_flux_isent_avg  = mrdnl_shum_flux_isent_avg  / dof
          rhum_isent_avg            = rhum_isent_avg              / dof
          sat_shum_isent_avg        = sat_shum_isent_avg          / dof
          zonal_sat_shum_flux_isent_avg = zonal_sat_shum_flux_isent_avg      / dof
          mrdnl_sat_shum_flux_isent_avg = mrdnl_sat_shum_flux_isent_avg      / dof
       endif
  
    endif

  end subroutine time_average_fields

! ----------------------------------------------------------------------------
  
  subroutine eddy_fields

!----------------------------- executable code  -------------------------------


    mrdnl_eddy_u_flux_avg    = mrdnl_u_flux_avg - u_avg * vcos_avg
    vrtcl_eddy_u_flux_avg    = vrtcl_u_flux_avg - w_avg * u_avg

    zonal_eddy_temp_flux_avg  = zonal_temp_flux_avg - u_avg * temp_avg
    mrdnl_eddy_temp_flux_avg  = mrdnl_temp_flux_avg - vcos_avg * temp_avg
    vrtcl_eddy_temp_flux_avg  = vrtcl_temp_flux_avg - w_avg * temp_avg    

    zonal_eddy_pot_temp_flux_avg  = zonal_pot_temp_flux_avg - u_avg * pot_temp_avg
    mrdnl_eddy_pot_temp_flux_avg  = mrdnl_pot_temp_flux_avg - vcos_avg * pot_temp_avg
    vrtcl_eddy_pot_temp_flux_avg  = vrtcl_pot_temp_flux_avg - w_avg * pot_temp_avg

    zonal_eddy_z_flux_avg    = zonal_z_flux_avg - u_avg * z_avg
    mrdnl_eddy_z_flux_avg    = mrdnl_z_flux_avg - vcos_avg * z_avg
    vrtcl_eddy_z_flux_avg    = vrtcl_z_flux_avg - w_avg * z_avg

    if (moisture) then
       zonal_eddy_shum_flux_avg  = zonal_shum_flux_avg - u_avg * shum_avg
       mrdnl_eddy_shum_flux_avg  = mrdnl_shum_flux_avg - vcos_avg * shum_avg
       vrtcl_eddy_shum_flux_avg  = vrtcl_shum_flux_avg - w_avg * shum_avg    

       zonal_eddy_pot_temp_e_flux_avg  = zonal_pot_temp_e_flux_avg - u_avg * pot_temp_e_avg
       mrdnl_eddy_pot_temp_e_flux_avg  = mrdnl_pot_temp_e_flux_avg - vcos_avg * pot_temp_e_avg
       vrtcl_eddy_pot_temp_e_flux_avg  = vrtcl_pot_temp_e_flux_avg - w_avg * pot_temp_e_avg

    end if

    if (isentrope) then     
     where (dens_isent_avg_nw .gt. 1e-6)
          mrdnl_eddy_pot_vor_flux_isent_avg = mrdnl_pot_vor_flux_isent_avg - vcos_isent_avg * pot_vor_isent_avg * 1 / dens_isent_avg_nw

     elsewhere
          mrdnl_eddy_pot_vor_flux_isent_avg =  0.0 
     endwhere

         if(moisture) then
          where (dens_isent_avg_nw .gt. 1e-6)
             zonal_eddy_shum_flux_isent_avg     = zonal_shum_flux_isent_avg - u_isent_avg * shum_isent_avg * 1 / dens_isent_avg_nw
             zonal_eddy_sat_shum_flux_isent_avg = zonal_sat_shum_flux_isent_avg - u_isent_avg * sat_shum_isent_avg * 1 / dens_isent_avg_nw
             mrdnl_eddy_shum_flux_isent_avg     = mrdnl_shum_flux_isent_avg - vcos_isent_avg * shum_isent_avg * 1 / dens_isent_avg_nw
             mrdnl_eddy_sat_shum_flux_isent_avg = mrdnl_sat_shum_flux_isent_avg - vcos_isent_avg * sat_shum_isent_avg * 1 / dens_isent_avg_nw
          elsewhere 
             zonal_eddy_shum_flux_isent_avg     = 0.0 
             zonal_eddy_sat_shum_flux_isent_avg = 0.0 
             mrdnl_eddy_shum_flux_isent_avg     = 0.0 
             mrdnl_eddy_sat_shum_flux_isent_avg = 0.0 
           endwhere
         endif
     endif


  end subroutine eddy_fields

! ----------------------------------------------------------------------------

  subroutine read_namelist

    integer :: unit

    unit = 20
    rewind(unit)
    open(unit, file='input.nml')
    read(unit, nml=main_list)
    rewind(unit)
    read(unit, nml=filename_list)

    close(unit)


  end subroutine read_namelist

end program offlinediag
  

