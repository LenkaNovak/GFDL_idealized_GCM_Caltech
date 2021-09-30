program offlinediag

!This has been adapted to extract the bare minimum for the MWM analysis
!
! sigma levs: [0.05501126, 0.0873321 , 0.13321935, 0.19526839, 0.27502227,
!       0.37219895, 0.48400947, 0.60478926, 0.72614911, 0.83775905,
!       0.92871842, 0.98928317]

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
       num_lev_out,         &   ! total number of levels of the output (set manually below and in the output.f90 file)
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
!  call isentropic_init(   PotTempMin,                      PotTempMax,        &
!                      pot_temp_coord,                d_pot_temp_coord )

! average (around latitude circle) sigma-coordinate fields (inc. surface)
call zonal_average_sigma

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
                            ps_grid,                       shum_grid)

     ! calculate virtual temperature consistent with the approximation 
     ! used in the fms code; for use in z, mont_x_grid, N^2, svol
     !if (virtual .and. moisture) then
     virtual_temp_grid = temp_grid*(1.0+shum_grid*(rvgas/rdgas-1.0))
     !else
     ! virtual_temp_grid = temp_grid
     !endif

call get_grid_fields(   u_grid,                         v_grid,        &
                        vor_grid,                       div_grid,      &
                        temp_grid,                      shum_grid,     &
                        ts_grid,                                       &
                        ps_grid,                                       &
                        w_grid,                         pot_temp_grid, &
                        pot_temp_e_grid,                pot_temp_e_sat_grid, &
                        p_half,                         p_full,           &
                        virtual_temp_grid,              sigma,            &
                        sat_mr_grid,                                      &
                        sat_shum_grid,                  rhum_grid)

    ! average (around latitude circle) sigma-coordinate fields (inc. surface)
    call zonal_average_sigma



  enddo ! end loop over times

  call time_average_fields(time-TimeIn)


  !if(eddy_fluxes) then
  !   call eddy_fields
  !endif

  num_lev_out = int(12)
  call create_output_file(   num_lat,                         num_lev_out,        &
                             num_lon                                   )

  call write_variables(      num_lat,                     num_fourier,        &
                             num_lev,                     time-TimeIn,        &
                             num_lon                                   )

  call close_output_file()
  
contains
  
! ----------------------------------------------------------------------------

  subroutine zonal_average_sigma

!                            --- local variables ---


!                            --- executable code ---

! Pressure and density

    ps_avg          = ps_avg           + ps_grid 
    ts_avg           = ts_avg            + ts_grid


! Zonal Winds

    u_avg            = u_avg             + sigma_avg_3d(u_grid, ps_grid)


! Meridional Winds

    v_avg               = v_avg             + sigma_avg_3d(v_grid, ps_grid)



! Vertical Winds

    w_avg            = w_avg             + sigma_avg_3d(w_grid, ps_grid)

! Temp Quantities

    temp_avg              = temp_avg            + sigma_avg_3d(temp_grid, ps_grid)
    virtual_temp_avg = virtual_temp_avg  + sigma_avg_3d(virtual_temp_grid, ps_grid)

! Potential Temp Quantities

    pot_temp_avg              = pot_temp_avg            + sigma_avg_3d(pot_temp_grid, ps_grid)



! Moisture Quantities

       pot_temp_e_avg      = pot_temp_e_avg        + sigma_avg_3d(pot_temp_e_grid, ps_grid)





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
  

