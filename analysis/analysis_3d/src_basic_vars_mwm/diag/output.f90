module output_mod

  use vars
  use typesizes, only: bytesizesOK
  use netcdf
  use constants_and_switches_mod, only: gspval

implicit none
private

public :: create_output_file, write_variables, close_output_file

  integer, parameter ::                                                       &
       nvars = 1002,       &  ! number of output fields (max)
       theta_switch = 1,   &  ! theta variable switch
       sigma_switch = 2,   &  ! sigma variable switch
       surface_switch = 3, &  ! surface variable switch
       zon_wave_switch = 4,&  ! zonal wavenumber variable switch
       spectra_switch = 5, &  ! spectral variable switch
       pdf_switch = 6,     &  ! pdf (lat, lev, bin)
       spectra_switch_3d = 7, &  ! 3d spectral (wn, wn, height)
       sigma_switch_3d = 8,   &
       surface_switch_3d = 9, &
       theta_switch_3d = 10,  &
       pdf_switch_3d = 12

  integer ::                                                                  &
       ncFileID,           &  ! file ID
       latDimID,           &  ! latitude file ID
       lonDimID,           &  ! longitude file ID
       sigmaDimID,         &  ! sigma dimension ID
       thetaDimID,         &  ! theta dimension ID
       zon_wavenDimID,     &  ! zonal wavenumber dimension ID 
       legendreDimID,      &  ! legendre dimension ID
       binDimID,           &  ! pdf bin dimension ID
       pbinDimID,          &  ! precip pdf bin dimension ID
       timeDimID,          &  ! time dimension ID
       timeVarID,          &  ! time variable ID
       latVarID,           &  ! latitude variable ID
       lonVarID,           &  ! longitude variable ID
       sigmaVarID,         &  ! sigma variable ID
       thetaVarID,         &  ! theta variable ID
       zon_wavenVarID,     &  ! zonal wavenumber variable ID
       legendreVarID,      &  ! legendre variable ID
       binVarID,           &  ! pdf bin variable ID
       pbinVarID,          &  ! precip pdf bin variable ID
       statVarID(nvars)       ! array of field variable IDs


    character(100) :: fileplace
contains
! #############################################################################

  subroutine create_output_file(                                              &
                             num_lat,                         num_lev_out,        &
                             num_lon                                          )

! creates output file and axis variables
    
!----------------------------- input arguments --------------------------------
    
    integer, intent(in) :: num_lat, num_lev_out, num_lon
    
!----------------------------- local variables --------------------------------

!----------------------------- executable code --------------------------------

    if(.not. byteSizesOK()) then
       stop 'open_and_initialize: Compiler does not appear to support required kinds of variables'
    end if

    ! -----------   create netCDF file and define variables  -------------
    ! create netcdf file
    call check(nf90_create(path = trim(OutputFileName),                       &
         cmode = nf90_share, ncid = ncFileID))

    ! global attributes
    call check(nf90_put_att(ncFileID, nf90_global, "title",                   &
         "OfflineDiag Analyses") )
       
       ! global attributes
    call check(nf90_put_att(ncFileID, nf90_global, "title",                   &
         "OfflineDiag Analyses") )


    ! Define axes
    ! -----------
    ! latitude
    call check(NF90_DEF_DIM(ncFileID, "lat", num_lat, latDimID) )
    call check(nf90_def_var(ncFileID, "lat", nf90_double, latDimID,           &
         latVarID) )  
    call check(nf90_put_att(ncFileID, latVarID,                               &
         "long_name", "latitude") )
    call check(nf90_put_att(ncFileID, latVarID,                               &
         "units", "degree") )
       
    ! longitude
    call check(NF90_DEF_DIM(ncFileID, "lon", num_lon, lonDimID) )
    call check(nf90_def_var(ncFileID, "lon", nf90_double, lonDimID,           &
         lonVarID) )  
    call check(nf90_put_att(ncFileID, lonVarID,                               &
         "long_name", "longitute") )
    call check(nf90_put_att(ncFileID, lonVarID,                               &
         "units", "degrees_east") )

    ! sigma level
    call check(NF90_DEF_DIM(ncFileID, "sigma", num_lev_out, sigmaDimID) )
    call check(nf90_def_var(ncFileID, "sigma", nf90_double, sigmaDimID,       &
         sigmaVarID))
    call check(nf90_put_att(ncFileID, sigmaVarID,                             &
         "long_name", "sigma level") )
    call check(nf90_put_att(ncFileID, sigmaVarID,                             &
         "units", "none") )



    ! time
    call check(NF90_DEF_DIM(ncFileID, "times", 1, timeDimID) )
    call check(nf90_def_var(ncFileID, "times", nf90_double, timeDimID,        &
         timeVarID) )  
    call check(nf90_put_att(ncFileID, timeVarID,                              &
         "long_name", "number of instants") )
       
    call initialize_variables
       
     end subroutine create_output_file
    
! #############################################################################

     subroutine initialize_variables

    ! Define variables
    ! ----------------

    call variable_init("u",                                                   &
         "Zonal Wind",                                                        & 
         "m/s",                                                               &
         sigma_switch_3d,                                                        &
         statVarID(5))


    call variable_init("v",                                                   &
         "Meridional Wind",                                                   &
         "m/s",                                                               &
         sigma_switch_3d,                                                        &
         statVarID(13))

    call variable_init("w",                                                   &
            "Vertical Wind; w = \dot\sigma",                                     &
            "1/s",                                                               &
            sigma_switch_3d,                                                        &
            statVarID(17))

    call variable_init("pot_temp",                                            &
            "Potential Temperature",                                             &
            "K",                                                                 &
            sigma_switch_3d,                                                        &
            statVarID(33))

    call variable_init("pot_temp_eqv",                                     &
            "Equivalent Potential Temperature",                               &
            "K",                                                              &
            sigma_switch_3d,                                                     &
            statVarID(61))
     
     call variable_init("z",                                                   &
            "Geopotential Height",                                               &
            "m",                                                                 &
            sigma_switch_3d,                                                        &
            statVarID(28))


     call variable_init("temp",                                                &
            "Temperature",                                                       &
            "K",                                                                 &
            sigma_switch_3d,                                                        &
            statVarID(22))
     ! end fridooo


      
    call check(nf90_enddef(ncfileID))
    
  end subroutine initialize_variables
  
! #############################################################################

  subroutine variable_init(name, long_name, units, switch, VarID)
  
    character(len=*), intent(in) :: name, long_name, units
    integer, intent(in) :: switch
    integer, intent(out) :: VarID

       if(switch == sigma_switch_3d) then
          call check(nf90_def_var(ncFileID, trim(name), nf90_double,          &
               (/ lonDimID, latDimID, sigmaDimID /), VarID ) )
       endif
       
       call check(nf90_put_att(ncFileID, VarID,                               &
            "long_name", trim(long_name)))
       call check(nf90_put_att(ncFileID, VarID,                               &
            "units", trim(units))) 

     end subroutine variable_init

! #############################################################################

  subroutine write_variables(num_lat,                     num_fourier,        &
                             num_lev_out,                       num_times,         &
                             num_lon                                    )

!----------------------------- input arguments --------------------------------

    integer, intent(in) ::                                                    &
         num_lat,          &   !
         num_lon,          &   !
         num_fourier,      &   !
         num_lev_out               !
    
    real, intent(in) ::                                                       &
         num_times

!----------------------------- local variables --------------------------------

    real    :: zon_waven(num_fourier+1), legendre(num_fourier+1)
    integer :: i, j, k

!----------------------------- executable code --------------------------------

!------------------          write coordinate axes         --------------------

    ! latitudes
    call check(nf90_put_var(ncFileID, latVarID, deg_lat,                      &
         start = (/ 1 /), count = (/ num_lat /) ) )

    ! longitudes
    call check(nf90_put_var(ncFileID, lonVarID, deg_lon,                      &
         start = (/ 1 /), count = (/ num_lon /) ) )

    ! sigma level
    call check(nf90_put_var(ncFileID, sigmaVarID, sigma,                      &
         start = (/ 1 /) ) )


    ! times
    call check(nf90_put_var(ncFileID, timeVarID, num_times ) )

    ! --------------             write statistics             -----------------

    call check(nf90_put_var(ncFileID, statVarID(5), u_avg(:,:,:),                                    &
             start = (/ 1, 1, 1 /), stride = (/ 1, 1, 1 /) ))



    call check(nf90_put_var(ncFileID, statVarID(13), v_avg(:,:,:),                                   &
         start = (/ 1, 1, 1 /), stride = (/ 1, 1, 1 /) ))



    call check(nf90_put_var(ncFileID, statVarID(17), w_avg(:,:,:),                                   &
        start = (/ 1, 1, 1 /), stride = (/ 1, 1, 1 /) ))


    call check(nf90_put_var(ncFileID, statVarID(33), pot_temp_avg(:,:,:),                            &
        start = (/ 1, 1, 1 /), stride = (/ 1, 1, 1 /)))


    call check(nf90_put_var(ncFileID, statVarID(61), pot_temp_e_avg(:,:,:),                       &
        start = (/ 1, 1, 1 /), stride = (/ 1, 1, 1 /) ))

     call check(nf90_put_var(ncFileID, statVarID(22), temp_avg(:,:,:),                            &
        start = (/ 1, 1, 1 /), stride = (/ 1, 1, 1 /)))


    call check(nf90_put_var(ncFileID, statVarID(28), z_avg(:,:,:),                       &
        start = (/ 1, 1, 1 /), stride = (/ 1, 1, 1 /) ))


    call check(nf90_sync(ncFileID))

  end subroutine write_variables

! #############################################################################

  subroutine close_output_file

    call check(nf90_close(ncFileID))

  end subroutine close_output_file

! #############################################################################

  subroutine check(status)
    
    ! checks error status after each netcdf, prints out text message each time
    !   an error code is returned. 
    
    integer, intent(in) :: status
    
    if(status /= nf90_noerr) then 
       write(*, *) trim(nf90_strerror(status))
    end if
  end subroutine check
  
end module output_mod









         
         










   





