module rayleigh_drag_forcing

! Rayleigh drag 

use     constants_mod, only: kappa, CP_air, RADIUS

use           fms_mod, only: error_mesg, fatal, file_exist,       &
                             open_namelist_file, check_nml_error, &
                             mpp_pe, mpp_root_pe, close_file,     &
                             write_version_number, stdlog

use  time_manager_mod, only: time_type
use  diag_manager_mod, only: register_diag_field, send_data
!
  implicit none
  
  real, private ::                     & ! 
       kf_days             =     1.,   & ! Rayleigh drag timescale
       kf,                             & ! 
       day                 = 86400. 

  real, public ::                      &
       sigma_b             = 0.85,     & ! extent of frictional `PBL'
       rc                  = 0.84,     & ! radius of the magnetic core
       H_lambda            = 100.0e3     ! scale height of the magnetic diffusivity      

  logical ::                           &
       do_energy_conserv_ray = .true., & 
             ! convert the energy produced by friction to heat
       variable_drag        = .true. 
             ! consider latitudinally varying drag produced by magnetic field
             
  namelist/rayleigh_drag_list/ sigma_b, rc, H_lambda,             &
              kf_days, do_energy_conserv_ray, variable_drag

  private rayleigh_drag_list

  character(len=128) :: version='$Id: rayleigh_drag_forcing.f90 $'
  character(len=128) :: tag='homemade'
  character(len=15)  :: mod_name = "rayleigh_drag_forcing"

! diagnostic IDs
  integer id_udt, id_vdt
!  integer  id_dissdt
  real :: missing_value = -1.e10
contains
  
  !-------------------------------------------------------------------------------

  subroutine rayleigh_drag_forcing_init(axes, Time)
    ! reads namelist file with parameters and echoes parameter values

    integer, intent(in) :: axes(4)
    type(time_type), intent(in) :: Time

    integer  unit, io, ierr

      if (file_exist('input.nml')) then
         unit = open_namelist_file ( )
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=rayleigh_drag_list, iostat=io, end=100)
            ierr = check_nml_error (io, 'rayleigh_drag_list')
         enddo
  100     call close_file (unit)
      endif

!     ----- write version info and namelist to log file -----

      call write_version_number (version,tag)
      if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=rayleigh_drag_list)

!     ---------------write info to screen--------------------

      if(mpp_pe()==mpp_root_pe()) then
         write(*, *)
         write(*, *)
         write(*, *) 'Rayleigh Drag in boundary layer.'
         write(*,10) 'Extent of frictional layer: sigma_b =', sigma_b
         write(*,10) 'Frictional timescale in days: ', kf_days
         write(*, *)
         write(*, *)
      endif

      kf = 1./(kf_days * day)
      
    
10  format(1x,a,1x,f8.2,1x,a)
20  format(1x,a,1x,e8.2,1x,a)


! register fields with diagnostic manager
      id_udt = register_diag_field( mod_name, 'udt_rd',  &
           axes(1:3), Time, 'zonal wind tendency', 'm/sec/sec',             &
           missing_value=missing_value)

      id_vdt = register_diag_field( mod_name, 'vdt_rd',  &
           axes(1:3), Time, 'meridional wind tendency', 'm/sec/sec',        &
           missing_value=missing_value)

  end subroutine rayleigh_drag_forcing_init
  
  !------------------------------------------------------------------

  subroutine compute_rayleigh_drag_forcing( is,       ie,     js,   je,    &
                   Time, delta_t,   lat,  dt_ug,    dt_vg,    ug,   vg,    &
                   p_half,   p_full,   ug_previous,      vg_previous,      &
                   p_half_previous,  p_full_previous,                      &
                   dt_tg, dissipative_heat)

    integer, intent(in)                     :: is, ie, js, je
    type(time_type), intent(in)            ::  Time
    real, intent(in)                       ::  delta_t
    real, intent(in), dimension (:,:,:)     :: ug, vg, p_half, p_full
    real, intent(in), dimension (:,:)       :: lat
    
    real, optional, dimension (:,:,:)       :: ug_previous, vg_previous, & 
                                        p_half_previous, p_full_previous
    
    real, intent(inout), dimension(:,:,:)   :: dt_ug, dt_vg, dt_tg
    real, intent(out),   dimension(:,:,:)   :: dissipative_heat


    if(present(ug_previous)) then
          call surface_drag(Time, delta_t, lat, is, js, ug_previous, vg_previous,    &
                            p_half_previous,  p_full_previous, dt_ug, dt_vg,    &
                            dt_tg, dissipative_heat)
    else
       call error_mesg('compute_rayleigh_drag_forcing',                         &
            'ug_previous is not present',fatal)
    end if

  end subroutine compute_rayleigh_drag_forcing

  !---------------------------------------------------------------------
  
  subroutine surface_drag(Time, delta_t, lat, is, js, ug, vg, p_half, p_full, & 
                         dt_ug, dt_vg, dt_tg, dissipative_heat )

    real,    intent(in)                       :: delta_t
    real,    intent(in),    dimension(:,:)    :: lat
    integer, intent(in)                       :: is, js
    real,    intent(in),    dimension (:,:,:) :: ug, vg, p_full, p_half
    real,    intent(inout), dimension(:,:,:)  :: dt_ug, dt_vg, dt_tg
    real,    dimension(size(ug,1),size(ug,2)) :: sigma, sigma_norm, sigma_max
    real,    intent(out),   dimension(:,:,:)  :: dissipative_heat 
    
    real, dimension(size(ug,1),size(ug,2),size(ug,3)):: dt_u_temp, dt_v_temp
    real, dimension(size(ug,2))  :: drag_coeff
    integer :: j, k, num_level 
    type(time_type), intent(in)            ::  Time
    logical :: used
    real :: half_delt, cp_inv

! ------------------------------------------------------------------------------
   cp_inv = 1.0/CP_air

   if (do_energy_conserv_ray) then
      dt_u_temp = dt_ug
      dt_v_temp = dt_vg
   endif

  ! pseudo=parameter for number of levels
    num_level = size(ug,3)

   if (variable_drag) then
      do j = 1, size(ug,2)
         if (cos(lat(1,j)) .le. rc) then
            drag_coeff(j) = kf
         else
           drag_coeff(j)  = kf * exp(-(cos(lat(1,j)) - rc)*RADIUS/H_lambda)
         endif
      end do
   else
      drag_coeff = kf
   endif


    do k = 1, num_level
       do j = 1, size(ug,2)
       sigma(:, j)  = p_full(:,j,k) / p_half(:,j,num_level+1)
       sigma_norm(:,j)   = (sigma(:,j) - sigma_b) / (1.0 - sigma_b)
       sigma_max(:,j)    = max(sigma_norm(:,j), 0.0)
       dt_ug(:,j,k) = dt_ug(:,j,k)     &
                        - drag_coeff(j) * sigma_max(:,j) * ug(:,j,k)
       dt_vg(:,j,k) = dt_vg(:,j,k)     &
                        - drag_coeff(j) * sigma_max(:,j) * vg(:,j,k)
       end do
    end do

    if (do_energy_conserv_ray) then
       dt_u_temp = dt_ug - dt_u_temp
       dt_v_temp = dt_vg - dt_v_temp
       dissipative_heat = - cp_inv * ( (ug + 0.5 * delta_t * dt_u_temp)*dt_u_temp   &
                                      +(vg + 0.5 * delta_t * dt_v_temp)*dt_v_temp )
       dt_tg = dt_tg + dissipative_heat
    else
       dissipative_heat = 0.0
     endif


    ! send data to diagnostic manager
    if(id_udt > 0) used = send_data(id_udt, dt_ug, Time,      &
                                    is, js)
    if(id_vdt > 0) used = send_data(id_vdt, dt_vg, Time,      &
                                    is, js)
!    if(id_dissdt > 0) used = send_data(id_dissdt, dissipative_heat, &
!                                    Time, is, js)    


  end subroutine surface_drag

 !-----------------------------------------------------------------------


  function surface_temperature(temp, p_full, p_half)

    ! estimate the surface temperature from the temperature on lowest
    ! full level by assuming that the potential temperature on the surface
    ! is equal to the potential temperature on the lowest full level

    implicit none

    real, intent(in) ::                                                    &
         temp(:, :, :),    & ! temperature field
         p_full(:, :, :),  & ! full-level pressure
         p_half(:, :, :)     ! half-level pressure

    real, dimension(size(temp, 1), size(temp, 2)) ::                       &
         surface_temperature

    integer :: num_level
    num_level = size(temp, 3)

    surface_temperature = temp(:, :, num_level)                            &
         * ( p_half(:, :, num_level+1) / p_full(:, :, num_level) )**kappa

  end function surface_temperature

end module rayleigh_drag_forcing
