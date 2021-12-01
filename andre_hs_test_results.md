# Andre test project 
- to compare FMS to the DG ClimateMachine 

## Setup:
- Replicate Andre's setup 
    - T63: 96x192 (~1.89 deg)
    - 30 levs
    - simtime: 1000 days
    - dry

- Other choices:
    - timestep: 900s bsed on [jablonowski williamson 06](https://www.semanticscholar.org/paper/A-baroclinic-wave-test-case-for-dynamical-cores-of-Jablonowski-Williamson/6801fad37b45c5f593c11f7aff2a1dcbc111c4d1/figure/3)
    - initial conditions - from rest (after enough time should yield same results)
    - boundary conditions: default is to use surface fluxes with u,v_sfc = 0, and t_surf

- dignostics
    - online (diag_table):
        - variables
        - output frequency
    - offline (analysis_3d) - split up the scripts


- helpful tips:
    - list_exp_sd: hs_control num_fourier 63 num_num_spherical 64 lon_max 192 lat_max 96 num_levels 30
    - 

## Tests
- control run
- 

## test log
- `hs_control_num_fourier_85_num_num_spherical_86_lon_max_256_lat_max_128_num_levels_30`
    - 
    - /groups/esm/lenka/GFDL_idealized_GCM_Caltech/exp/hs_control/run/out_err/slurm_17148050.err
    - Wall time 5h50
    - sim time 90 days
    - ntasks=1
    - tstep = 300s
    - turb                 = .true.,
- `hs_control_num_fourier_63_num_num_spherical_64_lon_max_192_lat_max_96_num_levels_30`
    - /groups/esm/lenka/GFDL_idealized_GCM_Caltech/exp/hs_control/run/out_err/slurm_18210859.out
    - Wall time 20h
    - sim time 600 days
    - ntasks=1
    - tstep = 300s
    - turb                 = .true.,
- `hs_control_num_fourier_63_num_num_spherical_64_lon_max_192_lat_max_96_num_levels_30`
    - /groups/esm/lenka/GFDL_idealized_GCM_Caltech/exp/hs_control/run/out_err/slurm_18221956.out
    - Wall time 87 mins = 1.45h
    - sim time 1000 days
    - ntasks=32
    - tstep = 300s
    - turb                 = .true.,
- `hs_control_num_fourier_63_num_num_spherical_64_lon_max_192_lat_max_96_num_levels_30_include_surface_flux_.false.`
    - /groups/esm/lenka/GFDL_idealized_GCM_Caltech/exp/hs_control/run/out_err/slurm_18229487.out
    - Wall time 85 mins = 1.41h
    - sim time 1000 days
    - ntasks=32
    - tstep = 300s
    - turb   
`hs_control_num_fourier_63_num_num_spherical_64_lon_max_192_lat_max_96_num_levels_30_include_surface_flux_.false._turb_.false.`
    - /groups/esm/lenka/GFDL_idealized_GCM_Caltech/exp/hs_control/run/out_err/slurm_18269206.out
    - Wall time 73 mins = 1.21h
    - sim time 1000 days
    - ntasks=32
    - tstep = 300s
    - turb = .false.  
- `hs_control_num_fourier_63_num_num_spherical_64_lon_max_192_lat_max_96_num_levels_30_include_surface_flux_.false._turb_.false._damping_coeff_0.0`
    - /groups/esm/lenka/GFDL_idealized_GCM_Caltech/exp/hs_control/run/out_err/slurm_18269863.out
    - Wall time 74 mins = 1.24h
    - sim time 1000 days
    - ntasks=32
    - tstep = 300s
    - turb = .false.  
- `hs_control_num_fourier_63_num_num_spherical_64_lon_max_192_lat_max_96_num_levels_30_include_surface_flux_.false._turb_.false._damping_coeff_0.0_dt_atmos_900`
    - /groups/esm/lenka/GFDL_idealized_GCM_Caltech/exp/hs_control/run/out_err/slurm_18270715.out
    - Wall time 24 mins
    - sim time 1000 days
    - ntasks=32
    - tstep = 900s
    - turb = .false.  
`hs_control_num_fourier_63_num_num_spherical_64_lon_max_192_lat_max_96_num_levels_30_include_surface_flux_.false._turb_.false._damping_coeff_0.0_dt_atmos_1800`
    - /groups/esm/lenka/GFDL_idealized_GCM_Caltech/exp/hs_control/run/out_err/slurm_18278526.out
    - Wall time 12 mins
    - sim time 1000 days
    - ntasks=32
    - tstep = 1800s
    - turb = .false.  
`hs_control_num_fourier_63_num_num_spherical_64_lon_max_192_lat_max_96_num_levels_30_include_surface_flux_.false._turb_.false._damping_coeff_0.0_dt_atmos_3600`
    -/groups/esm/lenka/GFDL_idealized_GCM_Caltech/exp/hs_control/run/out_err/slurm_18322614.out
    - crashed at d 57
`hs_control_num_fourier_63_num_num_spherical_64_lon_max_192_lat_max_96_num_levels_30_include_surface_flux_.false._turb_.false._damping_coeff_0.0_dt_atmos_2700`
    -/groups/esm/lenka/GFDL_idealized_GCM_Caltech/exp/hs_control/run/out_err/slurm_18324237.out
    - crashed at d 64
`hs_control_num_fourier_63_num_num_spherical_64_lon_max_192_lat_max_96_num_levels_30_include_surface_flux_.false._turb_.false._damping_coeff_0.0_dt_atmos_2400`
    -/groups/esm/lenka/GFDL_idealized_GCM_Caltech/exp/hs_control/run/out_err/slurm_18338302.out
    - crashed at d 70

# checking subroutines
- atmosphere.F90
    - INIT
        - nml overwrites
        &atmosphere_nml      
                **turb                 = .true.,**
                dry_model            = .true.,
                hs                   = .true.,
                **roughness_mom        = 0.05,**
                roughness_heat       = 0.05,
                roughness_moist      = 0.05,
                bucket               = .false./

        - spectral_dynamics_init
            - default_representation = 'spectral'
            - vert_advect_t = vert_advect_uv = default_advect_vert    = 'second_centered'
            - vert_difference_option = 'simmons_and_burridge'
            - vert_coord_option      = 'even_sigma'
            - damping_option         = 'resolution_dependent'
            - **damping_coeff       = 1.15740741e-4,** --> leave for now
            - all sponge coeffs = 0
            - default_hole_filling   = 'off'
            - **robert_coeff        = .04,**  # robert filter --> leave (default)
            - **do_energy_correction = do_mass_correction = .true.** --> leave for now
            - use_implicit = true
            
        - gcm_vert_diff_init
            - not doing much here

        - vert_turb_driver_init
    - MAIN  
        - **surface_flux**
                """
                real, intent(in),  dimension(:) :: &
                    t_atm,     q_atm_in,   u_atm,     v_atm,              &
                    p_atm,     z_atm,      t_ca,                          &
                    p_surf,    t_surf,     u_surf,    v_surf,  &
                    rough_mom, rough_heat, rough_moist,  gust
                real, intent(out), dimension(:) :: &
                    flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
                    dhdt_surf, dedt_surf,  dedq_surf, drdt_surf,          &
                    dhdt_atm,  dedq_atm,   dtaudv_atm,                    &
                    w_atm,     u_star,     b_star,    q_star,             &
                    cd_m,      cd_t,       cd_q
                """
        - **vert_turb_driver** --> switch off with turb = .false.
                """
                real, intent(in), dimension(:,:) :: frac_land,   &
                                                        u_star, b_star, rough
                real, intent(in), dimension(:,:,:) :: p_half, p_full, &
                                                        z_half, z_full, &
                                                        u, v, t, q, um, vm, tm, qm, &
                                                        udt, vdt, tdt, qdt
                real, intent(out),   dimension(:,:,:) :: diff_t, diff_m
                real, intent(out),   dimension(:,:)   :: gust
                """
                - do_molecular_diffusion = .false.
            - diffusivity
                - diffusivity_pbl 
                        """
                            real,    intent(in)  ,           dimension(:,:,:) :: t, u, v, z_half
                            real,    intent(in)  ,           dimension(:,:)   :: h, u_star, b_star
                            real,    intent(inout) ,         dimension(:,:,:) :: k_m, k_t
                        """
                    - mo_diff(..)
                - diffusivity_free = .false.
        - spectral_dynamics
        - spectral_diagnostics
        - **alpha_implicit** = 0.5 : centered implicit gravity wave scheme --> ok to leave - default

Refs
- idealized/src/atmos_spectral/documentation/spectral_core.pdf



Notes
- since triangular truncation: num_total_wavenumbers = num_spherical - 1

Q
- what is the implicit_correction?
    - adjusting p and T using linear_geopotential_1d

    vert_diff = implicit algorithm


    If the dependence of the surface heat flux and evaporation
       on surface temperature is treated implicitly, then interfaces 
       are provided for dividing the up and down sweeps of the 
       standard tridiagonal elimination into two parts.  One can 
       visualize this procedure as starting with a downward 
       sweep through the atmosphere, which continue through the 
       surface, then re-emerges and continues up the atmosphere.
       A sample program is included below to help the user undertand
       the approriate calling sequence

COupler:
The component models are coupled to allow implicit vertical
     diffusion of heat and moisture at the interfaces of the atmosphere,
     land, and ice models. As a result, the atmosphere, land, and ice
     models all use the same time step.  The atmospheric model has been
     separated into down and up calls that correspond to the down and up
     sweeps of the standard tridiagonal elimination.

      <a href="#flux_down_from_atmos">flux_down_from_atmos:</a>
 <dd>   Returns fluxes and derivatives corrected for the implicit
        treatment of atmospheric diffusive fluxes, as well as the
        increments in the temperature and specific humidity of the
        lowest atmospheric layer due to all explicit processes as well
        as the diffusive fluxes through the top of this layer.