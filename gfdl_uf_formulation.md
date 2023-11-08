# Universal Functions in GFDL

## General intro

# MOST Formulation
We aim to inform the GCM with unresolved vertical turbulent fluxes in the surface layer (SL, which forms the bottom ~ 10% of the planetary boundary layer). For momentum fluxes (or drag), the formulation is:

$$
\frac{\partial u}{\partial t} = ... - \frac{\partial \overline{w'u'}}{\partial z}
$$

and an analogous form is applicable for tracers.
- Note: the layer is assumed to be close to horizontally homogeneous, so only vertical variations (and fluxes) are of our primary interest. Vertical fluxes themselves are assumed constant in the SL, and the u and T profiles are assumed to be functions of surface stress, thermal fluxes and height only.

The SL fluxes can be modelled using characteristic scales:

$$
\overline{w'u'} \approx u^{*2} \approx - K_m \frac{\partial \overline{u}}{\partial z} \tag{1}
$$

$u^{*}$ is also referred to as friction velocity, and $K_m$ is the diffusivity.


If the layer neutrally buoyant, the diffusivity takes the following form:

$$
K_m = (\kappa z)^2 \frac{\partial \overline{u}}{\partial z}
$$

so when (1) is square rooted and integrated, we recover a logarithmic $\overline{u}(z)$ profile.

$$
u(z) = \frac{u_{\star}}{\kappa}\log (\frac{z}{z_{0m}})
$$

with $\kappa$ representing the dimensionless von-Karman constant that captures the slope of near-wall velocity (~0.4, obtained from the law-of-the-wall in boundary layer theory), and $z_{0m}$ the roughness length (whose value depends on the height of the obstacles on the surface. E.g. $z_{0m}=0.01$ m for grassland)
This neutral profile has been verified and justified extensively (following its inception as the "law of the wall" by von Karman in the 1930s) using observations and DNS, so it provides a solid basis (and a constraint) for our parameterization.

If the atmosphere is statically stable/unstable, the neutral log profile needs to be modified using a stability function:

$$
K_m = (kz)^2 \phi_m (z/L_{mo})^2 \frac{\partial \overline{u}}{\partial z}
$$

- where $\phi_m$ is called a stability function, with a so far unknown form, and $L_{mo}$ is the Monin Obukhov lengthscale which is a ratio of the friction and buoyancy characteristic scales:

$$
L_{mo} = \frac{u_{\star}^2}{\kappa b_{\star}},
$$


For GCM parameterizations it is again useful to integrate this expression to obtain a $u(z)$ profile. This enables the calculation of the surface fluxes from the states of the surface and the overlying atmosphere (e.g., the height of the first GCM model level). Our modified expression for the near-wall velocity profile is:

$$
\frac{\partial u}{\partial z} = \frac{u_{\star}}{\kappa z} \phi_{m}(\zeta),
$$

Collapsing $\zeta = \frac{z}{L_{mo}}$ into a dimensionless stability parameter, and defining integrated *stability correction functions* $\psi$ such that

$$
\psi_{m}(\zeta) = \int_{0}^{\zeta}{\frac{1- \phi_{m}(\zeta^{\prime})}{\zeta^{\prime}} d\zeta^{\prime}},
$$

we can express the velocity profile in a planetary boundary layer, using Taylor expansion around $\zeta = 0$, as

$$
u(z) = \frac{u_{\star}}{\kappa}\Big[\log (\frac{z}{z_{0m}}) - \psi_{m}(\frac{z}{L_{mo}}) + \psi_{m}(\frac{z_{0m}}{L_{mo}}) \Big]
$$

## Ref
- Garratt, J. R. "The Atmospheric Boundary Layer", Cambridge University Press, 1992.
-

##

### init:

b_stab = 1.0/rich_crit
r_crit = 0.95*rich_crit


   delta_b = grav*(pt0 - pt)/pt0
   rich    = - z*delta_b/(speed*speed + small)
   zz      = max(z,z0,zt,zq)



### mo_drag

computes the drag coefficients for momentum, heat, and moisture;
                also returns u_star (the friction velocity) and
                b_star (the buoyancy scale)


   delta_b = grav*(pt0 - pt)/pt0 ! buoyancy
   rich    = - z*delta_b/(speed*speed + small) ! richardson number
   zz      = max(z,z0,zt,zq)


Neutral
us/ speed = kappa * log(zz/z0)

Stable
us/ speed = kappa * log(zz/z0) - psi_m(z/Lmo) + psi_m(z0/Lmo)


mask_2 = mask .and. rich >= r_crit
-> min values for ustar and bstar

mask_1 = mask .and. rich <  r_crit
-> solve_zeta (OUT:  f_m, f_t, f_q)
    us   = max(vonkarm/fm, sqrt_drag_min)
    bs   = max(vonkarm/ft, sqrt_drag_min)
    qs   = max(vonkarm/fq, sqrt_drag_min)
    drag_m   = us*us
    drag_t   = us*bs
    drag_q   = us*qs
    u_star   = us*speed
    b_star   = bs*delta_b

    strt iter
        mo_derivative_m(OUT: phi_m)
            where (unstable)
            x     = (1 - 16.0*zeta  )**(-0.5)
            phi_m = sqrt(x)  ! phi_m = (1 - 16.0*zeta)**(-0.25)
            end where

            if(stable_option == 1) then

            where (stable)
                phi_m = 1.0 + zeta  *(5.0 + b_stab*zeta)/(1.0 + zeta)
            end where

        mo_integral_m(OUT: psi_m)
            where(unstable)

            x     = sqrt(1 - 16.0*zeta)
            x_0   = sqrt(1 - 16.0*zeta_0)

            x      = sqrt(x)
            x_0    = sqrt(x_0)

            x1     = 1.0 + x
            x1_0   = 1.0 + x_0

            num    = x1*x1*(1.0 + x*x)
            denom  = x1_0*x1_0*(1.0 + x_0*x_0)
            y      = atan(x) - atan(x_0)
            psi_m  = ln_z_z0 - log(num/denom) + 2*y

            end where

            if( stable_option == 1) then

            where (stable)
                psi_m = ln_z_z0 + (5.0 - b_stab)*log((1.0 + zeta)/(1.0 + zeta_0)) &
                + b_stab*(zeta - zeta_0)
            end where

        assess error

            zeta_0 = zeta/z_z0
            zeta_t = zeta/z_zt
            zeta_q = zeta/z_zq

            df_m  = (phi_m - phi_m_0)*rzeta
            df_t  = (phi_t - phi_t_0)*rzeta
            rich_1 = zeta*f_t/(f_m*f_m)
            d_rich = rich_1*( rzeta +  df_t/f_t - 2.0 *df_m/f_m)
            correction = (rich - rich_1)/d_rich
            corr = min(abs(correction),abs(correction/zeta))
            ! the criterion corr < error seems to work ok, but is a bit arbitrary
            !  when zeta is small the tolerance is reduced

        max_cor= maxval(corr)


        if(max_cor > error) then
            mask_1 = mask_1 .and. (corr > error)
            ! change the mask so computation proceeds only on non-converged points
            where(mask_1)
            zeta = zeta + correction
            end where
            cycle iter_loop



Byun 1990 https://journals.ametsoc.org/view/journals/apme/29/7/1520-0450_1990_029_0652_otasof_2_0_co_2.xml?tab_body=pdf


