# Numerical details

- [doc](/groups/esm/lenka/GFDL_idealized_GCM_Caltech/idealized/src/atmos_spectral/documentation/spectral_core.pdf) 


# Moisture 
- ignore in the dynamical core (obs mixing rations < 2%), except the effect of vapour on density $T_v$ where the mixing ratio is $m_w / m_d = R_d / R_w \approx 0.622$ so 
$$
p = p_d + p_w = (\rho R_d + \rho_w(R_w - R_d))T = \rho R_d (1 + q\underbrace{(R_w/R_d - 1)}_{\mu = 0.608} T = \rho R_d T_v
$$
where the specific humidity is $q = \rho_w/\rho$. 


# Prognostic Equations

## TD:
$$
\frsc{DT}{Dt} = \frac{RT_v \omega}{p c_p} + \frac{Q}{c+p}
$$

## Motion:

In spherical coordinates:
$$
\frsc{Du}{Dt} = 
$$

$$
\frsc{Dv}{Dt} = 
$$

$$
\frsc{Dw}{Dt} = 
$$

- hydrostatic ok (as long as aspect ratios not of order unity) and assume radius is constant = primitive equarions or traditional approximation
    - prognostic vars: $u,v,T,\rho$ 

- convert to pressure coordinstes:
    - 
    $$
    \frac{1}{\rho}  \nabla \cdot p \rightarrow \nabla \cdot \Phi  
    $$

$$
\frac{\partial p_s}{\partial t} = -\nabla \cdot \int_0^{p_s} \vec{v} dp
$$
so the lower BC is no longer a coordinate surface. In  a  spectral  model,  in  particular,  thereis no alternative to using a terrain following coordinate.  The simplest and traditional choice is Phillips’ sigma-coordinate,σ≡p/p. 

# Discretisation
- Simmons, A. J. and D. M. Burridge, 1981: 
- hybrid coord sigma and pressure with different weights, see doc above. 

# Timestepping
- The code uses a standardsemi-implicitleapfrog scheme, followed by a Robert-Asselin time filter
- Use a centered leapfrog scheme for explicit part of the operator, N_E, and an implicit scheme for the implicit part, N_I
$$
\frac{V^{i+1}-V^{i-1}}{2\Delta t} = N_E(V^i) + N_I(\alpha V^{i+1} + (1 - \alpha) V^{i-1})
$$
- assume that N_I os linear (otherwise expensive)
- \$\alpha =0.5$ = centered implicit; but one can also consider a fullybackward scheme (α= 1.0), for example, if one prefers to damp as well asslow down the gravity waves.
- We need first of all to define the parts of the equations that are to be partof the implicit operatorNI.  This is determined by linearizing the equationsfor a non-rotating atmosphere about a state of rest, in which orography isignored, surface pressure is uniform, and temperatures are a function of thevertical coordinate only
- Robert-Asselin filter - suppresses the computational model of the leapfrog scheme
