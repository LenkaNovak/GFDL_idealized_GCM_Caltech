# exec(open("test_data_albedo_comparison.py").read())

#modules
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import integrate
import h5py as h5
from scipy import ndimage

# global parameters:
_a = 6371000.0
_g = 9.81

# read in data_helper.py
exec(open("data_helper.py").read())

day_list = np.arange(1,361)
last_day_of_yr = "1440"# "%04d"%day_list[-1]

# control run
exp_name, run_name = ("uf_tests_hpc_albedo", "uf_tests_hpc_albedo_albedo02_1")
work_dir = "/central/scratch/elencz/caltech_fms_idealized/%s/%s/workdir"%(exp_name, run_name)
output_dir = "/central/scratch/elencz/caltech_fms_idealized/%s/%s/output"%(exp_name, run_name)
ppp_dir = output_dir + "/ppp_output/"
fname = "day%sh00.segment%s.nc"
lat_c, lon_c, sig_c, v_c, u_c, T_c, Th_c, uv_c, vT_c, TT_c, eke_c, u_vrtcl_flux_c, pot_temp_vrtcl_flux_c, z1_c, T1_c, u1_c, lh_flux_sfc_c, sh_flux_sfc_c, T_sfc_c = append_year_days_zonal_mean(ppp_dir+fname%(last_day_of_yr, "%s"), day_list)

# control run
exp_name, run_name = ("uf_tests_hpc_businger_albedo", "uf_tests_hpc_businger_albedo_albedo02rcrit15_1")
work_dir = "/central/scratch/elencz/caltech_fms_idealized/%s/%s/workdir"%(exp_name, run_name)
output_dir = "/central/scratch/elencz/caltech_fms_idealized/%s/%s/output"%(exp_name, run_name)
ppp_dir = output_dir + "/ppp_output/"
fname = "day%sh00.segment%s.nc"
lat, lon, sig, v, u, T, Th, uv, vT, TT, eke, u_vrtcl_flux, pot_temp_vrtcl_flux, z1, T1, u1, lh_flux_sfc, sh_flux_sfc, T_sfc = append_year_days_zonal_mean(ppp_dir+fname%(last_day_of_yr, "%s"), day_list)


# plotting helpers
def zonal_plot(lat, lev, var, title = '', subplot_no=(1,1,1), vrange = 20):
    a,b,c = subplot_no
    plt.subplot(a,b,c)
    f1 = plt.contourf(lat, lev, var, vrange)
    plt.colorbar(f1)
    plt.gca().invert_yaxis()
    plt.title(title)

# diagnostic helpers
def psi(V_xBAR,level_list_array,lat):
    _a, _g = (6371000. , 9.81) # Earth_radius, grav_constant
    latrad = lat*np.pi/180.
    psi = _a/_g * 2. * np.pi * np.cos(latrad)[None,None,:] * integrate.cumtrapz(V_xBAR,level_list_array,axis = 1)
    return psi

# U, T, meridional streamfunctions
plt.figure(figsize=(20,4))
zonal_plot(lat, sig, np.nanmean(u,axis=0), title = 'U', subplot_no = (1,4,1), vrange = np.arange(-32,33,4))
zonal_plot(lat, sig, np.nanmean(T,axis=0), title = 'T', subplot_no = (1,4,2), vrange = np.arange(170,330,10))
zonal_plot(lat, sig, np.nanmean(TT,axis=0), title = 'TT', subplot_no = (1,4,3), vrange = np.arange(-32,41,8))
pressure = sig * 100000.0 # units: Pa
zonal_plot(lat, sig[:-1], np.nanmean(psi(v, pressure ,lat),axis=0), title = 'meridional streamfunction', subplot_no = (1,4,4), vrange = np.arange(-8.8e10,8.9e10,0.8e10))
plt.savefig("circ1.png")

plt.figure(figsize=(20,4))
zonal_plot(lat, sig, np.nanmean(u,axis=0) -  np.nanmean(u_c,axis=0), title = 'U', subplot_no = (1,4,1))
zonal_plot(lat, sig, np.nanmean(T,axis=0) -  np.nanmean(T_c,axis=0), title = 'T', subplot_no = (1,4,2))
zonal_plot(lat, sig, np.nanmean(TT,axis=0) -  np.nanmean(TT_c,axis=0), title = 'TT', subplot_no = (1,4,3))
pressure = sig * 100000.0 # units: Pa
zonal_plot(lat, sig[:-1], np.nanmean(psi(v, pressure ,lat),axis=0) - np.nanmean(psi(v_c, pressure ,lat),axis=0), title = 'meridional streamfunction', subplot_no = (1,4,4))
plt.savefig("circ1_diff.png")

# vT, uv, u^2+v^2, TT
plt.figure(figsize=(20,4))
zonal_plot(lat, sig, np.nanmean(uv,axis=0), title = 'uv', subplot_no = (1,4,1) , vrange = np.arange(-66,67,6))
zonal_plot(lat, sig, np.nanmean(vT,axis=0), title = 'vT', subplot_no = (1,4,2) , vrange = np.arange(-18,19,2))
zonal_plot(lat, sig, np.nanmean(0.5*eke,axis=0), title = 'eke', subplot_no = (1,4,3), vrange = np.arange(0,330,15))
plt.savefig("circ2.png")

# diff
plt.figure(figsize=(20,4))
zonal_plot(lat, sig, np.nanmean(uv,axis=0) - np.nanmean(uv_c,axis=0), title = 'uv', subplot_no = (1,4,1) )
zonal_plot(lat, sig, np.nanmean(vT,axis=0) - np.nanmean(vT_c,axis=0), title = 'vT', subplot_no = (1,4,2) )
zonal_plot(lat, sig, np.nanmean(0.5*eke,axis=0) - np.nanmean(0.5*eke_c,axis=0), title = 'eke', subplot_no = (1,4,3))
plt.savefig("circ2_diff.png")


# turbulent surface fluxes
plt.figure(figsize=(20,4))
zonal_plot(lon, lat, (np.nanmean(lh_flux_sfc,axis=0)), title = 'latent (W/m2)', subplot_no = (1,4,1) , vrange = np.arange(-50,51,5))
zonal_plot(lon, lat, np.nanmean(sh_flux_sfc,axis=0), title = 'sensible (W/m2)', subplot_no = (1,4,2) , vrange = np.arange(-50,51,5))
plt.savefig("test_data_albedo_comparison.png")

# differences of the above figures
plt.figure(figsize=(20,4))
zonal_plot(lon, lat, (np.nanmean(lh_flux_sfc,axis=0) - np.nanmean(lh_flux_sfc_c,axis=0)), title = 'latent (W/m2)', subplot_no = (1,4,1) )
zonal_plot(lon, lat, np.nanmean(sh_flux_sfc,axis=0) - np.nanmean(sh_flux_sfc_c,axis=0), title = 'sensible (W/m2)', subplot_no = (1,4,2) )
plt.savefig("test_data_albedo_comparison_diff.png")






