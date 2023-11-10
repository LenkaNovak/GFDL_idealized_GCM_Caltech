# exec(open("test_data_cospectra.py").read())

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
exec(open("../helpers/data_helper.py").read())
exec(open("../helpers/cospectra_helper.py").read())

day_list = np.arange(1,361)
last_day_of_yr = "1440"# "%04d"%day_list[-1]

# control run
exp_name, run_name = ("uf_tests_hpc_albedo", "uf_tests_hpc_albedo_albedo02_1")
work_dir = "/central/scratch/elencz/caltech_fms_idealized/%s/%s/workdir"%(exp_name, run_name)
output_dir = "/central/scratch/elencz/caltech_fms_idealized/%s/%s/output"%(exp_name, run_name)
ppp_dir = output_dir + "/ppp_output/"
fname = "day%sh00.segment%s.nc"
lat_c, lon_c, sig_c, v_c, u_c, T_c, Th_c, uv_c, vT_c, TT_c, eke_c, u_vrtcl_flux_c, pot_temp_vrtcl_flux_c, z1_c, T1_c, u1_c, lh_flux_sfc_c, sh_flux_sfc_c, T_sfc_c, u_02_c, v_02_c = append_year_days_zonal_mean(ppp_dir+fname%(last_day_of_yr, "%s"), day_list)

# control run
exp_name, run_name = ("uf_tests_hpc_businger_albedo", "uf_tests_hpc_businger_albedo_albedo02rcrit15_1")
work_dir = "/central/scratch/elencz/caltech_fms_idealized/%s/%s/workdir"%(exp_name, run_name)
output_dir = "/central/scratch/elencz/caltech_fms_idealized/%s/%s/output"%(exp_name, run_name)
ppp_dir = output_dir + "/ppp_output/"
fname = "day%sh00.segment%s.nc"
lat, lon, sig, v, u, T, Th, uv, vT, TT, eke, u_vrtcl_flux, pot_temp_vrtcl_flux, z1, T1, u1, lh_flux_sfc, sh_flux_sfc, T_sfc, u_02, v_02 = append_year_days_zonal_mean(ppp_dir+fname%(last_day_of_yr, "%s"), day_list)

def get_phase_speed_spectra_v_latitude(u_02, v_02, day_list, lat, lon):
    t = np.arange(0, len(day_list))
    lat_sampled = lat[::15]
    for l4 in lat_sampled:
        la_i = np.where(lat ==l4)[0][0]
        v_lev = v_02[:,la_i,:] # t, lon
        u_lev = u_02[:,la_i,:]
        calculate_cospectra(v_lev, v_lev, lat, lon, t, dt_secs = 86400.0, fdir=".", fname="cospectra_lat_%00d_%s"%(l4, last_day_of_yr))
    nc = 60
    nm = int(len(lon) / 2)
    # get phase speed spectra
    K_wm = np.zeros((len(lat_sampled), int(2*nc), nm ))
    for i, l4 in enumerate(lat_sampled):
        K_wm[i, :, :], cp, m = phase_speed_spectra(l4, nc = nc, cmax = 40, fdir = ".", fname = "cospectra_lat_%00d_%s"%(l4, last_day_of_yr))
    return K_wm, cp, m, lat_sampled


K_wm, cp, m, lat_sampled = get_phase_speed_spectra_v_latitude(u_02, v_02, day_list, lat, lon)

K_wm_c, cp, m, lat_sampled = get_phase_speed_spectra_v_latitude(u_02_c, v_02_c, day_list, lat, lon)

cp_2sided = np.append(-np.flipud(cp),cp)

plt.figure()
plt.subplot(121)
c = plt.contourf(lat_sampled,cp_2sided,np.transpose(np.sum(K_wm,axis=-1)))
plt.plot(lat, np.mean(u_02, axis=(0,-1)))
plt.subplot(122)
c = plt.contourf(lat_sampled,cp_2sided,np.transpose(np.sum(K_wm_c,axis=-1)))
plt.plot(lat_c, np.mean(u_02_c, axis=(0,-1)))
plt.savefig("ps_spectrum.png")

plt.figure()
plt.subplot(121)
c = plt.contourf(lat_sampled,cp_2sided,np.transpose(np.sum(K_wm,axis=-1))- np.transpose(np.sum(K_wm_c,axis=-1)))
plt.plot(lat_c, np.mean(u_02_c, axis=(0,-1)))
plt.colorbar(c)
plt.savefig("ps_spectrum_diff.png")


    # c = plt.contourf(w,m,np.transpose(np.mean(K_wm,axis=0)), rrnge)
    # if cbar ==True:
    #     plt.colorbar(c)



