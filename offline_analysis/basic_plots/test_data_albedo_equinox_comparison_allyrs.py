# exec(open("test_data_albedo_equinox_comparison_allyrs.py").read())

#modules
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import integrate
import h5py as h5
from scipy import ndimage
import os
import dask.array as da

# global parameters:
_a = 6371000.0
_g = 9.81

# read in data_helper.py
exec(open("../helpers/data_helper.py").read())

day_list = np.arange(1,361)

yr_range = ["1440"]#, '1800'] # ["1080", "1440", '1800'] # "%04d"%day_list[-1]

save_variable_list = ['lat', 'lon', 'sig', 'v', 'u', 'T', 'Th1', 'uv', 'vT', 'TT', 'eke', 'u_vrtcl_flux', 'pot_temp_vrtcl_flux', 'z1', 'T1', 'u1', 'lh_flux_sfc', 'sh_flux_sfc', 'T_sfc', 'u_02', 'v_02']

# Save variables
def save_basic_diags(exp_name, run_name, day_list, yr_range, save_variable_list, h5_dir = "data/", h5_name = "basic_vars"):

    for last_day_of_yr in yr_range:

        work_dir = "/central/scratch/elencz/caltech_fms_idealized/%s/%s/workdir"%(exp_name, run_name)
        output_dir = "/central/scratch/elencz/caltech_fms_idealized/%s/%s/output"%(exp_name, run_name)
        ppp_dir = output_dir + "/ppp_output/"
        nc_name = "day%sh00.segment%s.nc"
        nc_path = ppp_dir+nc_name%(last_day_of_yr, "%s")
        if not os.path.exists(h5_dir):
            os.makedirs(h5_dir)

        h5_path = "%s/%s_%s.hdf5"%(h5_dir, h5_name, last_day_of_yr)

        # if h5 file already exists, skip
        if os.path.exists(h5_path):
            continue
        else:
            print("saving %s"%h5_path)

            lat, lon, sig, v, u, T, Th1, uv, vT, TT, eke, u_vrtcl_flux, pot_temp_vrtcl_flux, z1, T1, u1, lh_flux_sfc, sh_flux_sfc, T_sfc, u_02, v_02 = append_year_days_zonal_mean(nc_path, day_list)

            # save as only the tail name
            hfile = h5.File(h5_path,'w')
            for var in save_variable_list:
                hfile.create_dataset(var, data = eval(var),compression = 4)
            hfile.close()

# read in data
def read_basic_diags(last_day_of_yr, h5_dir = "data/", h5_name = "basic_vars"):

    # create empty arrays for all variables
    lat, lon, sig = [np.array([]) for i in range(3)]

    hfile = h5.File("%s/%s_%s.hdf5"%(h5_dir, h5_name, yr_range[0]),'r')
    lat = hfile['lat'][:]
    lon = hfile['lon'][:]
    sig = hfile['sig'][:]
    hfile.close()

    v, u, T, uv, vT, TT, eke = [np.zeros((360, len(sig), len(lat))) for i in range(7)]

    Th1, z1, T1, u1, lh_flux_sfc, sh_flux_sfc, T_sfc, u_02, v_02, u_vrtcl_flux, pot_temp_vrtcl_flux = [np.zeros((360, len(lat), len(lon))) for i in range(11)]

    # read in variables for all years
    if append_yrs == True:
        for yr in yr_range:
            hfile = h5.File("%s/%s_%s.hdf5"%(h5_dir, h5_name, yr),'r')
            v = np.append(v, hfile['v'][:], axis = 0)
            u = np.append(u, hfile['u'][:], axis = 0)
            T = np.append(T, hfile['T'][:], axis = 0)
            Th1 = np.append(Th1, hfile['Th1'][:], axis = 0)
            uv = np.append(uv, hfile['uv'][:], axis = 0)
            vT = np.append(vT, hfile['vT'][:], axis = 0)
            TT = np.append(TT, hfile['TT'][:], axis = 0)
            eke = np.append(eke, hfile['eke'][:], axis = 0)
            u_vrtcl_flux = np.append(u_vrtcl_flux, hfile['u_vrtcl_flux'][:], axis = 0)
            pot_temp_vrtcl_flux = np.append(pot_temp_vrtcl_flux, hfile['pot_temp_vrtcl_flux'][:], axis = 0)
            z1 = np.append(z1, hfile['z1'][:], axis = 0)
            T1 = np.append(T1, hfile['T1'][:], axis = 0)
            u1 = np.append(u1, hfile['u1'][:], axis = 0)
            lh_flux_sfc = np.append(lh_flux_sfc, hfile['lh_flux_sfc'][:], axis = 0)
            sh_flux_sfc = np.append(sh_flux_sfc, hfile['sh_flux_sfc'][:], axis = 0)
            T_sfc = np.append(T_sfc, hfile['T_sfc'][:], axis = 0)
            u_02 = np.append(u_02, hfile['u_02'][:], axis = 0)
            v_02 = np.append(v_02, hfile['v_02'][:], axis = 0)
            hfile.close()
        else:
            hfile = h5.File("%s/%s_%s.hdf5"%(h5_dir, h5_name, yr_range[0]),'r')
            v += hfile['v'][:] / len(yr_range)
            u += hfile['u'][:] / len(yr_range)
            T += hfile['T'][:] / len(yr_range)
            Th1 += hfile['Th1'][:] / len(yr_range)
            uv += hfile['uv'][:] / len(yr_range)
            vT += hfile['vT'][:] / len(yr_range)
            TT += hfile['TT'][:] / len(yr_range)
            eke += hfile['eke'][:] / len(yr_range)
            u_vrtcl_flux += hfile['u_vrtcl_flux'][:] / len(yr_range)
            pot_temp_vrtcl_flux += hfile['pot_temp_vrtcl_flux'][:] / len(yr_range)
            z1 += hfile['z1'][:] / len(yr_range)
            T1 += hfile['T1'][:] / len(yr_range)
            u1 += hfile['u1'][:] / len(yr_range)
            lh_flux_sfc += hfile['lh_flux_sfc'][:] / len(yr_range)
            sh_flux_sfc += hfile['sh_flux_sfc'][:] / len(yr_range)
            T_sfc += hfile['T_sfc'][:] / len(yr_range)
            u_02 += hfile['u_02'][:] / len(yr_range)
            v_02 += hfile['v_02'][:] / len(yr_range)
            hfile.close()

    return lat, lon, sig, v, u, T, Th1, uv, vT, TT, eke, u_vrtcl_flux, pot_temp_vrtcl_flux, z1, T1, u1, lh_flux_sfc, sh_flux_sfc, T_sfc, u_02, v_02

# save variables


exp_name, run_name = ("uf_tests_hpc_businger_albedo_eqnx", "uf_tests_hpc_businger_albedo_eqnx_albedo02rcrit40_0")
save_basic_diags(exp_name, run_name, day_list, yr_range, save_variable_list, h5_dir = "data/", h5_name = "basic_vars")

exp_name, run_name = ("uf_tests_hpc_albedo_eqnx", "uf_tests_hpc_albedo_eqnx_albedo02_1")
save_basic_diags(exp_name, run_name, day_list, yr_range, save_variable_list, h5_dir = "data/", h5_name = "basic_vars_c")

# bus run
lat, lon, sig, v, u, T, Th1, uv, vT, TT, eke, u_vrtcl_flux, pot_temp_vrtcl_flux, z1, T1, u1, lh_flux_sfc, sh_flux_sfc, T_sfc, u_02, v_02 = read_basic_diags(yr_range, h5_name = "basic_vars")

# control run
lat_c, lon_c, sig_c, v_c, u_c, T_c, Th1_c, uv_c, vT_c, TT_c, eke_c, u_vrtcl_flux_c, pot_temp_vrtcl_flux_c, z1_c, T1_c, u1_c, lh_flux_sfc_c, sh_flux_sfc_c, T_sfc_c, u_02_c, v_02_c  = read_basic_diags(yr_range, h5_name = "basic_vars_c")


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
plt.savefig("circ1_eqnx_allyrs.png")

plt.figure(figsize=(20,4))
zonal_plot(lat, sig, np.nanmean(u,axis=0) -  np.nanmean(u_c,axis=0), title = 'U', subplot_no = (1,4,1))
zonal_plot(lat, sig, np.nanmean(T,axis=0) -  np.nanmean(T_c,axis=0), title = 'T', subplot_no = (1,4,2))
zonal_plot(lat, sig, np.nanmean(TT,axis=0) -  np.nanmean(TT_c,axis=0), title = 'TT', subplot_no = (1,4,3))
pressure = sig * 100000.0 # units: Pa
zonal_plot(lat, sig[:-1], np.nanmean(psi(v, pressure ,lat),axis=0) - np.nanmean(psi(v_c, pressure ,lat),axis=0), title = 'meridional streamfunction', subplot_no = (1,4,4))
plt.savefig("circ1_eqnx_diff_allyrs.png")

# vT, uv, u^2+v^2, TT
plt.figure(figsize=(20,4))
zonal_plot(lat, sig, np.nanmean(uv,axis=0), title = 'uv', subplot_no = (1,4,1) , vrange = np.arange(-66,67,6))
zonal_plot(lat, sig, np.nanmean(vT,axis=0), title = 'vT', subplot_no = (1,4,2) , vrange = np.arange(-18,19,2))
zonal_plot(lat, sig, np.nanmean(0.5*eke,axis=0), title = 'eke', subplot_no = (1,4,3), vrange = np.arange(0,330,15))
plt.savefig("circ2_eqnx_allyrs.png")

# diff
plt.figure(figsize=(20,4))
zonal_plot(lat, sig, np.nanmean(uv,axis=0) - np.nanmean(uv_c,axis=0), title = 'uv', subplot_no = (1,4,1) )
zonal_plot(lat, sig, np.nanmean(vT,axis=0) - np.nanmean(vT_c,axis=0), title = 'vT', subplot_no = (1,4,2) )
zonal_plot(lat, sig, np.nanmean(0.5*eke,axis=0) - np.nanmean(0.5*eke_c,axis=0), title = 'eke', subplot_no = (1,4,3))
plt.savefig("circ2_eqnx_diff_allyrs.png")


# turbulent surface fluxes
plt.figure(figsize=(20,4))
zonal_plot(lon, lat, (np.nanmean(lh_flux_sfc,axis=0)), title = 'latent (W/m2)', subplot_no = (1,4,1) , vrange = np.arange(-50,51,5))
zonal_plot(lon, lat, np.nanmean(sh_flux_sfc,axis=0), title = 'sensible (W/m2)', subplot_no = (1,4,2) , vrange = np.arange(-50,51,5))
plt.subplot(1,4,3)
plt.plot(lat, np.nanmean(lh_flux_sfc,axis=(0, -1)), label = "bus")
plt.plot(lat, np.nanmean(lh_flux_sfc_c,axis=(0, -1)), label = "ctrl")
plt.legend()
plt.title('latent (W/m2)')
plt.subplot(1,4,4)
plt.plot(lat, np.nanmean(sh_flux_sfc,axis=(0, -1)), label = "bus")
plt.plot(lat, np.nanmean(sh_flux_sfc_c,axis=(0, -1)), label = "ctrl")
plt.legend()
plt.title('sensible (W/m2)')
plt.savefig("test_data_albedo_eqnx_comparison_allyrs.png")

# differences of the above figures
plt.figure(figsize=(20,4))
zonal_plot(lon, lat, (np.nanmean(lh_flux_sfc,axis=0) - np.nanmean(lh_flux_sfc_c,axis=0)), title = 'latent (W/m2)', subplot_no = (1,4,1) )
zonal_plot(lon, lat, np.nanmean(sh_flux_sfc,axis=0) - np.nanmean(sh_flux_sfc_c,axis=0), title = 'sensible (W/m2)', subplot_no = (1,4,2) )
plt.subplot(1,4,3)
plt.plot(lat, (np.nanmean(lh_flux_sfc,axis=(0,-1)) - np.nanmean(lh_flux_sfc_c,axis=(0,-1))))
plt.title('latent (W/m2)')
plt.subplot(1,4,4)
plt.plot(lat, np.nanmean(sh_flux_sfc,axis=(0,-1)) - np.nanmean(sh_flux_sfc_c,axis=(0,-1)))
plt.title('sensible (W/m2)')
plt.savefig("test_data_albedo_eqnx_comparison_diff_allyrs.png")

# stability
T_sfc_c[242,:,:] = T_sfc_c[241,:,:]
ri_bulk_c = 9.81 * (Th1_c[:, :,:] - T_sfc_c[:, :,:]) * z1_c[:,:,:]/(T_sfc_c[:, :,:]*u1[:,:,:]**2)
ri_bulk = 9.81 * (Th1[:, :,:] - T_sfc[:, :,:]) * z1[:,:,:]/(T_sfc[:, :,:]*u1[:,:,:]**2)

plt.figure(figsize=(20,4))
# histogram of ri_bulk
zonal_plot(lon, lat, (np.nanmean(ri_bulk_c,axis=0)), title = 'ri_bulk_c', subplot_no = (1,4,1) )
zonal_plot(lon, lat, np.nanmean(ri_bulk,axis=0), title = 'ri_bulk', subplot_no = (1,4,2) )
plt.subplot(1,4,3)
plt.hist(ri_bulk_c.flatten(), bins = np.arange(-50,50,1))
plt.ylabel('count')
# plt.ylim([0,0.2e6])
plt.subplot(1,4,4)
plt.hist(ri_bulk.flatten(), bins = np.arange(-50,50,1))
plt.xlabel('ri_bulk')
# plt.ylim([0,0.2e6])
plt.savefig("hist_albedo_eqnx_allyrs.png")

plt.figure(figsize=(20,4))
# histogram of ri_bulk
zonal_plot(lon, lat, (np.nanmean(ri_bulk_c,axis=0)), title = 'ri_bulk_c', subplot_no = (1,4,1) )
plt.plot([0,360], [70,70])
zonal_plot(lon, lat, np.nanmean(ri_bulk,axis=0), title = 'ri_bulk', subplot_no = (1,4,2) )
plt.subplot(1,4,3)
plt.hist(ri_bulk_c[:,119:].flatten(), bins = np.arange(-100,100,1))
plt.ylabel('count')
plt.xlabel('ri_bulk_c')
plt.ylim([0,0.02e6])
plt.subplot(1,4,4)
plt.hist(ri_bulk[:,119:].flatten(), bins = np.arange(-100,100,1))
plt.xlabel('ri_bulk')
plt.ylim([0,0.02e6])
plt.savefig("hist_albedo_eqnx_polar_allyrs.png")


plt.figure(figsize=(20,4))
# histogram of ri_bulk
zonal_plot(lon, lat, (np.nanmean(ri_bulk_c,axis=0)), title = 'ri_bulk_c', subplot_no = (1,4,1) )
plt.plot([0,360], [-75,-75])
zonal_plot(lon, lat, np.nanmean(ri_bulk,axis=0), title = 'ri_bulk', subplot_no = (1,4,2) )
plt.subplot(1,4,3)
plt.hist(ri_bulk_c[:,:10].flatten(), bins = np.arange(-100,100,1))
plt.ylabel('count')
plt.xlabel('ri_bulk_c')
plt.ylim([0,0.02e6])
plt.subplot(1,4,4)
plt.hist(ri_bulk[:,:10].flatten(), bins = np.arange(-100,100,1))
plt.xlabel('ri_bulk')
plt.ylim([0,0.02e6])
plt.savefig("hist_albedo_eqnx_polar_sh_allyrs.png")






