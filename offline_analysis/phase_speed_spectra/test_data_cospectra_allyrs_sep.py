# exec(open("test_data_cospectra_allyrs_sep.py").read())

#modules
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import integrate
import h5py as h5
from scipy import ndimage
import os

# global parameters:
_a = 6371000.0
_g = 9.81

# read in data_helper.py
exec(open("../helpers/data_helper.py").read())
exec(open("../helpers/cospectra_helper.py").read())

day_list = np.arange(1,361)

yr_range = ["3600", "3960", "4320", "4680", "5040", "5400", "5760", "6120", ] # [ "0720", "1080", "1440","1800"] # "%04d"%day_list[-1] # 1440 possibly corrupted for rcrit8 !
# yr_range = ["1440", "4320", "5400"]

save_variable_list = ['lat', 'lon', 'sig', 'v', 'u', 'T', 'Th1', 'uv', 'vT', 'TT', 'eke', 'u_vrtcl_flux', 'pot_temp_vrtcl_flux', 'z1', 'T1', 'u1', 'T_sfc', 'u_02', 'v_02']

# Save variables
def save_basic_diags(exp_name, run_name, day_list, yr_range, save_variable_list, h5_dir = "../basic_plots/data/", h5_name = "basic_vars"):

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
def read_basic_diags(yr, h5_dir = "../basic_plots/data/", h5_name = "basic_vars"):

    # create empty arrays for all variables
    lat, lon, sig = [np.array([]) for i in range(3)]

    hfile = h5.File("%s/%s_%s.hdf5"%(h5_dir, h5_name, yr_range[0]),'r')
    lat = hfile['lat'][:]
    lon = hfile['lon'][:]
    sig = hfile['sig'][:]
    hfile.close()

    # v, u, T, uv, vT, TT, eke = [np.zeros((360, len(sig), len(lat))) for i in range(7)]

    Th1, z1, T1, u1, lh_flux_sfc, sh_flux_sfc, T_sfc, u_02, v_02, u_vrtcl_flux, pot_temp_vrtcl_flux = [np.zeros((360, len(lat), len(lon))) for i in range(11)]

    # read in variables for all years
    hfile = h5.File("%s/%s_%s.hdf5"%(h5_dir, h5_name, yr),'r')
    u_02 = hfile['u_02'][:]
    v_02 = hfile['v_02'][:]
    hfile.close()

    # return lat, lon, sig, v, u, T, Th1, uv, vT, TT, eke, u_vrtcl_flux, pot_temp_vrtcl_flux, z1, T1, u1, lh_flux_sfc, sh_flux_sfc, T_sfc, u_02, v_02
    return lat, lon, sig, u_02, v_02

def get_phase_speed_spectra_v_latitude(u_02, v_02, day_list, lat, lon, lat_sampled, fdir = ".", nc = 60, cmax = 60, nm= 128):
    yr_range_name = "%s-%s"%(yr_range[0], yr_range[-1])
    t = np.arange(0, len(day_list))


    for l4 in lat_sampled:
        fname = "cospectra_lat_%00d_%s"%(l4, yr_range_name)
        la_i = np.where(lat ==l4)[0][0]
        v_lev = v_02[:,la_i,:] # t, lon
        u_lev = u_02[:,la_i,:]
        if os.path.exists("%s/%s.hdf5"%(fdir, fname)):
            continue
        else:
            calculate_cospectra(u_lev, v_lev, lat, lon, t, dt_secs = 86400.0, fdir=fdir, fname=fname)

    # get phase speed spectra
    K_wm = np.zeros((len(lat_sampled), int(2*nc), nm ))
    for i, l4 in enumerate(lat_sampled):
        fname = "cospectra_lat_%00d_%s"%(l4, yr_range_name)
        K_wm[i, :, :], cp, m = phase_speed_spectra(l4, nc = nc, cmax = cmax, fdir =fdir, fname = fname)
    return K_wm, cp, m, lat_sampled


# save variables
rctit = "40"

exp_name, run_name = ("uf_tests_hpc_businger_albedo_eqnx", "uf_tests_hpc_businger_albedo_eqnx_albedo02rcrit%s_0"%rctit)
save_basic_diags(exp_name, run_name, day_list, yr_range, save_variable_list, h5_dir = "../basic_plots/data_rcrit%s/"%rctit, h5_name = "basic_vars")

# exp_name, run_name = ("uf_tests_hpc_albedo_eqnx", "uf_tests_hpc_albedo_eqnx_albedo02_1")
# # exp_name, run_name = ("uf_tests_hpc_businger_albedo_eqnx", "uf_tests_hpc_businger_albedo_eqnx_albedo02rcrit8_0")
# save_basic_diags(exp_name, run_name, day_list, yr_range, save_variable_list, h5_dir = "../basic_plots/data/", h5_name = "basic_vars_c")

# ricurb
exp_name, run_name = ("uf_tests_hpc_albedo_eqnx_ricurb", "uf_tests_hpc_albedo_eqnx_ricurb_albedo02stOption2zetatrans8_0")
save_basic_diags(exp_name, run_name, day_list, yr_range, save_variable_list, h5_dir = "../basic_plots/data_Option2zetatrans8_0/", h5_name = "basic_vars")

for yr_i, yr in enumerate(yr_range):

    # bus run
    # lat, lon, sig, v, u, T, Th1, uv, vT, TT, eke, u_vrtcl_flux, pot_temp_vrtcl_flux, z1, T1, u1, lh_flux_sfc, sh_flux_sfc, T_sfc, u_02, v_02 = read_basic_diags(yr_range, h5_name = "basic_vars")
    lat, lon, sig, u_02, v_02 = read_basic_diags(yr, h5_dir = "../basic_plots/data_rcrit%s/"%rctit, h5_name = "basic_vars")

    # control run
    # lat_c, lon_c, sig_c, v_c, u_c, T_c, Th1_c, uv_c, vT_c, TT_c, eke_c, u_vrtcl_flux_c, pot_temp_vrtcl_flux_c, z1_c, T1_c, u1_c, lh_flux_sfc_c, sh_flux_sfc_c, T_sfc_c, u_02_c, v_02_c  = read_basic_diags(yr_range, h5_name = "basic_vars_c")
    # lat_c, lon_c, sig_c, u_02_c, v_02_c  = read_basic_diags(yr, h5_dir = "../basic_plots/data/", h5_name = "basic_vars_c")
    lat_c, lon_c, sig_c, u_02_c, v_02_c  = read_basic_diags(yr, h5_dir = "../basic_plots/data_Option2zetatrans8_0/", h5_name = "basic_vars")



    if yr_i == 0:
        lat_sampled_hemi = lat[lat>0][::5]
        lat_sampled= np.sort(np.append(lat_sampled_hemi, -lat_sampled_hemi))
        nc = 60
        nm = int(len(lon) / 2)
        K_wm = np.zeros((len(yr_range), len(lat_sampled), int(2*nc), nm ))
        K_wm_c = np.zeros((len(yr_range), len(lat_sampled), int(2*nc), nm ))


    K_wm[yr_i,:,:,:], cp, m, lat_sampled = get_phase_speed_spectra_v_latitude(u_02, v_02, day_list, lat, lon, lat_sampled, fdir ="data_rcrit%s/"%rctit, nc = nc, cmax = 60, nm= nm)
    K_wm_c[yr_i,:,:,:], cp, m, lat_sampled = get_phase_speed_spectra_v_latitude(u_02_c, v_02_c, day_list, lat, lon, lat_sampled, fdir ="data/", nc = nc, cmax = 60, nm= nm)


cp_2sided = np.append(-np.flipud(cp),cp)

plt.figure()
plt.subplot(121)
c = plt.contourf(lat_sampled,cp_2sided,np.transpose(np.mean(np.nansum(K_wm,axis=-1), axis=0)))
plt.colorbar(c)
plt.plot(lat, np.mean(u_02, axis=(0,-1)))
plt.subplot(122)
c = plt.contourf(lat_sampled,cp_2sided,np.transpose(np.mean(np.nansum(K_wm_c,axis=-1), axis=0)))
plt.colorbar(c)
plt.plot(lat_c, np.mean(u_02_c, axis=(0,-1)))
plt.savefig("ps_spectrum.png")

plt.figure()
plt.subplot(111)
c = plt.contourf(lat_sampled,cp_2sided,np.transpose(np.mean(np.nansum(K_wm,axis=-1), axis=0)) - np.transpose(np.mean(np.nansum(K_wm_c,axis=-1), axis=0)), np.arange(-15,15.1,0.1), cmap =  'RdBu_r')
plt.plot(lat_c, np.mean(u_02_c, axis=(0,-1)))
plt.colorbar(c)
plt.xlim([-85, 85])
plt.ylim([-20, 60])
plt.savefig("ps_spectrum_diff.png")

# hemispheric diff
plt.figure(figsize=(4,4))
plt.subplot(111)
K_wm_hemi = (K_wm[:,lat_sampled>0,:,:] - np.flip(K_wm[:,lat_sampled<0,:,:], axis=1)) / 2
K_wm_c_hemi = (K_wm_c[:,lat_sampled>0,:,:] - np.flip(K_wm_c[:,lat_sampled<0,:,:], axis=1)) / 2
c = plt.contourf(lat_sampled[lat_sampled>0],cp_2sided, (np.transpose(np.mean(np.nansum(K_wm_hemi,axis=-1), axis=0)) - np.transpose(np.mean(np.nansum(K_wm_c_hemi,axis=-1), axis=0))) / 2. / np.pi, np.arange(-67.5,70.1,5) / 2. / np.pi , cmap =  'RdBu_r')
plt.plot(lat[lat>0], np.mean(u_02[:,lat>0,:], axis=(0,-1)), color = 'k')
plt.colorbar(c)
plt.xlim([0, 85])
plt.ylim([-20, 60])
plt.xlabel("Latitude (deg)")
plt.ylabel("m/s")
plt.tight_layout()
plt.savefig("ps_spectrum_diff_hemi.png")




# sanity check
uv_c = np.mean((u_02_c - np.mean(u_02_c, axis=(-1))[:,:,None]) * (v_02_c - np.mean(v_02_c, axis=(-1))[:,:,None]), axis=-1)
uv =  np.mean((u_02 - np.mean(u_02, axis=(-1))[:,:,None]) * (v_02 - np.mean(v_02, axis=(-1))[:,:,None]), axis=-1)

plt.figure()
plt.subplot(121)
c = plt.plot(lat,(np.nanmean(uv_c,axis=-0)) )
c = plt.plot(lat,(np.nanmean(uv,axis=-0)) - (np.nanmean(uv_c,axis=-0)), label = "from eulerian")
plt.plot(lat_sampled,np.nansum((np.mean(np.nansum(K_wm,axis=-1), axis=0)) - (np.mean(np.nansum(K_wm_c,axis=-1), axis=0)), axis=-1) , label = "from ps spectrum")
plt.plot(lat_c, np.mean(u_02_c, axis=(0,-1)))
plt.xlim([-85, 85])
plt.legend()
plt.savefig("uv_diff.png")


    # c = plt.contourf(w,m,np.transpose(np.mean(K_wm,axis=0)), rrnge)
    # if cbar ==True:
    #     plt.colorbar(c)



