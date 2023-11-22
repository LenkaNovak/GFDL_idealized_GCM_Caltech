# exec(open("test_analytical.py").read())

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

yr_range = ["3600"] # [ "0720", "1080", "1440","1800"] # "%04d"%day_list[-1] # 1440 possibly corrupted for rcrit8 !
# yr_range = ["1440", "4320", "5400"]

save_variable_list = ['lat', 'lon', 'sig', 'v', 'u', 'T', 'Th1', 'uv', 'vT', 'TT', 'eke', 'u_vrtcl_flux', 'pot_temp_vrtcl_flux', 'z1', 'T1', 'u1', 'T_sfc', 'u_02', 'v_02']

# Save variables
# ... rm

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
        # if os.path.exists("%s/%s.hdf5"%(fdir, fname)):
        #     continue
        # else:
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

yr = yr_range[0]

# bus run
# lat, lon, sig, v, u, T, Th1, uv, vT, TT, eke, u_vrtcl_flux, pot_temp_vrtcl_flux, z1, T1, u1, lh_flux_sfc, sh_flux_sfc, T_sfc, u_02, v_02 = read_basic_diags(yr_range, h5_name = "basic_vars")
lat, lon, level_list, _, _ = read_basic_diags(yr, h5_dir = "../basic_plots/data_rcrit%s/"%rctit, h5_name = "basic_vars")

t = np.arange(0,360,1) * 60.*60.*24. # [s]
m_a = 5 # same wavenumber for both u and v
lon_per_s = 10. * 180. / 6371000. / np.pi # cp ~ 10 m/s at the equator
# since w = c m / a cos th;
u_02 = np.sin( ((lon[None, None,:] - (1* t[:, None, None] * lon_per_s)) * m_a * 2. * np.pi / 360.0) ) * (np.ones((len(lat)))*np.cos(lat/180*np.pi))[None, :,None]
v_02 = np.sin( ((lon[None, None,:] - (1* t[:, None, None] * lon_per_s)) * m_a * 2. * np.pi / 360.0) ) * (np.ones((len(lat)))*np.cos(lat/180*np.pi))[None, :,None]


lat_sampled_hemi = lat[lat>0][::15]
lat_sampled= np.sort(np.append(lat_sampled_hemi, -lat_sampled_hemi))
nc = 60
nm = int(len(lon) / 2)
K_wm = np.zeros((len(yr_range), len(lat_sampled), int(2*nc), nm ))

K_wm[0,:,:,:], cp, m, lat_sampled = get_phase_speed_spectra_v_latitude(u_02, v_02, day_list, lat, lon, lat_sampled, fdir ="data_analytical/", nc = nc, cmax = 60, nm= nm)

cp_2sided = np.append(-np.flipud(cp),cp)

plt.figure()
plt.contourf(cp_2sided, m, np.transpose(np.nanmean(K_wm[:,len(lat_sampled_hemi),:,:],axis=0)))
plt.xlabel("cp")
plt.ylabel("m")
plt.xlim([-30,30])
plt.ylim([0,20])
plt.savefig("an_test_ps_co_m.png") # check that cp = 10 and m = 5


plt.figure()
plt.subplot(121)
c = plt.contourf(lat_sampled,cp_2sided,np.transpose(np.mean(np.nansum(K_wm,axis=-1), axis=0)))
plt.colorbar(c)
plt.plot(lat, np.mean(u_02, axis=(0,-1)))
plt.savefig("an_test_ps_spectrum.png")


# hemispheric
plt.figure(figsize=(4,4))
plt.subplot(111)
K_wm_hemi = (K_wm[:,lat_sampled>0,:,:] + np.flip(K_wm[:,lat_sampled<0,:,:], axis=1)) / 2 # subtract if zontally asymmetric fluxes
c = plt.contourf(lat_sampled[lat_sampled>0],cp_2sided, (np.transpose(np.mean(np.nansum(K_wm_hemi,axis=-1), axis=0))), cmap =  'RdBu_r')
plt.plot(lat[lat>0], np.mean(u_02[:,lat>0,:], axis=(0,-1)), color = 'k')
plt.colorbar(c)
plt.xlim([0, 85])
plt.ylim([-20, 60])
plt.xlabel("Latitude (deg)")
plt.ylabel("m/s")
plt.tight_layout()
plt.savefig("an_test_ps_spectrum_hemi.png")

# sanity check
uv =  (np.mean((u_02 - np.mean(u_02, axis=(-1))[:,:,None]) * (v_02 - np.mean(v_02, axis=(-1))[:,:,None]), axis=-1))[:,lat>0][:,::15]

plt.figure()
plt.subplot(121)
c = plt.plot(lat_sampled_hemi,(np.nanmean(uv,axis=-0)) , label = "from eulerian")
kw_reduced = np.mean((np.nansum(np.nansum(K_wm,axis=-1), axis=1)), axis=0)[lat_sampled>0] #* m_a / (np.cos(lat_sampled_hemi * np.pi/180.))
plt.plot(lat_sampled_hemi, kw_reduced , label = "from ps spectrum")
plt.plot(lat, np.mean(u_02, axis=(0,-1)))
plt.xlim([-85, 85])
plt.legend()
plt.savefig("an_test_uv.png")


    # c = plt.contourf(w,m,np.transpose(np.mean(K_wm,axis=0)), rrnge)
    # if cbar ==True:
    #     plt.colorbar(c)



