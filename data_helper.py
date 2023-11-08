def print_ncinfo(dataset):
    for var in dataset.variables:
        var1 = dataset.variables[var]
        for name in var1.ncattrs():
            print(var, '-->',name, '=', getattr(var1,name))

def append_month_days(fn, mo_i, yr_i, hw, yr_list, ncinfo = False):
    """
    - extract momthly data from FMS daily output
    """
    print(fn)
    ncdata =Dataset(fn%(yr_list[yr_i], 1),'r')
    if ncinfo == True:
        print_ncinfo(ncdata)
    lat = ncdata.variables['lat'][:][::-1]
    lon = ncdata.variables['lon'][:]
    p = ncdata.variables['sigma'][:][:][::-1] * 1.0e5
    ncdata.close()
    U = np.zeros((30+hw*2, len(p), len(lat), len(lon)))
    V = np.zeros((30+hw*2, len(p), len(lat), len(lon)))
    T = np.zeros((30+hw*2, len(p), len(lat), len(lon)))
    Z = np.zeros((30+hw*2, len(p), len(lat), len(lon)))
    t = np.arange(mo_i*30-hw,(mo_i+1)*30+hw)

    for d_i, d in enumerate(t):
        if d < 0:
            ncdata = Dataset(fn%(yr_list[yr_i - 1], int(360. + d+1.0) ),'r')
        elif d > 359:
            ncdata = Dataset(fn%(yr_list[yr_i + 1], int(d+1.0 - 360.) ),'r')
        else:
            ncdata = Dataset(fn%(yr_list[yr_i], int(d+1.0) ),'r')
        U[d_i,:,:,:]=ncdata.variables['u'][:][:,:,:][::-1,::-1,:]
        V[d_i,:,:,:]=ncdata.variables['v'][:][:,:,:][::-1,::-1,:]
        T[d_i,:,:,:]=ncdata.variables['temp'][:][:,:,:][::-1,::-1,:]
        Z[d_i,:,:,:]=ncdata.variables['z'][:][:,:,:][::-1,::-1,:]
        ncdata.close()
    return U, V, T, Z, p, lat, lon, t


def append_year_days_onevar(fn, day_list, var_name, enable_print_ncinfo=False):
    ncdata =Dataset(fn%1,'r')
    if enable_print_ncinfo==True:
        print_ncinfo(ncdata)
    zon_waven = ncdata.variables['zon_waven'][:]
    legendre = ncdata.variables['legendre'][:]
    lat = ncdata.variables['lat'][:]
    lon = ncdata.variables['lon'][:]
    sig=ncdata.variables['sigma'][:]
    var_=ncdata.variables[var_name][:]
    ncdata.close()
    var = np.zeros((np.shape(var_)))[np.newaxis, ...]
    for d_i, d in enumerate(day_list):
        ncdata = Dataset(fn%d,'r')
        var=np.append(var, (ncdata.variables[var_name][:])[np.newaxis, ...], axis = 0)
        ncdata.close()

    return zon_waven, legendre, lat, lon, sig, var


def append_year_days_zonal_mean(fn, day_list, eddy_stats = True, enable_print_ncinfo=False):
    ncdata =Dataset(fn%1,'r')
    if enable_print_ncinfo==True:
        print_ncinfo(ncdata)
    lat = ncdata.variables['lat'][:]
    lon = ncdata.variables['lon'][:]
    sig=ncdata.variables['sigma'][:]
    ncdata.close()
    v = np.zeros((len(day_list), len(sig), len(lat)))
    u = np.zeros((len(day_list), len(sig), len(lat)))
    T = np.zeros((len(day_list), len(sig), len(lat)))
    Th = np.zeros((len(day_list), len(sig), len(lat)))
    uv = np.zeros((len(day_list), len(sig), len(lat)))
    vT = np.zeros((len(day_list), len(sig), len(lat)))
    TT = np.zeros((len(day_list), len(sig), len(lat)))
    eke = np.zeros((len(day_list), len(sig), len(lat)))
    # mom_drag_coeff = np.zeros((len(day_list), len(sig), len(lat)))
    # temp_drag_coeff  = np.zeros((len(day_list), len(sig), len(lat)))
    u_vrtcl_flux = np.zeros((len(day_list), len(lat), len(lon)))
    pot_temp_vrtcl_flux = np.zeros((len(day_list), len(lat), len(lon)))
    z_1 = np.zeros((len(day_list), len(lat), len(lon)))
    T_1 = np.zeros((len(day_list), len(lat), len(lon)))
    u_1 = np.zeros((len(day_list), len(lat), len(lon)))
    Th_1 = np.zeros((len(day_list), len(lat), len(lon)))
    lh_flux_sfc = np.zeros((len(day_list), len(lat), len(lon)))
    sh_flux_sfc = np.zeros((len(day_list), len(lat), len(lon)))
    T_sfc = np.zeros((len(day_list), len(lat), len(lon)))

    for d_i, d in enumerate(day_list):
        ncdata = Dataset(fn%d,'r')
        v[d_i,:,:]=np.mean(ncdata.variables['v'][:],axis=-1)
        u[d_i,:,:]=np.mean(ncdata.variables['u'][:],axis=-1)
        T[d_i,:,:]=np.mean(ncdata.variables['temp'][:],axis=-1)
        # Th[d_i,:,:]=np.mean(ncdata.variables['pot_temp'][:],axis=-1)
        T_sfc[d_i,:,:] = ncdata.variables['temp_sfc'][:][:,:]
        u_vrtcl_flux[d_i,:,:]=ncdata.variables['u_vrtcl_flux'][:][-1,:,:]
        pot_temp_vrtcl_flux[d_i,:,:]=ncdata.variables['pot_temp_vrtcl_flux'][:][-1,:,:]
        pot_temp_vrtcl_flux[d_i,:,:]=ncdata.variables['pot_temp_vrtcl_flux'][:][-1,:,:]
        lh_flux_sfc[d_i,:,:]=ncdata.variables['lh_flux_sfc'][:][:,:]
        sh_flux_sfc[d_i,:,:]=ncdata.variables['sh_flux_sfc'][:][:,:]
        z_1[d_i,:,:]=ncdata.variables['z'][:][-1,:,:]
        T_1[d_i,:,:]=ncdata.variables['temp'][:][-1,:,:]
        u_1[d_i,:,:]=ncdata.variables['u'][:][-1,:,:]
        Th_1[d_i,:,:]=ncdata.variables['pot_temp'][:][-1,:,:]
        if eddy_stats ==True:
            vST = ncdata.variables['v'][:] - v[d_i,:,:,None]
            uST = ncdata.variables['u'][:] - u[d_i,:,:,None]
            TST = ncdata.variables['temp'][:] - T[d_i,:,:,None]
            uv[d_i,:,:]=np.mean(uST*vST,axis=-1)
            vT[d_i,:,:]=np.mean(vST*TST,axis=-1)
            TT[d_i,:,:]=np.mean(TST*TST,axis=-1)
            eke[d_i,:,:]=np.mean(uST**2+vST**2,axis=-1)
        ncdata.close()

    return lat, lon, sig, v, u, T, Th_1, uv, vT, TT, eke, u_vrtcl_flux, pot_temp_vrtcl_flux, z_1, T_1, u_1, lh_flux_sfc, sh_flux_sfc, T_sfc