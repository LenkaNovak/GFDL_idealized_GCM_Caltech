import scipy.interpolate as si

def calculate_cospectra(u, v, lat, lon, t, dt_secs = 86400.0, fdir=".", fname="cospectra"):

    """
    Calculate cospectra for a given time series of u and v and save to hdf5 file.

    Input:
    - u, v: time series of zonal and meridional wind [m/s] (t, lon)
    - lat, lon: latitude and longitude [deg]
    - t: time [s]
    - dt_secs: time step [s]
    - fdir: directory to save file
    - fname: file name to save file
    """
    radlon = lon /180.*np.pi

    v = (v - np.mean(v, axis=-1)[:,None])# / np.std(v)
    u = (u - np.mean(u, axis=-1)[:,None])# / np.std(u)

    dlon = np.abs(radlon[1]-radlon[2])
    nt, nlon = np.shape(u)

    U = np.fft.fft( u, axis = 1 ) / nt / nlon
    V = np.fft.fft( v, axis = 1 ) / nt / nlon

    nm = int(nlon/2) # no. positive wavenumbers
    m =  np.fft.fftfreq( nlon, d = dlon )[:nm]

    CU = U[:, :nm].real # cosines
    SU = U[:, :nm].imag # sines
    CV = V[:, :nm].real
    SV = V[:, :nm].imag

    NFFT = 2**8 # 256 ...2^ by default
    nw = int(NFFT / 2 + 1) # no. positive frequencies

    #K_w,n arrays
    K_p = np.zeros( ( nw, nm ) )
    K_n = np.zeros( ( nw, nm ) )

    #Cross-spectra density using Welsch method
    for m_i in range( nm ): # for each wavenumber
        csd_CUCV, om = plt.csd( CU[:, m_i], CV[:, m_i], Fs = 1. / dt_secs, NFFT = NFFT, scale_by_freq = True, detrend = "linear")
        csd_SUSV, om = plt.csd( SU[:, m_i], SV[:, m_i], Fs = 1. / dt_secs, NFFT = NFFT, scale_by_freq = True, detrend = "linear")
        csd_CUSV, om = plt.csd( CU[:, m_i], SV[:, m_i], Fs = 1. / dt_secs, NFFT = NFFT, scale_by_freq = True, detrend = "linear")
        csd_SUCV, om = plt.csd( SU[:, m_i], CV[:, m_i], Fs = 1. / dt_secs, NFFT = NFFT, scale_by_freq = True, detrend = "linear")
        K_p[:, m_i] = csd_CUCV.real + csd_SUSV.real + csd_CUSV.imag - csd_SUCV.imag
        K_n[:, m_i] = csd_CUCV.real + csd_SUSV.real - csd_CUSV.imag + csd_SUCV.imag
            #Don't need factor 4 from Hayashi eq4.11, since Fourier co-efficients are 1/2 as large due to only retaining positive wavenumbers

        #Combine
        K_combine = np.zeros( ( nw * 2, nm ) )
        K_combine[:nw, :] = K_n[::-1, :]   #for the convolution
        K_combine[nw:, :] = K_p[:, :]

        width = 5.0
        smooth = True
        if smooth == True:
            x = np.linspace( -nw / 2, nw / 2., nw )
            gauss_filter = np.exp( -x ** 2 / (2. * width ** 2 ) )
            gauss_filter /= sum( gauss_filter )
            for i in range( nm ):
                K_combine[:, i] = np.convolve( K_combine[:, i], gauss_filter, 'same' )

    ## Save variables
    trunc = 1 # 149
    hfile = h5.File("%s/%s.hdf5"%(fdir, fname),'w')
    hfile.create_dataset('w', data = np.append(-np.flipud(range(nw)),range(nw))[:],compression = 4)
    hfile.create_dataset('m', data = range(nm)[:],compression = 4)
    hfile.create_dataset('time', data = t,compression = 4)
    hfile.create_dataset('omega', data = np.append(-np.flipud(om),om),compression = 4)
    hfile.create_dataset('K_wm', data = K_combine[:, :],compression = 4)
    hfile.create_dataset('lat', data = lat,compression = 4)
    hfile.close()

def load_cospectra_file(fdir, fname):
    file_path = "%s/%s.hdf5"%(fdir, fname)
    hfile = h5.File(file_path, 'r+')
    m = hfile['m'][:]
    w = hfile['omega'][:]
    K_wm = hfile['K_wm'][:][np.newaxis,:,:]
    time = hfile['time'][:][np.newaxis,:]
    lat = hfile['lat'][:]
    return m, w, K_wm, time, lat

def calPhaseSpeedSpectrum( P_p, P_n, m, w, l4, cmax, nc, m_min = 1, m_max = 20 ):
	"""
    from nick lutsko's website:
	Calculate space-time co-spectra, following method of Hayashi (1971)

	Input:
	  P_p - spectra for positive phase speeds
	  P_n - spectra for negative phase speeds
	  m - wavenumbers
	  w - frequencies
      l4 - latitude
	  cmax - maximum phase speed
	  nc - size of phase speed grid
	  m_min - lowest wave number to sum over
      m_max - highest wave number to sum over

	Output:
	  P_cp - spectra for positive phase speeds
	  P_cn - spectra for negative phase speeds
	  C * lon_unit / time_unit - phase speeds
	"""
	if m_max < m_min:
		print("WARNING: highest wavenumber smaller than lowest wavenumber")

	j = len( m )
	t = len( w )

	#Make phase speed grid
	C = np.linspace(0., cmax, nc)

	#K_n,c arrays
	P_cp = np.zeros( ( nc, j ) )
	P_cn = np.zeros( ( nc, j ) )

	#Interpolate
	for i in range( m_min, m_max ):
		#Make interpolation functions
        # k = 2pi / L = 2pi / (2pi / (m) * a cos th) = m / a cos th
        # c = w / k = w a cos th / m
		f1 = si.interp1d(w / m[i] * 6371000. * np.cos(l4*np.pi/180.), P_p[:, i], 'linear' )
		f2 = si.interp1d(w / m[i] * 6371000. * np.cos(l4*np.pi/180.), P_n[:, i], 'linear' )

		#interp1d doesn't handle requested points outside data range well, so just zero out these points
		k = -1
		for j in range( len(C) ):
			if C[j] > (max(w) / m[i] * 6371000. * np.cos(l4*np.pi/180.) ):
				k = j
				break
		if k == -1:
			k = len( C )

		ad1 = np.zeros( nc )
		ad1[:k] =  f1( C[:k]  )
		ad2 = np.zeros( nc )
		ad2[:k] =  f2( C[:k] )

		#Interpolate
		P_cp[:, i] = ad1 #* m[i]
		P_cn[:, i] = ad2 #* m[i]

 	#Sum over all wavenumbers
	return P_cp, P_cn, C





def phase_speed_spectra(l4, nc = 60, cmax = 20, fdir = ".", fname = "cospectra"):

    m, w, K_wm_av, time, lat = load_cospectra_file(fdir, fname)

    K_wm_av = K_wm_av[0,:,:]
    m_max = int(nc/2)
    w_rad = w *2. * np.pi
    nw = int(len(w) / 2)

    cp_cm_av = np.zeros(( int(2*nc), len(m) ))

    # nw = 129

    K_n = K_wm_av[:nw, :][::-1,:]
    K_p = K_wm_av[nw:, :]
    P_Cp, P_Cn, cp = calPhaseSpeedSpectrum(K_p, K_n, m, w_rad[nw:], l4, cmax, nc, int(1), m_max)
    cp_cm_av[:nc, :] = P_Cn[::-1,:] #negative phase speeds
    cp_cm_av[nc:, :] = P_Cp[:,:] #positive phase speeds

    return cp_cm_av, cp, m



# # calculate cospectra for several latitudinal bands at 200 hPa
# U_all, V_all, T_all, Z_all, p, lat, lon, t = e5.load_plev_data(ncfile_name)

# p4 = 200
# for l4 in lat[::15]:

#     la_i = np.where(lat ==l4)[0][0]
#     p_i = np.where(np.array(level_list) == p4)[0][0]

#     v = V_all[:,p_i,la_i,:] # t, lon
#     u = U_all[:,p_i,la_i,:]



#     calculate_cospectra(u, v, lat, lon, t, dt_secs = 86400.0, fdir=".", fname="cospectra_lat_%00d"%l4)

# # calculate phase speed spectra for several latitudinal bands at 200 hPa
# # cp_cm_av = np.zeros((len(yr_list_array), len(lat_sampled), len(mo_list_array), int(2*nps), len(m) ))
# nc = 60
# nm = len(lon)
# cp_cm_av = np.zeros((len(lat_sampled), int(2*nc), len(nm) ))
# for l4 in lat[::15]:

#     cp_cm_av, cp, m = phase_speed_spectra(nc = nc, cmax = 20, fdir = ".", fname = "cospectra_lat_%00d"%l4)




# c = plt.contourf(np.arange(-nc,nc),m,np.transpose(cp_cm_av))
#     # if cbar ==True:
#     #     plt.colorbar(c)



