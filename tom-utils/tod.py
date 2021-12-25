


# Utilities for time-ordered data. Mostly for ACT. 


def rot_2d(a):
    return np.array([[np.cos(a),np.sin(a)],[-np.sin(a),np.cos(a)]])

def xy_to_ae(dx,dy,ca=0,ce=0):

y,z = np.matmul()




# I want to standardize what data object I work on. I now use:
# a dict with keys time, azim, elev, data, tags, x, y


def dict_from_enact(tod):






import numpy as np
import scipy as sp

def configure(TOD):
    
    TOD.x, TOD.y  = TOD.point_offset[:,0], TOD.point_offset[:,1]
    TOD.time, TOD.azim, TOD.elev = TOD.boresight

    reg_data, reg_rms, is_reg = regularize(TOD.tod)

    TOD.d_t = np.gradient(TOD.time).mean()
    TOD.azv = np.gradient(TOD.azim) / np.gradient(TOD.time)
    TOD.n_t = len(TOD.time)
    TOD.reg_data = reg_data
    TOD.g = is_reg.copy()

    return TOD

def get_kernel(n,kind='triangle'):
    if kind == 'triangle':
        kernel = np.r_[np.linspace(0,1,n+1)[1:],np.linspace(1,0,n+1)[1:-1]]
        return kernel / kernel.sum()
    if kind == 'flat':
        kernel = np.r_[np.linspace(1,1,n+1)[1:],np.linspace(1,1,n+1)[1:-1]]
        return kernel / kernel.sum()

def downsample(DATA,dsr,axis=-1,kind='triangle'):
    if not dsr > 1:
        return DATA
    
    shaped = True if not len(DATA.shape) > 1 else False
    if shaped:
        DATA = DATA[None,:]
        
    kernel = get_kernel(n=dsr,kind=kind)
    n_kern = len(kernel)

    starts = np.arange(0, DATA.shape[-1] - n_kern, dsr)
    len_ds = len(starts)
    
    ds_data = np.concatenate([np.sum(DATA[:,s:s+n_kern]*kernel[None,:],axis=1)[:,None] for s in starts],axis=1)

    return ds_data if not shaped else ds_data[0]



transfer_gc = lambda d,g_pars=[1],c_pars=[0] : np.poly1d(g_pars)(np.linspace(0,1,len(d)))*d + np.poly1d(c_pars)(np.linspace(0,1,len(d)))

def fit_gc(d,ref,g_order=1,c_order=1,dsr=1):
    
    init_g_pars = np.r_[np.zeros(g_order),1]
    init_c_pars = np.r_[np.zeros(g_order),0]
    
    import scipy as sp
    from scipy import optimize
    
    def transfer(d,*pars):
        pt_ = np.linspace(0,1,len(d))
        return np.poly1d(pars[:g_order+1])(pt_)*d + np.poly1d(pars[g_order+1:])(pt_)
    
    pars,cpars = sp.optimize.curve_fit(transfer,d,ref,
                                       p0=np.r_[init_g_pars,init_c_pars])
    
    return transfer(d,*pars), pars[:g_order+1], pars[g_order+1:]



def regularize(DATA_ARRAY,method='fit',fit_dsr=16,gain_order=1,cnst_order=1):
    
    if method == 'fit':
        
        REG_DATA = DATA_ARRAY.copy()
        ref_data = np.nanmedian(DATA_ARRAY[:,::fit_dsr],axis=0)
        ref_time = np.arange(DATA_ARRAY.shape[-1])
        
        fit_fun  = lambda d, *pars : np.poly1d(pars[:gain_order+1])(d[1])*d[0] + np.poly1d(pars[gain_order+1:])(d[1])
        cal_pars = []
        p0 = [*np.zeros(gain_order),1,*np.zeros(cnst_order),1]
        for i, D in enumerate(DATA_ARRAY):
            cal_pars.append(sp.optimize.curve_fit(fit_fun,np.c_[ref_time[::fit_dsr],D[::fit_dsr]].T,ref_data,p0=p0)[0])
            REG_DATA[i] = fit_fun(np.c_[ref_time,D].T,*cal_pars[-1])

        tot_rms = (REG_DATA-REG_DATA.mean(axis=0)[None,:]).std(axis=-1)
        med_rms = np.percentile(tot_rms,q=25)
        is_regular = (tot_rms < 2 * med_rms)
        
    return REG_DATA, tot_rms, is_regular

def split_dict(dict, durations=[]):

    key = 'azim'
        
    v = np.gradient(dict[key]) 

    # flag wherever the scan velocity changing direction (can be more general)
    flags  = np.r_[0,np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0],len(v)-1]
    splits = np.array([[s,e] for s,e in zip(flags,flags[1:]) if (e-s) > 101]).astype(int)

    # compiles sub-scans that cover the TOD
    sub_splits = splits.copy()
    for i,(s,e) in enumerate(splits):
        split_dur = dict['time'][e] - dict['time'][s]
        for pld in durations:                
            if pld > split_dur:
                continue
            sub_n = 2 * int(split_dur / pld) + 1
            sub_s = np.linspace(s,e-int(pld/self.d_t),sub_n).astype(int)

            for sub_s in np.linspace(s,e-int(np.ceil(pld/self.d_t)),sub_n).astype(int):
                self.sub_splits = np.r_[self.sub_splits,np.array([sub_s,sub_s+int(pld/self.d_t)])[None,:]]
                
    return self

def taper_window(N,end_prop=5e-2):
    
    n_taper = int(N * end_prop)
    return np.r_[.5*(1-np.cos(np.linspace(0,np.pi,n_taper))),np.ones(N-2*n_taper),.5*(1-np.cos(np.linspace(np.pi,0,n_taper)))]




def get_clusters(self,n_clusters=8,filter_type='band-pass', freq_cutoff=[1e-2,1e1], order=3):
    
    from sklearn.cluster import KMeans

    points = np.vstack([self.x[self.g],self.y[self.g]]).T
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(points)

    self.cluster_id = kmeans.labels_
    self.clust_data = np.concatenate([self.reg_data[self.g][self.cluster_id==i].mean(axis=0)[None,:] for i in np.sort(np.unique(self.cluster_id))],axis=0)
    self.clust_x = np.array([self.x[self.g][self.cluster_id==i].mean(axis=0) for i in np.sort(np.unique(self.cluster_id))])
    self.clust_y = np.array([self.y[self.g][self.cluster_id==i].mean(axis=0) for i in np.sort(np.unique(self.cluster_id))])
    self.clust_n = np.array([np.sum(self.cluster_id==i) for i in np.sort(np.unique(self.cluster_id))])
    
    self.outer_x = np.subtract.outer(self.clust_x,self.clust_x)
    self.outer_y = np.subtract.outer(self.clust_y,self.clust_y)
    
    if filter_type == 'band-pass':
        self.filtered_clust_data = band_pass(self.clust_data,*freq_cutoff,1/self.d_t,order=order)
    if filter_type == 'low-pass':
        self.filtered_clust_data = low_pass(self.clust_data,freq_cutoff,1/self.d_t,order=order)
    if filter_type == 'high-pass':
        self.filtered_clust_data = high_pass(self.clust_data,freq_cutoff,1/self.d_t,order=order)

    return self

def get_pl(self,n_clusters=6, sub_scan_durs=[], max_lag=5, order=3, method='outer'):  
    # a list of sub-scan durations to sample the data, e.g. [10,20]
    
    # TOD is a enact-defined object. 
    
    band_pass = lambda data, lc, hc, fs, order : sp.signal.filtfilt(*sp.signal.butter(order,[2*lc/fs,2*hc/fs],btype='band'),data,axis=-1)
    high_pass = lambda data, c, fs, order : sp.signal.filtfilt(*sp.signal.butter(order,2*c/fs,btype='highpass'),data,axis=-1)
    low_pass  = lambda data, c, fs, order : sp.signal.filtfilt(*sp.signal.butter(order,2*c/fs,btype='lowpass'),data,axis=-1)
    

    self.max_di = int(max_lag / self.d_t)
    self.di_sampler = np.linspace(-self.max_di,self.max_di,2*self.max_di+1)
    
    self.nr_us = 1024
    self.us_di_sampler = np.linspace(-self.max_di,self.max_di,2*self.nr_us+1)

    self = get_sub_splits(self,sub_scan_durs)
    self = get_clustered(self,n_clusters=n_clusters)

    # use Fourier methods to compute the pair-lag, for each sub-split and each cluster pair
    self.pl = np.ones((self.sub_splits.shape[0], n_clusters, n_clusters))


    for i_spl,(s,e) in enumerate(self.sub_splits):

        sub_data    = high_pass(self.filtered_clust_data[:,s:e],c=2/((e-s)*self.d_t),fs=1/self.d_t,order=order)*taper_window(e-s,end_prop=1e-1)[None,:]
        #sub_data    = self.filtered_clust_data[:,s:e]*taper_window(e-s,end_prop=1e-1)[None,:]
        ft_sub_data = np.fft.fft(sub_data, axis=-1) 

        for i_det in range(self.pl.shape[2]):
            for j_det in range(i_det+1):
                
                self.dijt = np.real(np.fft.ifft(ft_sub_data[i_det]*np.conj(ft_sub_data[j_det])))
                self.dijt_adj = sp.interpolate.interp1d(self.di_sampler,np.r_[self.dijt[-self.max_di:],self.dijt[:self.max_di+1]],kind='cubic')(self.us_di_sampler)
                self.pl[i_spl,i_det,j_det] = (self.dijt_adj.argmax() - self.nr_us) * self.max_di / self.nr_us

        self.pl[i_spl] += - self.pl[i_spl].T
        self.pl[i_spl] *= self.d_t
        
    return self


        
        
def fit_pl(self,min_vel=np.radians(.1),max_vel=np.radians(5),max_lag=5,method='hybrid'):

    #lag_flat = lambda dz, vx, vy : - (vx*np.real(dz) + vy*np.imag(dz)) / np.square(np.abs(vx+1j*vy) + 1e-4) 
    norm_lag_flat = lambda oa, vx, vy : (vx*np.cos(oa) + vy*np.sin(oa)) / np.square(np.abs(vx+1j*vy)) 

    self.vel_pars = np.zeros((*self.pl.shape[:2],2))        
    self.ss = {}; 

    pa_ba_ss_keys = ['tod_id','ba','pa']
    float_ss_keys = ['time_c', 'dt', 'azim_c', 'azim_v', 'da', 'elev_c', 'de', 'nt', 'lf_x_v', 'lf_y_v', 'co_r_2', 'prop_ok']
    array_ss_keys = ['time','azim','pl','cm']

    self.ss['dir'] = np.array(['none' for x in range(self.pl.shape[0])])
    self.ss['clust_x'] = self.clust_x
    self.ss['clust_y'] = self.clust_y

    for key in list(float_ss_keys):
        self.ss[key] = np.nan * np.zeros(self.pl.shape[:1],dtype=float)
    for key in list(array_ss_keys):
        self.ss[key] = np.nan * np.zeros(self.pl.shape[:1],dtype=object)

    DZ = np.subtract.outer(self.clust_x,self.clust_x) + 1j*np.subtract.outer(self.clust_y,self.clust_y)
    DN = np.outer(self.clust_n,self.clust_n) 
    bounds = [[-max_vel,-max_vel],[max_vel,max_vel]]

    for i_spl,(s,e) in enumerate(self.sub_splits):

        OA, OR = np.angle(DZ), np.abs(DZ)
        PL  = self.pl[i_spl]
        NPL = PL / (OR + 1e-16)
        NPL[np.isnan(NPL) | ~np.isfinite(NPL)] = 0
        USE = (np.abs(NPL) > 0) & (np.abs(PL) < max_lag) & (np.abs(NPL) < 1 / min_vel)

        prop_ok = USE.sum() / len(USE.ravel())
        lin_fit = lambda x,a,b : a*x+b

        ss_dt = 5e-2
        ss_time  = np.arange(self.time[s],self.time[e-1],ss_dt)
        ss_time += (self.time[e-1] - ss_time[-1]) / 2
        ss_nt = len(ss_time)

        ss_azim = sp.interpolate.interp1d(self.time[s:e],self.azim[s:e])(ss_time)
        ss_elev = sp.interpolate.interp1d(self.time[s:e],self.elev[s:e])(ss_time)
        ss_data = sp.interpolate.interp1d(self.time[s:e],self.clust_data[:,s:e].mean(axis=0))(ss_time)

        (da, azim_c), az_cpars = sp.optimize.curve_fit(lin_fit,ss_time-np.mean(ss_time),ss_azim,p0=[0,np.mean(ss_azim)])
        (de, elev_c), el_cpars = sp.optimize.curve_fit(lin_fit,ss_time-np.mean(ss_time),ss_elev,p0=[0,np.mean(ss_elev)])
        azim_v = (np.gradient(ss_azim) / np.gradient(ss_time)).mean()

        pl = PL[np.triu_indices(n=PL.shape[0],k=1)]

        lf_x_v, lf_y_v, r_2 = np.nan, np.nan, np.nan
        if prop_ok >= .5:

            vx0, vy0 = - np.mean(self.azv[s:e]*np.cos(self.elev[s:e])), 0
            (lf_x_v, lf_y_v), cpars = sp.optimize.curve_fit(norm_lag_flat,
                                                OA[USE],NPL[USE],
                                                p0=[vx0,vy0],sigma=np.abs(NPL*OR**2)[USE]**(-.5),bounds=bounds,maxfev=10000)

            TSS = np.square(NPL[USE]).sum()
            RSS = np.square(NPL[USE] - norm_lag_flat(OA[USE],lf_x_v, lf_y_v)).sum()
            co_r_2 = RSS / TSS

        else:
            co_r_2 = np.nan
            pass
        
        direction = 's'
        if np.degrees(azim_v) < -1e-1 : direction = 'l'
        if np.degrees(azim_v) > +1e-1 : direction = 'r'

        pb_ss = {'time_c':np.mean(ss_time), 'dt':ss_dt, 'azim_c':azim_c, 'azim_v':azim_v, 'da':da, 'elev_c':elev_c, 'de':de, 
                'nt':ss_nt, 'dir':direction, 'lf_x_v':lf_x_v, 'lf_y_v':lf_y_v, 'co_r_2':co_r_2, 'prop_ok':prop_ok,
                 'pl':pl, 'cm':ss_data, 'time':ss_time, 'azim':ss_azim}

        for key in list(pb_ss):
            self.ss[key][i_spl] = pb_ss[key]
            
    for key in list(array_ss_keys):
        self.ss[key] = list(self.ss[key])
        
    self.nss = len(self.ss['dir'])
            
    return self

def get_apex(self):
    
    import os, re, glob
    import pandas as pd
    apex_dict = {}
    for fn in glob.glob('/projects/ACT/mhasse/depots/actpol_shared/aux_data/apex_weather/targets/*'):
        
        #print(self.ss['time_c'].mean())
        #print(fn)
        key, s, e = re.findall('APEX_([a-z]+)_([0-9]+)_([0-9]+)',fn)[0]
        if int(s) < self.ss['time_c'].mean() < int(e):
            df = pd.read_csv(fn,index_col=0,header=None,delim_whitespace=True)
            ts_, val_ =  np.array(df.index), np.array(df[df.columns[0]])
            in_range  = (self.ss['time_c'].min() < ts_) & (self.ss['time_c'].max() > ts_)
            if not in_range.sum() > 2:
                val_ *= np.nan
            apex_dict[f'APEX_{key}'] = sp.interpolate.interp1d(ts_,val_,kind='cubic')(self.ss['time_c'])
            
    self.ss['apex_pwv'] = apex_dict[f'APEX_radiometer']
    wz = -apex_dict['APEX_windspeed']*np.exp(1j*(np.radians(apex_dict['APEX_winddirection'])-self.ss['azim_c']))
    self.ss['apex_wx'], self.ss['apex_wy'] = np.imag(wz), np.real(wz)
    self.ss['apex_temp'] = apex_dict['APEX_temperature']

    return self

def get_spectra(self):
    
    fmids = np.geomspace(5e-2,2e1,256)
    rfreq = np.exp(np.gradient(np.log(fmids)).mean())
    fbins = np.r_[fmids[0] / np.sqrt(rfreq), fmids * np.sqrt(rfreq)]

    self.ss['power'] = np.zeros((len(self.sub_splits),len(fmids)))
    self.ss['cm_power'] = np.zeros((len(self.sub_splits),len(fmids)))
    for i,(s,e) in enumerate(self.sub_splits):

        f,ps = sp.signal.periodogram(sp.signal.detrend(self.reg_data[:,s:e]),fs=1/self.d_t,window='hanning'); mps = ps.mean(axis=0)
        bmps = np.exp(sp.stats.binned_statistic(f,np.log(mps),bins=fbins)[0])
        self.ss['power'][i] = bmps
        
        f,cmps = sp.signal.periodogram(sp.signal.detrend(self.reg_data[:,s:e].mean(axis=0)),fs=1/self.d_t,window='hanning')
        bcmps  = np.exp(sp.stats.binned_statistic(f,np.log(cmps),bins=fbins)[0])
        self.ss['cm_power'][i] = bcmps
        
    return self