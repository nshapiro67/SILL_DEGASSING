import pandas as pd
import matplotlib.pyplot as plt
#------------------- mathematical finctions: https://numpy.org
import numpy as np
from obspy import read

import sys
#-----------------------------------------------



#-----------------------------------------------
# We use for solid:
# density=3.0   bulk mod.=80.e+3  shear mod.=50.e+3    visc=inf.
#        Vp [m/s]      =  6992.05896412151105
#         Vs [m/s]      =  4082.48290463863032

# Melt:
# dens.=2.8   bulk mod. =10.e+3  shear mod.=10.e+0   visc=1.0e-3 or visc=1.e+3 Pa.s

#         Vp [m/s]      =  1891.08182692328774
#         Vs [m/s]      =  59.7614304667196876
#-----------------------------------------------



#-----------------------------------------------
# function to smooth signal
#-----------------------------------------------
def smooth1d(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        print("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        print("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        print("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),x,mode='same')
    return y
#--------------------------------------------------------------


#-----------------------------------------------
# function to differentiate and smooth
#-----------------------------------------------
def smooth_differentiate(data, dt, nsmooth):
    ndata = np.size(data)
    dersm = np.zeros(ndata)
    dersm[0:ndata-1] = smooth1d(np.diff(data),nsmooth)/dt
    return dersm
#-----------------------------------------------


#-----------------------------------------------
# cosine right tapering
#-----------------------------------------------
def right_cos_taper(x,x1,x2):
    if x<x1:
        return 1.
    elif x>x2:
        return 0.
    else:
        return (1 + np.cos(2*np.pi*(x-x1)/(x2-x1)))/2
#-----------------------------------------------
    
#-----------------------------------------------
# cosine right tapering for array
#-----------------------------------------------
def arr_right_cos_taper(x,x1,x2):
    ct = x.copy()
    idl = x<x1
    idr = x>x2
    idm = (x>=x1) & (x<=x2)
    ct[idl] = 1
    ct[idr] = 0
    ct[idm] = (1 + np.cos(np.pi*(x[idm]-x1)/(x2-x1)))/2
    return ct
#-----------------------------------------------


#-----------------------------------------------
# function to compute Fourier spectra
#-----------------------------------------------
def signal_fft1d(sig,dt):
    npt = np.size(sig)
    spe = np.fft.fft(sig)
    freq = np.fft.fftfreq(npt,dt)
    sp_amp = np.sqrt(spe.real**2+spe.imag**2)
    sp_pha = np.arctan2(spe.imag, spe.real)
    npt_spe = int(npt/2)
    return npt_spe, dt*sp_amp[0:npt_spe],sp_pha[0:npt_spe],freq[0:npt_spe]
#--------------


#-----------------------------------------------
# function to compute Fourier spectra
#-----------------------------------------------
def signal_fft_version2(sig,dt):
    npt = np.size(sig)
    spa = np.absolute(np.fft.fft(sig))/np.sqrt(npt)
    df = 1/(npt*dt)
    freq = np.arange(npt)*df
    return spa, freq


#-----------------------------------------------
# function to read Melnik et al., 2020 files
#-----------------------------------------------
def read_data_res_spec(fname):
    data = pd.read_csv(fname+'res.csv')
    dd = data.values
    t = dd[:,0]
    r = dd[:,1]
    pg = dd[:,2]
    pm = dd[:,3]
    mco = dd[:,4]
    v = dd[:,7]
    
    data = pd.read_csv(fname+'spec.csv')
    dd = data.values
    f = dd[:,0]
    amp = dd[:,1]
    
    return t,r,pg,pm,mco,v,f,amp
#--------------------------------------





#----------------- closing all previous figureds
plt.close("all")

#---------------
#   crustal parameters for teh Green's function
#---------------
Vs = 3500
Vp = np.sqrt(3)*Vs
rho = 2900
dist = 40000

normS = (4*np.pi*rho*dist*Vs**3)
normP = (4*np.pi*rho*dist*Vp**3)

#---------------
#   elastic modulus in the source regions
#----------------
mu = 5e10
Kbulk = 8e10
lambd = (3*Kbulk-2*mu)/3.




filt_len = 0
#=======================!!!!!!!!!!!!!!!!!
# nmodel = 'N4D2'
# filt_len = 1

nmodel = 'N1A'
filt_len = 2

# nmodel = 'N3B'

# nmodel = 'N1R100'
# filt_len = 1



#---------------------
#  reading potency file
#----------------------
tmp = np.loadtxt('POTENCY/potency_'+nmodel+'.dat')

time = tmp[:,0]
delta_t = time[2]-time[1]


xx = tmp[:,1]
yy = tmp[:,2]
zz = tmp[:,3]
xy = tmp[:,4]
xz = tmp[:,5]
yz = tmp[:,6]



#---  STF
stf = zz

dstf = smooth_differentiate(stf, delta_t, 1000)
ddstf = smooth_differentiate(dstf, delta_t, 1000)

if filt_len == 1 :
    ddstf = ddstf*arr_right_cos_taper(time,0.6,0.65)
elif filt_len == 2 :
    ddstf = ddstf*arr_right_cos_taper(time,0.75,0.85)
else:
    ddstf = ddstf*arr_right_cos_taper(time,0.9,0.99)



#---------- Potency -> Moment
# coordinates:
# 0 - X
# 1 - Y
# 2 - Z

M0 =  np.zeros([3,3,np.size(zz)])

M0[0,0,:] = lambd*(xx+yy+zz) + 2*mu*xx
M0[1,1,:] = lambd*(xx+yy+zz) + 2*mu*yy
M0[2,2,:] = lambd*(xx+yy+zz) + 2*mu*zz
M0[0,1,:] = 2*mu*xy
M0[0,2,:] = 2*mu*xz
M0[1,2,:] = 2*mu*yz
M0[1,0,:] = M0[0,1,:]
M0[2,0,:] = M0[0,2,:]
M0[2,1,:] = M0[1,2,:]






M0scal = np.max(M0[0,0,:]+M0[1,1,:]+M0[2,2,:])

Mw = (2/3)*(np.log10(M0scal)-9.05)

print('moment magnitude ',Mw)



#---------------------
# moment rate

dM0 =  np.zeros([3,3,np.size(zz)])

dM0[0,0,:] = smooth_differentiate(M0[0,0,:], delta_t, 1000)
dM0[0,1,:] = smooth_differentiate(M0[0,1,:], delta_t, 1000)
dM0[0,2,:] = smooth_differentiate(M0[0,2,:], delta_t, 1000)
dM0[1,1,:] = smooth_differentiate(M0[1,1,:], delta_t, 1000)
dM0[1,2,:] = smooth_differentiate(M0[1,2,:], delta_t, 1000)
dM0[2,2,:] = smooth_differentiate(M0[2,2,:], delta_t, 1000)
dM0[1,0,:] = dM0[0,1,:]
dM0[2,0,:] = dM0[0,2,:]
dM0[2,1,:] = dM0[1,2,:]




#----- kroneker delta
kr_delta = np.zeros([3,3])
kr_delta[0,0] = 1
kr_delta[1,1] = 1
kr_delta[2,2] = 1


#------- ray angles (degrees)
theta = 30      # polar angle
phi = 90        # azimuth
theta_rad = theta*np.pi/180.
phi_rad = phi*np.pi/180.


#------ direction cosines
gamma = np.zeros(3)
gamma[0] = np.sin(theta_rad)*np.sin(phi_rad)
gamma[1] = np.sin(theta_rad)*np.cos(phi_rad)
gamma[2] = np.cos(theta_rad)


#------- displacements
up =  np.zeros([3,np.size(zz)])
us =  np.zeros([3,np.size(zz)])

for n in range (0,3):
    for p in range (0,3):
        for q in range (0,3):
            up[n] = up[n] + gamma[n]*gamma[p]*gamma[q]*dM0[p,q,:]/normP
            us[n] = us[n] + (gamma[n]*gamma[p]-kr_delta[n,p])*gamma[q]*dM0[p,q,:]/normS



#-------  velocities
vp =  np.zeros([3,np.size(zz)])
vs =  np.zeros([3,np.size(zz)])

for n in range (0,3):
    vp[n] = smooth_differentiate(up[n], delta_t, 1000)
    if filt_len == 1 :
        vp[n] = vp[n]*arr_right_cos_taper(time,0.6,0.65)
    elif filt_len == 2 :
        vp[n] = vp[n]*arr_right_cos_taper(time,0.75,0.85)
    else :
        vp[n] = vp[n]*arr_right_cos_taper(time,0.8,0.9)

    vs[n] = smooth_differentiate(us[n], delta_t, 1000)
    if filt_len == 1 :
        vs[n] = vs[n]*arr_right_cos_taper(time,0.6,0.65)
    elif filt_len == 2 :
        vs[n] = vs[n]*arr_right_cos_taper(time,0.75,0.85)
    else :
        vs[n] = vs[n]*arr_right_cos_taper(time,0.8,0.9)





nlsig = int(45/delta_t)

timel = np.arange(nlsig)*delta_t

vz = np.zeros(nlsig)

nfirstP = int(dist/Vp/delta_t)
nfirstS = int(dist/Vs/delta_t)

for i in range (0,np.size(zz)-2):
    vz[i+nfirstP] = vz[i+nfirstP] + vp[0,i]

for i in range (0,np.size(zz)-2):
    vz[i+nfirstS] = vz[i+nfirstS] + vs[0,i]



#--------------------------- real data# st = read('/Users/nshapiro/My Drive/MY_PAPERS/MELNIK_DEGASSING_SOURCE/MODEL_DATA/lgne.sac')
st = read('lgne.sac')


vzdata = st[0].data
delta_t_data = st[0].stats.delta
tdata = np.arange(st[0].stats.npts)*st[0].stats.delta



nspe, spamp2, sppha, fr2 = signal_fft1d(vz,delta_t)
nspe_data, spamp_data, sppha_data, fr_data = signal_fft1d(vzdata,delta_t_data)


spamp_s = smooth1d(spamp2,45)
spamp_data_s = smooth1d(spamp_data,45)



#---------------------------- plotting comparison
plt.figure(figsize=[9,6])

plt.subplots_adjust(left=None, bottom=0.07, right=None, top=.93, hspace=0.3, wspace=0.3)


ax1 = plt.subplot(221)
col='0.7'
ax1.plot(time,zz,color= col,linewidth=3,label = 'ZZ')
ax1.plot(time,xx,'--',linewidth=3,color= col,label = 'XX')
ax1.plot(time,yy,'--',linewidth=3,color= col,label = 'YY')
ax1.plot(time,xy,':',linewidth=3,color= col,label = 'XY')
ax1.plot(time,xz,':',linewidth=3,color= col,label = 'XZ')
ax1.plot(time,yz,':',linewidth=3,color= col,label = 'YZ')
ax1.set_ylabel('seismic potency ($\mathregular{m^3}$)',color=col)
ax1.tick_params(axis='y', labelcolor=col)
# plt.legend(loc='upper left',fontsize=8)
ax1.text(-0.1,1.05,'a)',weight='bold',fontsize=16,transform=ax1.transAxes)

plt.xlim(0,1)

ax12 = ax1.twinx()
col='k'
ax12.plot(time,M0[2,2,:],color=col,linewidth=1,label = 'ZZ')
ax12.plot(time,M0[0,0,:],'--',color=col,linewidth=1,label = 'XX')
ax12.plot(time,M0[1,1,:],'--',color=col,linewidth=1,label = 'YY')
ax12.plot(time,M0[0,1,:],':',color=col,linewidth=1,label = 'XY')
ax12.plot(time,M0[0,2,:],':',color=col,linewidth=1,label = 'XZ')
ax12.plot(time,M0[1,2,:],':',color=col,linewidth=1,label = 'YZ')
ax12.set_ylabel('seismic moment (Nm)')

plt.xlim(0,1)

plt.xlabel('time (s)')
ax12.legend(loc=[.85,.5],fontsize=6)



ax3 = plt.subplot(222)
ax3.plot(time,ddstf,'k')
plt.xlim(0,1)
plt.xlabel('time (s)')
plt.ylabel('potency 2-nd derivative ($\mathregular{m^3/s^2}$)')
ax3.text(-0.1,1.05,'b)',weight='bold',fontsize=16,transform=ax3.transAxes)
ax3.yaxis.set_ticks_position("right")
ax3.yaxis.set_label_position("right")



ax5 = plt.subplot(223)

ax5.plot(tdata,vzdata,'0.8',label = 'data')
ax5.plot(timel,vz*1e6,'k', label = 'model')

plt.xlabel('time (s)')
plt.ylabel('ground velocity (micron/s)')
plt.legend(loc='upper right',fontsize=8)
ax5.text(-0.1,1.05,'c)',weight='bold',fontsize=16,transform=ax5.transAxes)


ax6 = plt.subplot(224)

ax6.semilogx(fr_data,spamp_data_s/np.max(spamp_data_s), '0.8', label = 'data')
ax6.semilogx(fr2,spamp_s/np.max(spamp_s), 'k', label = 'model')
plt.xlim(0.1,100)
plt.xlabel('frequency (Hz)')
plt.ylabel('normalized Fourier amplitude')
plt.legend(loc='upper left',fontsize=8)
ax6.text(-0.1,1.05,'d)',weight='bold',fontsize=16,transform=ax6.transAxes)
ax6.yaxis.set_ticks_position("right")
ax6.yaxis.set_label_position("right")

                 
plt.show()
