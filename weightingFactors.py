import numpy as np
import matplotlib.pyplot as plt

from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel, Gaussian1DKernel
from astropy.convolution import convolve
import matplotlib.gridspec as gridspec
from astropy.cosmology import Planck15
from spectral_cube import SpectralCube
from spectral_cube import Projection
from astropy.io import fits
import matplotlib
import matplotlib.patheffects as patheffects

from matplotlib.ticker import NullFormatter, StrMethodFormatter


orange = '#ff9500'#(1,0.584,0)
blue =  '#007aff'  #(0,.478,1) blue
green = '#4cd964'
red = '#ff3b30'
grey = '#8e8e93'   #(142./255,142./255,147./255)


def Dlpen(redshift, giveAnswerInMeters = False):
    from numpy import sqrt
    redshift += 0.0
    Mpc = 3.0857e22
    Dl = Planck15.luminosity_distance(redshift).value
    if giveAnswerInMeters:
        return Dl*Mpc
    else:
        return Dl



def blackBody(v,T):
	h = 6.626e-34
	c = 3.0e8
	k = 1.38e-23
	from numpy import exp
	return (2*h*(v*v*v))/(c*c*exp((h*v)/(k*T)) - 1)

def Tdust(Tzero,z,Beta):
	# Dust temperature equivalent of a CMB-heated source
	return (((Tzero)**(Beta+4)) + ((2.73)**(Beta+4))*((1+z)**(Beta+4) -1))**(1/(4+Beta))

def funFintrinsic(Fobs,freq,Tzero,z,Beta):
	# Convert the observed flux into the flux as if the source were at z = 0
	bcmb = tomModelIso(freq,1,z,0.001*(1+z),Beta)
	bdust = tomModelIso(freq,1,z,Tdust(Tzero,z,Beta),Beta)
	return Fobs/(1 - (bcmb/bdust))


def funFobs(Fintrinsic,freq,Tzero,z,Beta):
	# Convert the intrinsic flux into the flux as it would be observed at the correct redshift
	bcmb = tomModelIso(freq,1,z,0.001*(1+z),Beta)
	bdust = tomModelIso(freq,1,z,Tdust(Tzero,z,Beta),Beta)
	return Fintrinsic*(1 - (bcmb/bdust))



def giveKappa(freq,fstar=1900e9, kappastar=10.41,Beta=2.03):
	return kappastar*((freq/fstar)**Beta)


def Fv(freq_obs,z,logMdust,Tzero,beta):
	freq = freq_obs*(1+z)
	Mdust = 10**logMdust
	kappa = giveKappa(freq,Beta=beta)  ## cm2/g
	dl = Dlpen(z,giveAnswerInMeters=True) ## m
	TdustVal = Tdust(Tzero,z,beta) ## K
	g = (1+z)/(dl**2) ## 1/m2
	Tcmb = (1+z)*2.7255
	BBval = (blackBody(freq,TdustVal)-blackBody(freq,Tcmb)) ## W / sr / m2 / Hz
	convUnit = (1.0e-4)*(1.988e33)*(1e26) ## cm2-m2 ; g-Mstar ; Jy - W 
	return convUnit*g*Mdust*kappa*BBval


sourceName = [
'REBELS-24',
'GS_z9_3',
'ID4590 Band 7',
'ID4590 Band 5',
'GS-z10-1',
'GHZ1',
'JD1 Band 5',
'JD1 Band 7',
'SPT0615-JD',
'GN-z11',
'COS-z12-1',
'GHZ2 Band 6',
'GHZ2 Band 8',
'GS-z14 Band 7',
'GS-z14 Band 5',
'GS-z11 Band 7'
]



fitsfiles = [
'REBELS24cont.fits',
'GS_z9_3_cont.fits',
'ID4590_cont2.fits',
'ID4590_cont1.fits',
'GS-z10-1_cont.image.fits',
'GHZ1_GLz11cont.fits',
'JD1_cont_cii.fits',
'JD1_cont_new.fits',
'SPT0615-JDcont.fits',
'GNz11.fits',
'COS-z12-1mfs.fits',
'GHZ2_mfs_new.fits',
'GHZ2z12_b8_cont.fits',
'gs14_mfs-complete.fits',
'GSz14_targets_cont158.image.fits',
'GS-z11.fits'
]




stellarMass = np.array([
8.89,# 'REBELS24cont.fits',
9.19,# 'GS_z9_3_cont.fits',
7.78,# 'ID4590_cont2.fits',
7.78,# 'ID4590_cont1.fits',
8.17,# 'GS-z10-1_cont.image.fits',
9.2,# 'GHZ1_GLz11cont.fits',
8.2,# 'JD1_cont_cii.fits',
8.2,# 'JD1_cont_new.fits',
7.47,# 'SPT0615-JDcont.fits',
8.73,# 'GNz11.fits',
9.60,# 'COS-z12-1mfs.fits',
9.05,# 'GHZ2_mfs_new.fits',
9.05,# 'GHZ2z12_b8_cont.fits',
8.7,# 'gs14_mfs-complete.fits',
8.7,# 'GSz14_targets_cont158.image.fits'
8.67#GSz11
	])

redshiftsValues = np.array([8.21,#REBELS24cont.fits,
8.228,#GS_z9_3_cont.fits,
8.496,#ID4590_cont2.fits,
8.496,#ID4590_cont1.fits,
9.433,#GS-z10-1_cont.image.fits,
9.875,#GHZ1_GLz11cont.fits,
9.110,#JD1_cont_cii.fits,
9.110,#JD1_cont_new.fits,
10.16,#SPT0615-JDcont.fits,
10.603,
12.25,#COS-z12-1mfs.fits,
12.333,#GHZ2_mfs_new.fits,
12.333,#GHZ2z12_b8_cont.fits,
14.178,
14.178,
11.58
])

magnificationFactor = [
1,# 'REBELS24cont.fits'
1,# 'GS_z9_3_cont.fits'
8.69,# 'ID4590_cont2.fits'
8.69,# 'ID4590_cont1.fits'
1,# 'GS-z10-1_cont.image.fits'
1,# 'GHZ1_GLz11cont.fits'
11.5,# 'JD1_cont_cii.fits'
11.5,# 'JD1_cont_new.fits'
120,# 'SPT0615-JDcont.fits'
1,# 'GNz11.fits'
1,# 'COS-z12-1mfs.fits'
1.3,# 'GHZ2_mfs_new.fits'
1.3,# 'GHZ2z12_b8_cont.fits'
1,# 'gs14_mfs-complete.fits'
1,# 'GSz14_targets_cont158.image.fits'
1
]

freq = np.array([202,# 'REBELS24cont.fits'
349,# 'GS_z9_3_cont.fits'
360,# 'ID4590_cont2.fits'
200,# 'ID4590_cont1.fits'
307,# 'GS-z10-1_cont.image.fits'
248,# 'GHZ1_GLz11cont.fits'
335,# 'JD1_cont_cii.fits'
188,# 'JD1_cont_new.fits'
307,# 'SPT0615-JDcont.fits'
161,# 'GNz11.fits'
258,# 'COS-z12-1mfs.fits'
248,# 'GHZ2_mfs_new.fits'
427,# 'GHZ2z12_b8_cont.fits'
224,# 'gs14_mfs-complete.fits'
125,# 'GSz14_targets_cont158.image.fits'
281])#GSz11


sourceMask = np.array([
1,# 'REBELS24cont.fits'
1,# 'GS_z9_3_cont.fits'
1,# 'ID4590_cont2.fits'
1,# 'ID4590_cont1.fits'
1,# 'GS-z10-1_cont.image.fits'
1,# 'GHZ1_GLz11cont.fits'
1,# 'JD1_cont_cii.fits'
1,# 'JD1_cont_new.fits'
1,# 'SPT0615-JDcont.fits'
1,# 'GNz11.fits'
1,# 'COS-z12-1mfs.fits'
1,# 'GHZ2_mfs_new.fits'
1,# 'GHZ2z12_b8_cont.fits'
1,# 'gs14_mfs-complete.fits'
1,# 'GSz14_targets_cont158.image.fits'
1#GSz11
	])

radius = np.array([
100,# 'REBELS24cont.fits'
100,# 'GS_z9_3_cont.fits'
280,# 'ID4590_cont2.fits'
280,# 'ID4590_cont1.fits'
110,# 'GS-z10-1_cont.image.fits'
500,# 'GHZ1_GLz11cont.fits'
332,# 'JD1_cont_cii.fits'
332,# 'JD1_cont_new.fits'
100,# 'SPT0615-JDcont.fits'
64,# 'GNz11.fits'
420,# 'COS-z12-1mfs.fits'
105,# 'GHZ2_mfs_new.fits'
105,# 'GHZ2z12_b8_cont.fits'
260,# 'gs14_mfs-complete.fits'
260,# 'GSz14_targets_cont158.image.fits'
80#GSz11
	])


betaUV = np.array([
-1.5,# 'REBELS24cont.fits'
-2,# 'GS_z9_3_cont.fits'
-1.7,# 'ID4590_cont2.fits'
-1.7,# 'ID4590_cont1.fits'
-2.54,# 'GS-z10-1_cont.image.fits'
-1.79,# 'GHZ1_GLz11cont.fits'
-2.2,# 'JD1_cont_cii.fits'
-2.2,# 'JD1_cont_new.fits'
-2.7,# 'SPT0615-JDcont.fits'
-2.36,# 'GNz11.fits'
-1.78,# 'COS-z12-1mfs.fits'
-2.46,# 'GHZ2_mfs_new.fits'
-2.46,# 'GHZ2z12_b8_cont.fits'
-2.2,# 'gs14_mfs-complete.fits'
-2.2,# 'GSz14_targets_cont158.image.fits'
-2.18#GSz11
	])

def giveOpticalMdust(beta,radius):
	betazero = 4.43
	betazero = 2.63*1.99
	tau = betazero +1.99*beta
	mval = (4/3)*(1/(2.17e-8))*tau*((radius/1000)**2)
	return mval


def giveOpticalMdustOld(beta,radius):
	betazero = 4.43
	betazero = 2.23*1.99
	tau = betazero +1.99*beta
	mval = (4/3)*(1/(2.17e-8))*tau*((radius/1000)**2)
	return mval


TdustVal=50
BetaValAssumed = 1.5
sensitivityLimit = []
dustMassLimit = []
dustToStellarMassLimit = []



for i in range(len(fitsfiles)):
	data = fits.open('cont/'+fitsfiles[i])
	datalen = data[0].data.shape[0]
	val1 = -1*np.percentile(data[0].data[np.isnan(data[0].data) == False],15.865)#
	val2 = np.percentile(data[0].data[np.isnan(data[0].data) == False],84.135)#
	val3 = np.nanstd(data[0].data[int(datalen/2 - 25):int(datalen/2 + 25),int(datalen/2 - 25):int(datalen/2 + 25)])
	stdtemp = min(val1,val2,val3)
	sensitivityLimit.append(stdtemp)
	mdust = (1/magnificationFactor[i])*1/Fv(freq[i]*1e9,redshiftsValues[i],0,TdustVal,BetaValAssumed)
	dustMassLimit.append(mdust)
	print(fitsfiles[i])
	mopt = (1/magnificationFactor[i])*giveOpticalMdust(betaUV[i],radius[i])
	moptold = (1/magnificationFactor[i])*giveOpticalMdustOld(betaUV[i],radius[i])
	print(mopt*1e-6)
	print(moptold*1e-6)
	if stellarMass[i] == -99:
		dustToStellarMassLimit.append(-99)
	else:
		dustToStellarMassLimit.append(mdust/(10**stellarMass[i]))
	try:
		BMAJ = str(np.round(data[0].header['BMAJ']*3600,2))
		BMIN = str(np.round(data[0].header['BMIN']*3600,2))
		BPA = str(np.round(data[0].header['BPA'],1))
		OBJECT = str((data[0].header['OBJECT'],1))
	except:
		1+1


for i in range(len(fitsfiles)):
	try:
		contFile = 'cont/'+fitsfiles[i]
		cont = SpectralCube.read(contFile)
		cont = cont[0]
	except:
		contFile = 'cont/'+fitsfiles[i]
		hdul = fits.open(contFile)
		cont = Projection.from_hdu(hdul)
	target_header = cont.header
	if abs(0.1/3600 - target_header['CDELT1']) > 0.01/3600:
		target_header['CDELT1'] = 0.1/3600
		target_header['CDELT2'] = 0.1/3600
		resampledCont = cont.reproject(target_header)
		resampledCont.write('reprojectCont/'+fitsfiles[i],overwrite=True)




for i in range(len(fitsfiles)):
	try:
		contFile = 'reprojectCont/'+fitsfiles[i]
		cont = SpectralCube.read(contFile)
		cont = cont[0]
	except:
		contFile = 'reprojectCont/'+fitsfiles[i]
		hdul = fits.open(contFile)
		cont = Projection.from_hdu(hdul)
	contcopy = cont.copy()
	dustMass = contcopy*dustMassLimit[i]
	dustMass.write('dustMassFits/'+fitsfiles[i],overwrite=True)
	if dustToStellarMassLimit[i] > 0:
		contcopy = cont.copy()
		dustMass = contcopy*dustToStellarMassLimit[i]
		dustMass.write('dustToStellarMassFits/'+fitsfiles[i],overwrite=True)
































