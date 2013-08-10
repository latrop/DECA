#!/usr/bin/python
# -*- coding:  cp1251 -*-

import sys
import math
import numpy as np
from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import NullFormatter
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
from numpy import fft
import pyfits
import re

DECA_PATH = os.path.split(os.path.dirname(__file__))[0]

sys.path.append(DECA_PATH)
sys.path.append(DECA_PATH + '/ini_modules')
sys.path.append(DECA_PATH + '/analysis_modules')
sys.path.append(DECA_PATH + '/output_modules')
sys.path.append(DECA_PATH + '/prep_modules')

import setup
import sextr
import image
import sky_est

SN_limit = setup.SN_limit
star_cl = setup.star_cl
gal_cl = setup.gal_cl
gal_infl_coef = setup.gal_infl_coef
bright_coeff = setup.bright_coeff
kron_coef = setup.kron_coef

box_psf = setup.box_psf
beta = setup.beta
window = setup.window

tmp_out = sys.stdout

#**************** Analytical Functions*****************
def Imof(r,fwhm,I0,beta):
#Moffat function
	alpha = fwhm/(2.*sqrt(2**(1./beta) -1.))
	return I0*(1.+(fabs(r)/alpha)**2)**(-beta)

def Igau(r,fwhm,I0):
	c = fwhm/2.35482
	return I0*np.exp(-fabs(r)**2/(2.*c*c))


#*********** Convolution***********
def convolution(f,h):
	if type(f)!=ndarray: f = np.array(f)
	if type(h)!=ndarray: h = np.array(h)
	g = fft.irfft( fft.rfft(f)*fft.rfft(h) )
	return g

#*********** Deconvolution*********
def deconvolution(g,h):
	if type(g)!=ndarray: g = np.array(g)
	if type(h)!=ndarray: f = np.array(f)
	f = fft.irfft( fft.rfft(g)/fft.rfft(h) )
	return f


def deconv_gauss(inten,r,fwhm=1.3):
	m0 = 26.
	mstar = 12.
	lstar= 10**(0.4*(m0-mstar))
	c = fwhm/2.35482
	I0 = lstar/(c*sqrt(2.0*math.pi))
	Ideconv = deconvolution(inten,Igau(r,fwhm,I0)/Igau(r,fwhm,I0).sum())
	return Ideconv

def deconv_moffat(inten,r,fwhm=1.3,beta=4.765):
	m0 = 26.
	mstar = 2.
	lstar= 10**(0.4*(m0-mstar))
	c = fwhm/2.35482
	I0 = lstar/(c*sqrt(2.0*math.pi))
	Ideconv = deconvolution(inten,Imof(r,fwhm,I0,beta)/Imof(r,fwhm,I0,beta).sum())
	return Ideconv



# Function to create a psf-image with GALFIT

# ********************************* Simulation with GALFIT *********************************
def mod_psf(fwhm_psf,M_psf,ell,PA,m0,pix2secx,pix2secy):
	nx_psf = box_psf
	ny_psf = box_psf
	if nx_psf%2==0:
		xc_psf = int(nx_psf/2. + 1)
	else:
		xc_psf = int(nx_psf/2. + 0.5)
	if ny_psf%2==0:
		yc_psf = int(ny_psf/2. + 1)
	else:
		yc_psf = int(ny_psf/2. + 0.5)

	#========================================
	#pix2secx = pix2sec
	#pix2secy = pix2sec
	#========================================


	f = open(r"modelPSF.txt", "w") 
	sys.stdout = f
	print "\n==============================================================================="
	print "# IMAGE and GALFIT CONTROL PARAMETERS"
	print "A) none                # Input data image (FITS file)"
	print "B) psf.fits         # Output data image block"
	print "C) none                # Sigma image name (made from data if blank or none)" 
	print "D) none                # Input PSF image and (optional) diffusion kernel"
	print "E) 1                   # PSF fine sampling factor relative to data" 
	print "F) none                # Bad pixel mask (FITS image or ASCII coord list)"
	print "G) none                # File with parameter constraints (ASCII file)" 
	print "H) 1    %i   1    %i   # Image region to fit (xmin xmax ymin ymax)" % (nx_psf, ny_psf)
	print "I) %.3f    %.3f        # Size of the convolution box (x y)" % (0, 0)
	print "J) %.3f              # Magnitude photometric zeropoint" % (m0)
	print "K) %.3f    %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (pix2secx,pix2secy)
	print "O) regular             # Display type (regular, curses, both)"
	print "P) 1                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n"

	print "# INITIAL FITTING PARAMETERS"
	print "#"
	print "#   For object type, the allowed functions are:" 
	print "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat," 
	print "#       ferrer, powsersic, sky, and isophote." 
	print "#"  
	print "#   Hidden parameters will only appear when they're specified:"
	print "#       C0 (diskyness/boxyness)," 
	print "#       Fn (n=integer, Azimuthal Fourier Modes),"
	print "#       R0-R10 (PA rotation, for creating spiral structures)."
	print "#" 
	print "# -----------------------------------------------------------------------------"
	print "#   par)    par value(s)    fit toggle(s)    # parameter description" 
	print "# -----------------------------------------------------------------------------\n"

	if window=='gauss':
		print "# Gaussian function\n"
		print "0) gaussian           # object type"
		print "1) %.3f  %.3f  0 0   # position x, y        [pixel]" % (xc_psf,yc_psf)
		print "3) %.3f       0        # total magnitude" % (M_psf)     
		print "4) %.3f       0        #   FWHM               [pixels]" % (fwhm_psf)
		print "9) %.3f        0       # axis ratio (b/a)" % (1.-ell)  
		print "10) %.3f         0       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA)
		print "Z) 0                  # leave in [1] or subtract [0] this comp from data"
		print "\n================================================================================"

	if window=='moffat':
		print "# Moffat function\n"
		print "0) moffat           # object type"
		print "1) %.3f  %.3f  0 0   # position x, y        [pixel]" % (xc_psf,yc_psf)
		print "3) %.3f       0        # total magnitude" % (M_psf)     
		print "4) %.3f       0        #   FWHM               [pixels]" % (fwhm_psf)
		print "5) %.3f        0       # powerlaw" % (beta)
		print "9) %.3f        0       # axis ratio (b/a)" % (1.-ell)  
		print "10) %.3f         0       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA)
		print "Z) 0                  # leave in [1] or subtract [0] this comp from data"
		print "\n================================================================================"          

	sys.stdout = tmp_out
	f.close()
	os.chmod(r"modelPSF.txt",0777)
	subprocess.call("galfit modelPSF.txt", shell=True)

def SDSS_psf(file_in,xc,yc):
	coords = str(xc)+".0 "+str(yc)+".0"
	#psf = "read_PSF psField.fit 2 "+ str(xc)+str(yc)+" psf.fit"
	psf = "./read_PSF psField.fit 2 "+ coords + " psf.fit"
	#subprocess.call("read_PSF psField-001336-2-0051.fit 5 500.0 600.0 foo.fit", shell=True)
	subprocess.call(psf, shell=True)
	#try:
	#	from pyraf import iraf
	#	iraf.imarith(operand1="psf.fit",op="-",operand2="1000.0",result="psf1.fit")
	#except:
	sky_est.imarith_func("psf.fit","-",1000.,"psf1.fit")
	exit()
	os.remove(r"psf.fit")
	os.rename(r"psf1.fit","psf.fits")


def gal_clean_ima(file_segm,file_gal,file_gal_clean,xc,yc):
		hdulist1 = pyfits.open(file_segm, do_not_scale_image_data=True)
		img1 = hdulist1[0].data
		hdulist2 = pyfits.open(file_gal, do_not_scale_image_data=True)
		img2 = hdulist2[0].data

		(dimy,dimx) = img1.shape
		shutil.copy(file_gal,file_gal_clean) 
		hdulist3 = pyfits.open(file_gal_clean, do_not_scale_image_data=True, mode='update')
		img3 = hdulist3[0].data

		for i in range(dimy):
			for k in range(dimx):
				if img1[i,k]>0 and img1[i,k]!=img1[yc,xc]:
					img3[i,k] = 0.
		hdulist3.flush()

def noise(image,sky_level):
	#iraf.imarith(image,'-',sky_level,"field1.fits")
	#try:
	#	from pyraf import iraf
	#	iraf.imarith(operand1=image,op="-",operand2=sky_level,result="field1.fits")
	#except:
	sky_est.imarith_func(image,"-",sky_level,"field1.fits")
	hdulist1 = pyfits.open("fon.fits", do_not_scale_image_data=True)
	img1 = hdulist1[0].data
	fon = fabs(np.std(img1))
	return fon

def StarPSF_extr(image,xc,yc,sky_level):
	box_psf1 = box_psf - 1
	if box_psf1%2==0:	
		xmin = int((xc -1) - int(box_psf1/2))
		xmax = int((xc -1) + int(box_psf1/2))
		ymin = int((yc -1) - int(box_psf1/2))
		ymax = int((yc -1) + int(box_psf1/2))
	else:
		xmin = int((xc -0.5) - box_psf1/2.)
		xmax = int((xc -0.5) + box_psf1/2.)
		ymin = int((yc -0.5) - box_psf1/2.)
		ymax = int((yc -0.5) + box_psf1/2.)

	gal_clean_ima("segm.fits","field2.fits","field_clean.fits",xc,yc)		# earlier was field1.fits

	#file_in = "field_clean.fits" + '['+str(xmin)+':'+str(xmax)+','+str(ymin)+':'+str(ymax)+']'
	#iraf.imcopy(file_in,"psf.fits")
	if os.path.exists('psf.fits')==True:
		os.remove('psf.fits')
	sky_est.imcopy_func("field_clean.fits","psf.fits",int(xmin),int(ymin),int(xmax),int(ymax))
	#exit()


def signal_to_noise(m_star,R,image,sky_level,GAIN,m0):
	#http://iraf.net/phpBB2/viewtopic.php?t=87000
	F_star = 10**(0.4*(m0-m_star))
	N = math.pi*R*R
	F_backgr = fabs(noise(image,sky_level))
	SN = F_star*GAIN/(sqrt( F_star*GAIN + N*F_backgr*GAIN + N*R*R )) 
	return SN





def GoodStars_find(image_file,file_out,GAIN,m0,sky_level):
	#number,x_image,y_image,flux_radius,mag_auto,a_image,b_image,theta_image,ellipticity,kron_radius,backgr,class_star = loadtxt(file_out, usecols=[0,5,6,9,12,13,14,15,16,17,18,19], unpack=True, skiprows = 18)
	number,x_image,y_image,x_world,y_world,XPEAK_IMAGE,YPEAK_IMAGE,XPEAK_WORLD,YPEAK_WORLD,flux_radius25,flux_radius50,flux_radius75,flux_radius99,mag_auto,a_image,b_image,theta_image,ellipticity,kron_radius,backgr,class_star = loadtxt(file_out, usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], unpack=True, skiprows = 18)
	print class_star


	mag = np.median(mag_auto) - bright_coeff
	xc = []
	yc = []
	r = []
	for k in range(0,len(number)):
		if class_star[k] < gal_cl:
			xc.append(x_image[k])
			yc.append(y_image[k])
			r.append(kron_radius[k]*max(a_image[k],b_image[k])*gal_infl_coef)
	print "Number of galaxies is %i" % (len(r))



	xc_star = []
	yc_star = []
	re_star = []
	m_star = []
	PA_star = []
	ell_star = []
	backgr_level = []
	StoN = []
	R = []

	R_kron = a_image*kron_radius
	SN = signal_to_noise(mag_auto,R_kron,image_file,sky_level,GAIN,m0)
	print max(SN)
	ii = []

	for k in range(0,len(number)):
		if class_star[k] > star_cl and SN[k]>SN_limit:
			ii.append(k)


	for k in range(len(ii)):
		kk = 0
		for l in range(len(r)):
 			if sqrt( (x_image[ii[k]] - xc[l])**2 + (y_image[ii[k]] - yc[l])**2 ) > r[l]:
				kk = kk+1
				if kk==len(r):
					xc_star.append(x_image[ii[k]])
					yc_star.append(y_image[ii[k]])
					re_star.append(flux_radius50[ii[k]])
					m_star.append(mag_auto[ii[k]])
					PA_star.append(theta_image[ii[k]])
					ell_star.append(ellipticity[ii[k]])
					backgr_level.append(backgr[ii[k]])
					R.append(R_kron[ii[k]])
					StoN.append(SN[ii[k]])

	print "Number of good PSF stars is %i" % (len(re_star))
	return xc_star,yc_star,re_star,m_star,PA_star,ell_star,backgr_level,R,StoN,len(re_star)


def header(file_image,file_out,bad_pix_mask,xc,yc,m0,pix2secx,pix2secy):
	print "\n==============================================================================="
	print "# IMAGE and GALFIT CONTROL PARAMETERS"
	print "A) %s                # Input data image (FITS file)" % (file_image)
	print "B) %s  	              # Output data image block" % (file_out)
	print "C) none                # Sigma image name (made from data if blank or none)" 
	print "D) none                # Input PSF image and (optional) diffusion kernel"
	print "E) 1                   # PSF fine sampling factor relative to data" 
	print "F) %s                # Bad pixel mask (FITS image or ASCII coord list)" % (bad_pix_mask)
	print "G) none                # File with parameter constraints (ASCII file)" 
	print "H) %i    %i   %i    %i   # Image region to fit (xmin xmax ymin ymax)" % (xc-int(floor(box_psf/2.)),xc+int(floor(box_psf/2.)),yc-int(floor(box_psf/2.)),yc+int(floor(box_psf/2.)))
	print "I) 0    0        # Size of the convolution box (x y)"
	print "J) %.3f              # Magnitude photometric zeropoint" % (m0)
	print "K) %.3f    %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (pix2secx,pix2secy)
	print "O) regular             # Display type (regular, curses, both)"
	print "P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n"

	print "# INITIAL FITTING PARAMETERS"
	print "#"
	print "#   For object type, the allowed functions are:" 
	print "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat," 
	print "#       ferrer, powsersic, sky, and isophote." 
	print "#"  
	print "#   Hidden parameters will only appear when they're specified:"
	print "#       C0 (diskyness/boxyness)," 
	print "#       Fn (n=integer, Azimuthal Fourier Modes),"
	print "#       R0-R10 (PA rotation, for creating spiral structures)."
	print "#" 
	print "# -----------------------------------------------------------------------------"
	print "#   par)    par value(s)    fit toggle(s)    # parameter description" 
	print "# -----------------------------------------------------------------------------\n"

def moffat (xc,yc,mtot,fwhm,ell,PA):
	 print "Moffat function"
 	 print "0) moffat           # object type"
	 print "1) %.3f  %.3f  1 1  # position x, y        [pixel]" % (xc,yc)
	 print "3) %.3f       1       # total magnitude" % (mtot)     
	 print "4) %.3f       1       #   FWHM               [pixels]" % (fwhm)
	 print "5) %.3f        1       # powerlaw" % (beta) 
	 print "9) %.3f      1       # axis ratio (b/a)" % (1.-ell)
	 print "10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA)
	 print "Z) 1                  # leave in [1] or subtract [0] this comp from data?"

def gauss (xc,yc,mtot,fwhm,ell,PA):
	 print "Gaussian function"
 	 print "0) gaussian           # object type"
	 print "1) %.3f  %.3f  1 1  # position x, y        [pixel]" % (xc,yc)
	 print "3) %.3f       1       # total magnitude" % (mtot)     
	 print "4) %.3f       1       #   FWHM               [pixels]" % (fwhm)
	 print "9) %.3f      1       # axis ratio (b/a)" % (1.-ell)
	 print "10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA)
	 print "Z) 1                  # leave in [1] or subtract [0] this comp from data?"

def sky(sky_level,dskyX,dskyY):
	 print "# sky\n"
 	 print "0) sky                # Object type"
 	 print "1) %.3f       1       # sky background       [ADU counts]" % (sky_level)
 	 print "2) %.3f      0       # dsky/dx (sky gradient in x)" % (dskyX)
	 print "3) %.3f      0       # dsky/dy (sky gradient in y)" % (dskyY)
	 print "Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"

def PSF_model_sph(survey,image_file,m0,EXPTIME,pix2secX,pix2secY,GAIN,sky_level):
	xc_star1,yc_star1,re_star1,m_star1,PA_star1,ell_star1,backgr_level1,R1,SN1,Nstars = GoodStars_find(image_file,'field.cat',GAIN,m0,sky_level)
	PA_star1 = image.pa_arr(PA_star1)
	xc_star = []; yc_star = []; re_star = []; m_star = []; PA_star = []; ell_star = []; backgr_level = []; R = []; SN = []


	if len(xc_star1)>10:
		s = random.sample(range(len(xc_star1)), 10)
		for k in range(len(s)):
			xc_star.append(xc_star1[s[k]])
			yc_star.append(yc_star1[s[k]])
			re_star.append(re_star1[s[k]])
			m_star.append(m_star1[s[k]])
			PA_star.append(PA_star1[s[k]])
			ell_star.append(ell_star1[s[k]])
			backgr_level.append(backgr_level1[s[k]])
			R.append(R1[s[k]])
			SN.append(SN1[s[k]])


	hdulist = pyfits.open(image_file)
	scidata = hdulist[0].data
	ny,nx = scidata.shape
	mtot_star=[]
	fwhm_star=[]
	beta_star=[]
	theta_star=[]
	Ell_star=[]
	chi2_star=[]


	for k in range(len(xc_star)):
		if (xc_star[k]-int(box_psf/2.))>0. and (xc_star[k]+int(box_psf/2.))<nx and (yc_star[k]-int(box_psf/2.))>0. and (yc_star[k]+int(box_psf/2.))<ny:
			sextr.bad_pix_mask("segm.fits","bad_pix_stars.txt",xc_star[k],yc_star[k])
			f = open(r"modelIN.txt", "w") 
			sys.stdout = f
			header(image_file,'psf.fits',"bad_pix_stars.txt",xc_star[k],yc_star[k],m0 - 2.5*log10(EXPTIME),pix2secX,pix2secY)
			if window=='moffat':	moffat (xc_star[k],yc_star[k],m_star[k],3.3,ell_star[k],PA_star[k])
			if window=='gauss':	gauss (xc_star[k],yc_star[k],m_star[k],3.3,ell_star[k],PA_star[k])
			#sky(backgr_level[k],0.0,0.0)
			sys.stdout = tmp_out
			f.close()
			os.chmod(r"modelIN.txt",0777)
			subprocess.call("galfit modelIN.txt", shell=True)
			hdulist = pyfits.open('psf.fits')
			prihdr = hdulist[2].header
			mtot_star.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_MAG']))[0])
			fwhm_star.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_FWHM']))[0])
			if window=='moffat':	beta_star.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_C']))[0])
			theta_star.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_PA']))[0])
			Ell_star.append(1.-map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_AR']))[0])
			chi2_star.append(prihdr['CHI2NU'])
			os.remove(r"modelIN.txt") 
			os.remove(r"fit.log") 
			os.remove(r"galfit.01") 
			os.remove(r"bad_pix_stars.txt")
 			os.remove(r"psf.fits")
	theta = np.median(theta_star)
	ell = np.median(Ell_star)
	fwhm = np.median(fwhm_star)
	mtot = np.median(mtot_star)
	if window=='moffat':	beta = np.median(beta_star)
	mod_psf(fwhm,mtot,ell,theta,m0 - 2.5*log10(EXPTIME),pix2secX,pix2secY)
	return fwhm,Nstars



def PSF_model_best(survey,image_file,m0,EXPTIME,pix2secX,pix2secY,GAIN,sky_level):
	xc_star,yc_star,re_star,m_star,PA_star,ell_star,backgr_level,R,SN,Nstars = GoodStars_find(image_file,'field.cat',GAIN,m0,sky_level)
	PA_star = image.pa_arr(PA_star)
	ii = []
	SNN = []
	hdulist = pyfits.open(image_file)

	scidata = hdulist[0].data


	ny,nx = scidata.shape
	for k in range(len(xc_star)):
		if (xc_star[k]-int(box_psf/2.))>0. and (xc_star[k]+int(box_psf/2.))<nx and (yc_star[k]-int(box_psf/2.))>0. and (yc_star[k]+int(box_psf/2.))<ny:
			ii.append(k)
			SNN.append(SN[k])
	
	iii = ii[SNN.index(max(SNN))]
	sextr.bad_pix_mask("segm.fits","bad_pix_stars.txt",xc_star[iii],yc_star[iii])
	f = open(r"modelIN.txt", "w") 
	sys.stdout = f
	header(image_file,'psf.fits',"bad_pix_stars.txt",xc_star[iii],yc_star[iii],m0 - 2.5*log10(EXPTIME),pix2secX,pix2secY)
	if window=='moffat':	moffat (xc_star[iii],yc_star[iii],m_star[iii],3.3,ell_star[iii],PA_star[iii])
	if window=='gauss':	gauss (xc_star[iii],yc_star[iii],m_star[iii],3.3,ell_star[iii],PA_star[iii])
	#sky(backgr_level[iii],0.0,0.0)
	sys.stdout = tmp_out
	f.close()
	os.chmod(r"modelIN.txt",0777)
	subprocess.call("galfit modelIN.txt", shell=True)
	hdulist = pyfits.open('psf.fits')
	prihdr = hdulist[2].header
	mtot_star = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_MAG']))[0]
	fwhm_star = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_FWHM']))[0]
	if window=='moffat':	beta_star = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_C']))[0]
	theta_star = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_PA']))[0]
	Ell_star = 1.-map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_AR']))[0]
	chi2_star = prihdr['CHI2NU']
	os.remove(r"modelIN.txt") 
	os.remove(r"fit.log") 
	os.remove(r"galfit.01") 
	os.remove(r"bad_pix_stars.txt")
 	os.remove(r"psf.fits")
	mod_psf(fwhm_star,mtot_star,Ell_star,theta_star,m0 - 2.5*log10(EXPTIME),pix2secX,pix2secY)
	return fwhm_star*sqrt(pix2secX*pix2secY),Nstars
def PSF_star_best(survey,image_file,m0,pix2secX,pix2secY,GAIN,sky_level,EXPTIME):
	#iraf.imarith(image_file,'-',sky_level,"field1.fits")
	hdulist = pyfits.open(image_file)
	scidata = hdulist[0].data
	ny,nx = scidata.shape
	xc_star,yc_star,re_star,m_star,PA_star,ell_star,backgr_level,R,SN,Nstars = GoodStars_find(image_file,'field.cat',GAIN,m0,sky_level)
	PA_star = image.pa_arr(PA_star)
	ii = []
	SNN = []
	for k in range(len(xc_star)):
		if (xc_star[k]-int(box_psf/2.))>0. and (xc_star[k]+int(box_psf/2.))<nx and (yc_star[k]-int(box_psf/2.))>0. and (yc_star[k]+int(box_psf/2.))<ny:
			ii.append(k)
			SNN.append(SN[k])
	
	iii = ii[SNN.index(max(SNN))]
	StarPSF_extr(image_file,xc_star[iii],yc_star[iii],sky_level)

	sextr.bad_pix_mask("segm.fits","bad_pix_stars.txt",xc_star[iii],yc_star[iii])
	f = open(r"modelIN.txt", "w") 
	sys.stdout = f
	header(image_file,'psf_model.fits',"bad_pix_stars.txt",xc_star[iii],yc_star[iii],m0 - 2.5*log10(EXPTIME),pix2secX,pix2secY)
	if window=='moffat':	moffat (xc_star[iii],yc_star[iii],m_star[iii],3.3,ell_star[iii],PA_star[iii])
	if window=='gauss':	gauss (xc_star[iii],yc_star[iii],m_star[iii],3.3,ell_star[iii],PA_star[iii])
	#sky(backgr_level[iii],0.0,0.0)
	sys.stdout = tmp_out
	f.close()
	os.chmod(r"modelIN.txt",0777)
	subprocess.call("galfit modelIN.txt", shell=True)
	hdulist = pyfits.open('psf_model.fits')
	prihdr = hdulist[2].header
	mtot_star = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_MAG']))[0]
	fwhm_star = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_FWHM']))[0]
	if window=='moffat':	beta_star = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_C']))[0]
	theta_star = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_PA']))[0]
	Ell_star = 1.-map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_AR']))[0]
	chi2_star = prihdr['CHI2NU']
	os.remove(r"modelIN.txt") 
	os.remove(r"fit.log") 
	os.remove(r"galfit.01") 
	os.remove(r"bad_pix_stars.txt")
 	os.remove(r"psf_model.fits")
	return fwhm_star*sqrt(pix2secX*pix2secY),Nstars
def PSF_extract(xc,yc,filter_name):
	if filter_name=='u':
		f_n = 0
	if filter_name=='g':
		f_n = 1
	if filter_name=='r':
		f_n = 2
	if filter_name=='i':
		f_n = 3
	if filter_name=='z':
		f_n = 4
	shutil.move('psf.fits','psf1.fit')

	coords = str(xc)+".0 "+str(yc)+".0"
	psf = "read_PSF psf1.fit "+str(f_n)+ " " +str(xc) + " " + str(yc) + " psf.fits"
	subprocess.call(psf, shell=True)
	#exit()
	sky_est.imarith_func('psf.fits','-',1000.,'psf.fits')
	os.remove(r"psf1.fit")
	#exit()
	
