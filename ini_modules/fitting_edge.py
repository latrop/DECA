#!/usr/bin/python
# -*- coding:  cp1251 -*-

#*** Common modules ***
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
from scipy.odr.odrpack import *
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
import time
import pyfits
import re
from scipy import special
#from pyraf import iraf

DECA_PATH = os.path.split(os.path.dirname(__file__))[0]

sys.path.append(DECA_PATH)
sys.path.append(DECA_PATH + '/ini_modules')
sys.path.append(DECA_PATH + '/analysis_modules')
sys.path.append(DECA_PATH + '/output_modules')
sys.path.append(DECA_PATH + '/prep_modules')

#*** Colour fonts ***
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

#*** USED DECA MODULES ***
import sextr
import changer
import psf
import image 
#import surveys
import ell
import GALFIT
import INITIALS
import out
import res
import sky_est
import INITIALS_edge_on
import setup
#import superellipse
#import edge_soph
#import iter2
#import bulge_pic
import constraints
#import fitting

tmp_out = sys.stdout
import setup

file_out_ini = setup.file_out_ini
file_inimodels_txt = setup.file_inimodels_txt
way_ini = setup.way_ini
one_dec = setup.one_dec

file_weight = 'weight.fits'
file_constraints = 'constr.txt'
file_galfit_input = 'input.txt'
file_field = 'field.fits'
file_psf_field = 'psField.fit'	# for SDSS
file_field1 = 'field1.fits'
file_field2 = 'field2.fits'
file_gal = 'galaxy.fits'
file_sex_field = 'field.cat'
file_sex_galaxy = 'galaxy.cat'
file_gal_clean = 'galaxy_clean.fits'
file_gal_clean_new = 'galaxy_clean1.fits'
file_segm = 'segm.fits'
file_badpix = 'badpix.txt'
file_inter = 'interact.cat'
file_psf = 'psf.fits'
file_psf_model = 'psf_model.fits'
file_nonpar = 'nonpar.txt'
file_galfit_outimage = 'galfit_out.fits'
file_out_txt = setup.file_out_txt

C1 = setup.C1
C2 = setup.C2
Ccr = setup.Ccr
El_cr = setup.El_cr
ellip_cr = setup.ellip_cr



def add_mask_dust(zmin1,zmin2,xc,rmax,nx):
	f_dust = open(file_badpix, "a") 
	for i in range(zmin1,zmin2,1):
		for k in range(int(floor(xc-rmax/2.)),int(ceil(xc+rmax/2.)),1):
			print >>f_dust, '%i\t%i' % (k+1,i+1)
			#print k+1,i+1
	#print int(floor(xc-rmax/2.)),int(ceil(xc+rmax/2.))
	
	f_dust.close()	
	#exit()

def line(B, x):
    return (1.0857/B[0])*fabs(x) + B[1]

def disk_edge(B,z):
	#*** For edge-on disk SB in mag/arcsec^2 (along z-axis).  Sech^2-law. ***
	# B[0] = I0d
	# B[1] = z0 
	z = np.array(z)	
	return B[0] * ( 1. / (np.cosh(fabs(z)/B[1]))**2   )

def bump_find(r,I,rmax,m0,pix2sec,backgr_dev,fwhm):


	# Function to find a bulge 
	rmin_bulge_coeff = setup.rmin_bulge
	rmin_bulge = rmin_bulge_coeff*fwhm/pix2sec

	r_int = []
	Int = []

	for k in range(len(r)):
		if I[k]>backgr_dev and r[k]>rmin_bulge:
			r_int.append(r[k])
			Int.append(I[k])

	rr = []
	mag = []


	for k in range(len(r_int)):
		if r_int[k]>rmax/2. - (rmax/2.)/4. and  r_int[k]<rmax- (rmax/2.)/4.:
			rr.append(r_int[k]*pix2sec)
			mag.append(-2.5*log10(Int[k]) + m0 + 5.*log10(pix2sec))


	dataToFit = RealData(rr, mag)
	mo = Model(line)
	fitting = ODR(dataToFit, mo, [1.,1.])
	fitting.set_job()
	fit_result = fitting.run()
	h = fit_result.beta[0]/pix2sec  # now in pixels
	m0d = fit_result.beta[1]	# in mag/arcsec^2

	r_int = np.array(r_int)
	Int = np.array(Int)

	#plt.figure(102)
	#plt.plot(np.array(r_int)*pix2sec,-2.5*log10(Int) + m0+ 5.*log10(pix2sec),'o',color='black',markersize=3)
	#plt.plot(rr,mag,'o',color='blue',markersize=6)
	#plt.plot(np.array(r_int)*pix2sec,line([h*pix2sec,m0d],np.array(r_int)*pix2sec),'r',color='blue',lw=2)
	#plt.plot(r_int*pix2sec,m0+5.*log10(pix2sec)-2.5*log10(Int),'^')
	#ylim(30,18)




	def Iexp(r,h,I0d):
		# Exponential disc
		Iexp = I0d * np.exp( -fabs(r)/h )
		return Iexp


	I_bulge = Int # - Iexp(r_int,h,10**(0.4*(m0 - m0d + 5.*log10(pix2sec))))
	max_bul=max(I_bulge)
	I_bulge = I_bulge/max(I_bulge)
	r1 = []
	I_bulge1 = []

	for k in range(len(I_bulge)):
		if I_bulge[k]>0.2:
			I_bulge1.append(Int[k])
			r1.append(r_int[k])
	I_bulge1 = np.array(I_bulge1)
	r1 = np.array(r1)



	try:
		dataToFit = RealData(r1, I_bulge1)
		mo = Model(disk_edge)
		fitting = ODR(dataToFit, mo, [max(I_bulge1),max(r1)/2.])
		fitting.set_job()
		fit_result = fitting.run()
		A = fit_result.beta[0]
		B = fit_result.beta[1]

	except:
		'Fitting to find bump parameters failed!'



	#plt.plot(np.array(r1)*pix2sec,-2.5*log10(I_bulge1) + m0+ 5.*log10(pix2sec),'o',color='black',markersize=5)
	#plt.plot(np.array(r1)*pix2sec,-2.5*log10(disk_edge([A,B*pix2sec], r1*pix2sec)) + m0+ 5.*log10(pix2sec),'r',color='blue',lw=4)
	#plt.show()
	#exit()
	def func1(rr):
		return -2.5*log10(disk_edge([A,B*pix2sec], rr)) + m0+ 5.*log10(pix2sec) - line([h*pix2sec,m0d],rr)

	try:
		rmin = (sp.optimize.bisect(func1,0,rmax*pix2sec))/pix2sec	# in pix
		if rmin>rmin_bulge and rmin<rmax:
			print 'There is a bump! rmin=%.3f' % (rmin)
			bump = 0	# There is a bump
		else:
			print 'There is no significant bump! rmin=%.3f' % (rmin)
			bump = 1	# There is no significant bump

	except :
			print 'Could not find the bump parameters!'
			bump = 1	# Maybe there is a bump	


	return bump



def cut_major(file_ima,xc_ima,yc_ima,theta):
		#1. Profile cut along the major axis
		if theta<0.:	theta = 360+theta
		f = open("cuts.cl", "w") 
		sys.stdout = f
		print "#Cuts"
		print "set imtype=fits"
		print "pvector %s xc=%i yc=%i theta=%.1f vec_output= \"profileima.fits\" out_type=\"image\"" % (file_ima,xc_ima,yc_ima,theta)
		print "logout"
		sys.stdout = tmp_out
		f.close()

def major_prof(m0,rmax,backgr_dev,pix2sec):
		# Function to find a bulge domination radius rmin
		hdulist = pyfits.open("profileima.fits")
		scidata = hdulist[0].data
		prihdr = hdulist[0].header
		nx = int(prihdr['NAXIS1'])
		xc = list(scidata).index(max(scidata))

		r_int = []
		Int = []
		
		for k in range(0, nx):
			if scidata[k]>2.*backgr_dev:			# !!!!!!!!!!!!!!!
				Int.append(scidata[k]) # In case if there is sky subtraction in these pixels
				r_int.append(k - xc + 0.5)		# add 0.5pix shift (BUDDA)


		#plt.plot(r, Int,'r-',color='red', lw=3)
		#plt.show()

		rr = []
		mag = []

		for k in range(len(r_int)):
			if r_int[k]>rmax/2. - (rmax/2.)/4. and  r_int[k]<rmax - (rmax/2.)/4.:
				rr.append(r_int[k]*pix2sec)
				mag.append(-2.5*log10(Int[k]) + m0 + 5.*log10(pix2sec))
		#print len(r_int)
		dataToFit = RealData(rr, mag)
		mo = Model(line)
		fitting = ODR(dataToFit, mo, [1.,1.])
		fitting.set_job()
		fit_result = fitting.run()
		h = fit_result.beta[0]/pix2sec # now in pixels
		m0d = fit_result.beta[1]	# in mag/arcsec^2



		# Averaging:
		rr = []
		II = []
		for k in range(len(r_int)):
			if r_int[k]<0.:
				rr.append(int(floor(r_int[k])))
				II.append(Int[k])
			if r_int[k]>0.:
				rr.append(int(ceil(r_int[k])))
				II.append(Int[k])
		rrr = []
		III = []

		k_left = rr.index(-1)
		k_right = rr.index(1)

		for x in arange(k_right,len(rr),1):		
			for y in arange(0,k_left,1):
				if fabs(rr[x])==fabs(rr[y]):
					rrr.append(rr[x])
					III.append((II[x]+II[y])/2.)
		#plt.plot(rr,II)
		#plt.plot(rrr,III)

		#plt.show()
		#exit()


		return rrr, III, h, m0d 







def magBulge_meb (reb,magBulge,n):
	nu = 1.9987*n-0.3267
	fn = n*special.gamma(2.0*n)/(nu**(2.0*n))
	m0b = magBulge + 2.5*(log10(fn) + log10(2.0*math.pi*reb*reb))
	meb = m0b+1.0857*nu
	return meb

def magDisk_m0d(magDisk,h):
	return magDisk + 2.5*log10(2.*math.pi) + 5.*log10(h)



def way_corrs(model,flux_radius,mag_auto,C31,ellip,a_image,FWHM,pix2sec):
	reb0 = 99999.
	meb0 = 99999.
	m0d0 = 99999.
	magBul0 = 99999.
	magDisk0 = 99999.
	n0 = 99999.
	ellb0 = 99999.
	cb0 = 99999.
	h0 = 99999.
	z00 = 99999.
	elld0 = 99999.
	cd0 = 99999.
	rtr0 = 99999.

	if model=='ser': 
		reb0 = flux_radius
		magBul0 = mag_auto
		if C31>Ccr:	n0 = 3.8
		else:	n0 = 2.
		ellb0 = ellip
		cb0 = 0.
		meb0 = magBulge_meb (reb0*pix2sec,magBul0,n0)

	if model=='exp':
		h0 = flux_radius/1.678
		magDisk0 = mag_auto
		#m0d0 = 99999.
		z00 = 99999.
		elld0 = ellip
		cd0 = 0.
		rtr0 = a_image
		m0d0 = magDisk_m0d(magDisk0,h0*pix2sec)

	if model=='edge':
		h0 = flux_radius/1.678
		magDisk0 = mag_auto
		#m0d0 = 99999.
		z00 = h0/4.
		elld0 = 99999.
		cd0 = -1.
		rtr0 = a_image
		m0d0 = magDisk_m0d(magDisk0,h0*pix2sec) + 2.5*log10(z00/h0)	

	if model=='exp+ser':
		#reb0 = 2.*FWHM
		reb0 = 0.056*a_image*(C31-2.)
		def BT_c31(C31):
			return 0.18*(C31-2.5)

		if C31>Ccr:
			n0 = 3.0
			magBul0 = mag_auto - 2.5*log10(BT_c31(C31))
			magDisk0 = mag_auto -2.5*log10(1.-BT_c31(C31))
		else:
			n0 = 1.5
			magBul0 = mag_auto - 2.5*log10(BT_c31(C31))
			magDisk0 = mag_auto -2.5*log10(1.-BT_c31(C31))
			
		ellb0 = 0.
		cb0 = 0.
		h0 = a_image/3.4
		z00 = 99999.
		elld0 = ellip
		cd0 = 0.
		rtr0 = a_image
		meb0 = magBulge_meb (reb0*pix2sec,magBul0,n0)
		m0d0 = magDisk_m0d(magDisk0,h0*pix2sec)

	if model=='edge+ser':
		def BT_ell(x):
			A = 0.493
			B = 0.871	
			C = 0.830
			if x>=B:
				return 0.05
			else:
				return A*log10(-x+B) + C


		if ellip<El_cr:
			n0 = 3.0
			if BT_ell(ellip)>0. and BT_ell(ellip)<=1.:
				magBul0 = mag_auto - 2.5*log10(BT_ell(ellip))
				magDisk0 = mag_auto -2.5*log10(1.-BT_ell(ellip))
			else:
				magBul0 = mag_auto - 2.5*log10(0.05)
				magDisk0 = mag_auto -2.5*log10(1.-0.05)				
		else:
			n0 = 1.5
			if BT_ell(ellip)>0. and BT_ell(ellip)<=1.:
				magBul0 = mag_auto - 2.5*log10(BT_ell(ellip))
				magDisk0 = mag_auto -2.5*log10(1.-BT_ell(ellip))
			else:
				magBul0 = mag_auto - 2.5*log10(0.05)
				magDisk0 = mag_auto -2.5*log10(1.-0.05)		
			
		ellb0 = 0.
		cb0 = 0.
		h0 = a_image/5.
		z00 = h0/4.
		reb0 = h0*0.3
		elld0 = 99999.
		cd0 = -1.
		rtr0 = a_image
		meb0 = magBulge_meb (reb0*pix2sec,magBul0,n0)
		m0d0 = magDisk_m0d(magDisk0,h0*pix2sec) + 2.5*log10(z00/h0)

	#return m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0
	return m0d0,99999.,h0,z00,elld0,cd0,rtr0,meb0,99999.,reb0,n0,ellb0,cb0




def def_model(ellip,C31,bump):
	model = []
	if ellip<=ellip_cr and bump==0:
		model.append('exp+ser')
		model.append('ser')
		model.append('exp')

	if ellip<=ellip_cr and bump==1:
		model.append('exp+ser')
		model.append('exp')
		model.append('ser')

	if ellip>ellip_cr and ellip<0.7 and bump==0:
		model.append('exp+ser')
		model.append('edge+ser')
		model.append('exp')		
		model.append('ser')

	if ellip>ellip_cr and ellip<0.7 and bump==1:
		model.append('exp')
		model.append('ser')
		model.append('exp+ser')		
		model.append('edge+ser')

	if ellip>=0.7 and bump==0:
		model.append('edge+ser')
		model.append('edge')
		model.append('exp+ser')
		model.append('exp')

	if ellip>=0.7 and bump==1:
		model.append('edge')
		model.append('edge+ser')
		model.append('exp+ser')
		model.append('exp')		
	#print model			
	return model

def filters(f):
	if f=='J' or f=='H' or f=='K' or f=='Ks' or f=='NIR' or f=='i' or f=='r':
		return 0
	else:
		return 1


def ini_load(i,pix2sec):
		# !!! Surface Brightness should be given in mag/arcsec^2 and sizes in arcsec !!!
		# Then sizes will be transformed in pixels!!! 
		if os.path.exists('initials.py')==False and os.path.exists('initials.dat')==False:
			sys.exit(bcolors.FAIL+ 'The input file with initials is not found!'+ bcolors.ENDC)
		elif os.path.exists('initials.py')==True:
			print bcolors.OKBLUE+ 'The program is using your initial parameters from the file initials.py!' + bcolors.ENDC
			import initials

			m0d0 = initials.m0_d
			magDisk0 = initials.magDisk
			if initials.h_d!=99999.:
				h0 = (initials.h_d)/pix2sec
			else:
				h0 = 99999.
			if initials.z0!=99999.:
				z00 = (initials.z0)/pix2sec
			else:
				z00 = 99999.
			elld0 = initials.ell_d
			cd0 = initials.c_d
			if initials.r_tr!=99999.:
				rtr0 = (initials.r_tr)/pix2sec
			else:
				rtr0 = 99999.
		
			meb0 = initials.me_bul
			magBul0 = initials.magBul
			if initials.re_bul!=99999.:
				reb0 = (initials.re_bul)/pix2sec
			else:
				reb0 = 99999.
			n0 = initials.n_bul
			ellb0 = initials.ell_bul
			cb0 = initials.c_bul
			return m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0


		elif os.path.exists('initials.dat')==True:
			print bcolors.OKBLUE+ 'The program is using your initial parameters from the file initials.dat!' + bcolors.ENDC
			#*** PARAMETERS FROM THE INPUT FILE deca_input.dat as a list of objects ***
			with open('initials.dat', 'r') as f:
				lines = f.readlines()
				num_lines = len([l for l in lines if l.strip(' \n') != ''])
			n_it = num_lines

			#	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14
			#	#	filter	S0d	mag_d	h_d	z0_d	ell_d	c_d	r_tr	me_bul	mag_bul	re_bul	n_bul	ell_bul	c_bul	

			m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0 = loadtxt('initials.dat', usecols=[2,3,4,5,6,7,8,9,10,11,12,13,14],dtype=float, unpack=True,skiprows=0)
			numbs = loadtxt('initials.dat', usecols=[0],dtype=int, unpack=True,skiprows=0)
			filts = loadtxt('initials.dat', usecols=[1],dtype=str, unpack=True,skiprows=0)

			if h0[i]!=99999.:
				h0[i] = h0[i]/pix2sec
			
			if z00[i]!=99999.:
				z00[i] = z00[i]/pix2sec

			if rtr0[i]!=99999.:
				rtr0[i] = rtr0[i]/pix2sec

			if reb0[i]!=99999.:
				reb0[i] = reb0[i]/pix2sec

			return m0d0[i],magDisk0[i],h0[i],z00[i],elld0[i],cd0[i],rtr0[i],meb0[i],magBul0[i],reb0[i],n0[i],ellb0[i],cb0[i]



			




def  main(ff_log,number,xc,yc,kron_r,flux_radius,mag_auto,a_image,b_image,PA,NOISE,pix2sec,m0,FWHM,EXPTIME,ini_use,find_model,model,ellip,C31,sky_level,file_out_txt,file_out_pdf,filter_name,colour_name,colour_value,Aext,z,D,galaxy_name,find_psf,rmax,mode):
	if not os.path.exists("./pics/%i" % (number)):
		os.makedirs("./pics/%i" % (number))


	hdulist = pyfits.open(file_gal_clean)
	scidata = hdulist[0].data
	ny,nx = np.shape(scidata)


	ff = open('models.txt', 'w')
	f_ini = open('ini.txt', 'w')

	# Find out if there is a bump (i.e. bulge) in the SB profile:

	try:
		r,I = loadtxt('ellipse.txt', usecols=[1,2], unpack=True, skiprows = 6)
	except:
		cut_major(file_gal_clean,xc,yc,PA)
		os.chmod(r"cuts.cl",0777)
		subprocess.call("cl < cuts.cl -o", shell=True)
		r,I,h,m0d = major_prof(m0,rmax,NOISE,pix2sec)

	bump = bump_find(r,I,rmax,m0,pix2sec,NOISE,FWHM)

	

	#model = def_model(ellip,C31,bump)
	model = ['edge+ser','edge']

	status=1
	Model = model[0]
	k = 0
	way = way_ini	
	good = 1
	BT0 = -1
	rad_ratio0 = 2
	nbul0 = -1
	deltaSB0 = -1
	out_nan=0
	k_iter = 1

	if mode==4:	good=0
	NIR = filters(filter_name) # NIR = 0
	NIR = 0
	change_rmax = 0

	print 'Models to be checked: %s' % (model)
	print >>ff_log, 'Models to be checked: %s' % (model)
	n_first = 0
	while (status!=0 and k<len(model) and n_first==0)  or (k<2):
		print 'Model %s is now considering!' % (Model)
		if ini_use=='NO':
			if way==2:
				try:
					print 'Way to find initial parameters: 2'
					print >>ff_log, 'Way to find initial parameters: 2'
					if Model=='exp' or Model=='exp+ser' or Model=='ser':
						m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0 = INITIALS.ini_prof('azim',xc,yc,PA,m0,pix2sec,NOISE,a_image*kron_r,FWHM,Model,ellip,C31)
						good = 0
						if math.isnan(float(m0d0))==True or n0<0. or math.isnan(float(meb0))==True or reb0<=1.:
							good=1
						
					elif Model=='edge' or Model=='edge+ser':
						if NIR==0:
							m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0,xc,yc,incl,zmin1,zmin2 = INITIALS_edge_on.main_edge(file_gal_clean,xc,yc,kron_r,a_image,b_image,NOISE,pix2sec,m0,FWHM,Model,ellip,C31,NIR)
							good = 0
						else:
							good = 1
						if math.isnan(float(m0d0))==True or n0<0. or math.isnan(float(meb0))==True or reb0<=1.:
							good=1

				except:
					print bcolors.FAIL+ 'This way failed!'+ bcolors.ENDC
					print >>ff_log, '  This way failed!'
					good=1
					NIR = 1

			

			if way==3 or good!=0:
				try:
					print 'Way to find initial parameters: 3'
					print >>ff_log, 'Way to find initial parameters: 3'
					if Model=='exp' or Model=='exp+ser' or Model=='ser':
						m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0 = INITIALS.ini_prof('cut',xc,yc,PA,m0,pix2sec,NOISE,a_image*kron_r,FWHM,Model,ellip,C31)
						good = 0
					elif Model=='edge' or Model=='edge+ser':
						#m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0 = INITIALS.ini_prof('azim',xc,yc,PA,m0,pix2sec,NOISE,a_image*kron_r,FWHM,Model,ellip,C31)
						m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0,xc,yc,incl,zmin1,zmin2 = INITIALS_edge_on.main_edge(file_gal_clean,xc,yc,kron_r,a_image,b_image,NOISE,pix2sec,m0,FWHM,Model,ellip,C31,NIR)
						good = 0
					if math.isnan(float(m0d0))==True or n0<0. or math.isnan(float(meb0))==True or reb0<=1.:
						good=1
				except:
					good=1
					print bcolors.FAIL+ 'This way failed!'+ bcolors.ENDC
					print >>ff_log, '  This way failed!'
			#exit()
			if way==4 or good!=0:
				print 'Way to find initial parameters: 4'
				print >>ff_log, 'Way to find initial parameters: 4'
				m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0 = way_corrs(Model,flux_radius,mag_auto,C31,ellip,a_image*kron_r,FWHM,pix2sec)
				#print m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0,'okkkkkkkk'
				good = 0

		else:
			print 'Given initial parameters will be used!'
			print >>ff_log, 'Given initial parameters will be used!'
			way = 1
			cod = '0'
			m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0 = ini_load(number-1,pix2sec)

		if mode==4:
			if not os.path.exists("%s" % (file_inimodels_txt)):
				ff_ini = open(file_inimodels_txt,'w')
				print >>ff_ini, '#\tm0d0\tmd0\th0\tz0\telld0\tcd0\trtr0\tmeb0\tmbul0\treb0\tn0\tellb0\tcb0\tmodel\tway'
				print >>ff_ini, '\tmas\tmag\tarcsec\tarcsec\t\t\tarcsec\tmas\tmag\tarcsec\t\t\t\t\t\t'
				ff_ini.close()

			ff_ini = open(file_inimodels_txt, 'a')
			print >>ff_ini, '%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%s' % (number,m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0,Model,way)
			ff_ini.close()
			exit()

		if mode==5:
			if not os.path.exists("%s" % (file_inimodels_txt)):
				ff_ini = open(file_inimodels_txt,'w')
				print >>ff_ini, '#\tm0d0\tmd0\th0\tz0\telld0\tcd0\trtr0\tmeb0\tmbul0\treb0\tn0\tellb0\tcb0\tmodel\tway'
				print >>ff_ini, '\tmas\tmag\tarcsec\tarcsec\t\t\tarcsec\tmas\tmag\tarcsec\t\t\t\t\t\t'
				ff_ini.close()

			ff_ini = open(file_inimodels_txt, 'a')
			print >>ff_ini, '%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%s' % (number,m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0,Model,way)
			ff_ini.close()
			return k_iter-1,xc,yc,m0d0,h0,99999.,z00,rtr0,xc,yc,meb0,reb0,99999.,n0,cb0,1.,Model,1
		if Model=='exp' or Model=='ser' or Model=='exp+ser' or change_rmax==1:
			constraints.constr(file_constraints,Model,0,FWHM,pix2sec,h0,rmax)	# Creating constraints file
		else:
			constraints.constr(file_constraints,Model,1,FWHM,pix2sec,h0,rmax)	# Creating constraints file



		if setup.dust_mask==0:
			if filters(filter_name)==1:	# SHOULD BE 1
				#gap = max([1.,z00/4.])
				#zmin1 = min([int(floor(yc - gap)),zmin1])
				#zmin2 = min([int(ceil(yc + gap)),zmin2])			
				add_mask_dust(zmin1,zmin2,xc,rmax,nx)
				#print 'here1',zmin1,zmin2
				#exit()
			else:
				gap = max([1.,z00/6.])
				zmin1 = int(floor(yc - gap))
				zmin2 = int(ceil(yc + gap))
				add_mask_dust(zmin1,zmin2,xc,rmax,nx)
				#print 'here2',zmin1,zmin2
				#exit()
	

		#cb0 = 0.	# 99999.
		#exit()
		if change_rmax == 1:	rmax=99999.; rtr = rmax
	
		status = GALFIT.dec_galfit(find_psf,m0-2.5*log10(EXPTIME),pix2sec,nx,ny,xc,yc,m0d0,magDisk0,h0,z00,elld0,cd0,rmax,meb0,magBul0,reb0,n0,ellb0,cb0,sky_level,PAd=90.,PAb=90.)

		if status==0:
			try:
				#print Model
				xc_d,yc_d,m0_d,h,q_d,z0,rtr,xc_bul,yc_bul,me_bul,re_bul,q_bul,n_bul,c_bul,chi2 = res.results(galaxy_name,Model,'0',file_out_txt,file_out_pdf,filter_name,colour_name,colour_value,Aext,z,pix2sec,m0,D,status,tables='tog',gal_number=0,out=1)
				if Model=='exp+ser':
					BT = res.BT_opt(m0_d,h,q_d,z0,0.,me_bul,re_bul,n_bul,q_bul,c_bul,m0)
					BT0 = BT
					rad_ratio0 = 1.678*h/re_bul
					nbul0 = n_bul
					nu = 1.9987*n_bul-0.3267
					deltaSB0 = m0_d - (me_bul-1.0857*nu)
				elif Model=='edge+ser':
					BT = res.BT_opt(m0_d,h,q_d,z0,-1,me_bul,re_bul,n_bul,q_bul,c_bul,m0)
					BT0 = BT
					rad_ratio0 = 1.678*h/re_bul
					nbul0 = n_bul
					nu = 1.9987*n_bul-0.3267
					deltaSB0 = m0_d - (me_bul-1.0857*nu)
				elif Model=='ser':
					BT = 1.
				else:
					BT = 0.
				#print BT0

			except:
				print >>ff, '%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%i' % (k_iter,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,Model,1)
			try:
				print >>ff, '%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%i' % (k_iter,xc_d,yc_d,m0_d,h,q_d,z0,rtr,xc_bul,yc_bul,me_bul,re_bul,q_bul,n_bul,c_bul,chi2,Model,status)
				shutil.move(file_galfit_outimage,'./pics/%i/out_%i.fits' % (number,k_iter))
				shutil.move('subcomps.fits','./pics/%i/sub_%i.fits' % (number,k_iter))
				os.remove('galfit.01')
			except:
				print 'The model crashed!'
				out_nan = 1
				ini_use=='NO'
				chi2=99999.

		else:
			print >>ff, '%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%i' % (k_iter,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,Model,status)
		print >>f_ini, '%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s' % (k_iter,m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0,way)
		if one_dec=='YES' and status==0:
			if status==0:
				return k_iter,xc_d,yc_d,m0_d,h,q_d,z0,rtr,xc_bul,yc_bul,me_bul,re_bul,q_bul,n_bul,c_bul,chi2,Model,status
			else:
				ffff = open(file_out_txt,'a')
				print >>ffff, "%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%i\t%i" % (number,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,Model,1,0)
				ffff.close()
				exit()

		#print status,change_rmax,k,k_iter,BT0
		#exit()

		

		if status==0 and BT0>0.06 and rad_ratio0>1.3 and BT0<0.75:
			k = 4 # Exit from the cicle

		if status==0 and BT0>=0.75 and nbul0>1.5:
			k = 4

		if status==0 and ((BT0>=0. and BT0<=0.06) or rad_ratio0<1.3 or (BT0>=0.75 and nbul0>0.93 and nbul0<1.07) or k_iter>3):
			Model = model[1]
			k = k + 1
			if k==2: 	k = 3 # Exit from the cicle

		if status==0 and out_nan==1:
			out_nan=0
			Model = model[1]
			k = k + 1

		if status==0 and Model=='edge' and k_iter>2:
			k = 4

		if setup.inner_trunc==1:
			#print rtr,h,change_rmax
			#print 'I am here!%%%'		
			if (status!=0 or rtr<2.*h) and change_rmax==0:
				#print 'I am here!&&&'
				way = 2
				change_rmax = 1
				n_first = 0
				status = 1			
				good = 1
				Model = model[0]
				k = 1
		else:
			if status!=0 and change_rmax==0 and k+1<len(model):
				way = 2
				change_rmax = 1
				n_first = 0
				status = 1			
				good = 1
				Model = model[0]

		if status!=0 and k_iter>3:
			k = 4


		'''
		if status==0 and chi2>2.1 and way<=4:
			way = way + 1
			ini_use = 'NO'
			status = 1			
			good = 1
		'''

		if status!=0 and change_rmax==1 and k_iter!=1:
			Model = model[1]
			k = k + 1
			if k==2: 	k = 5 # Exit from the cicle




		k_iter = k_iter + 1
		#print 'status=%i' % (status)

		'''
		found_mod = 0
		if (status!=0 or chi2>2.1) and k+1<len(model):
			ini_use = 'NO'
			status = 1
			if change_rmax==1: way = way + 1
			good = 1
			#k = k + 1
			change_rmax = 1
			#if k<len(model):	Model = model[k]
			print 'Here 1!'
		print 'status',status,way,k,change_rmax
		

		
		if (BT0>=0.75 and deltaSB0>2.) or (rad_ratio0<1.3 and deltaSB0>0.5):
				ini_use = 'NO'
				status = 1
				way = 2
				good = 1
				n_first = 0
				model = [Model,'ser']
				k = k+1
				Model = model[1]
				found_mod = 1

		if ((BT0>=0. and BT0<=0.01) or (rad_ratio0<1.3 and deltaSB0<0.5) or (BT0>=0.75 and nbul0>0.93 and nbul0<1.07)) and ellip<ellip_cr:
				ini_use = 'NO'
				status = 1
				way = 2
				good = 1
				n_first = 0
				model = [Model,'exp']
				k = k+1
				Model = model[1]
				found_mod = 2

		if ((BT0>=0. and BT0<=0.01) or (rad_ratio0<1.3 and deltaSB0<0.5) or (BT0>=0.75 and nbul0>0.93 and nbul0<1.07)) and ellip>=ellip_cr:
				ini_use = 'NO'
				status = 1
				way = 2
				good = 1
				n_first = 0
				model = [Model,'edge']
				k = k+1
				Model = model[1]
				found_mod = 3
		

		if (status!=0 and way==4 and k+1<=len(model)):
			k = k + 1
			if k<len(model):	Model = model[k]
			way = 2
			good = 1
			print 'Here 2!'


		if status==0 and k+1<=2 and BT0<0.75 and rad_ratio0>=1.3:
			n_first = 1
			k = k + 1
			if k<len(model):	Model = model[k]
			way = 2
			good = 1
			print 'Here 3!'
		k_iter = k_iter + 1
		'''
	ff.close()
	f_ini.close()

	with open('models.txt', 'r') as ff:
		lines = ff.readlines()
		num_lines = len([l for l in lines if l.strip(' \n') != ''])
	n_it = num_lines
	if n_it>1:
		xc_d,yc_d,m0_d,h,q_d,z0,rtr,xc_bul,yc_bul,me_bul,re_bul,q_bul,n_bul,c_bul,chi2 = loadtxt('models.txt', usecols=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], unpack=True, dtype=float)
		Model,status = loadtxt('models.txt', usecols=[16,17], unpack=True, dtype=str)
		numb = loadtxt('models.txt', usecols=[0], unpack=True, dtype=int)
		chi2 = list(chi2)
		kk = chi2.index(min(chi2))
		Model = list(Model)

		if setup.inner_trunc==1:
			if kk==0 and rtr[kk]<2.*h[kk] and chi2[1]<2. and Model[1]=='edge+ser':
				kk = 1

		if k==3 or n_bul[kk]>6. or h[kk]/re_bul[kk]<1.3:
			try:
				kkk = Model.index('edge')
				if status[kkk]=='0':	kk = kkk
				#print 'we are here1'
			except:
				ll = 0

		#print Model[kk],ellip,chi2[Model.index('edge+ser')]
		#exit()

		try:
			if Model[kk]=='edge' and ellip<0.7 and chi2[0]<1.5 and Model[0]=='edge+ser':
				kk = 0
		except:
			kuku = 1
		try:
			if Model[kk]=='edge' and ellip<0.7 and chi2[1]<1.5 and Model[1]=='edge+ser':
				kk = 1
				#print 'we are here111'			
		except:
			kuku = 1

		#print kk,Model[kk]


		m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0 = loadtxt('ini.txt', usecols=[1,2,3,4,5,6,7,8,9,10,11,12,13], unpack=True, dtype=float)
		way0 = loadtxt('ini.txt', usecols=[14], unpack=True, dtype=str)
		k_iter = loadtxt('ini.txt', usecols=[0], unpack=True, dtype=int)

		fd = open(file_out_ini,'a')
		print >>fd, '%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s' % (number,m0d0[kk],magDisk0[kk],h0[kk],z00[kk],elld0[kk],cd0[kk],rtr0[kk],meb0[kk],magBul0[kk],reb0[kk],n0[kk],ellb0[kk],cb0[kk],way0[kk])
		fd.close()
		return numb[kk],xc_d[kk],yc_d[kk],m0_d[kk],h[kk],q_d[kk],z0[kk],rtr[kk],xc_bul[kk],yc_bul[kk],me_bul[kk],re_bul[kk],q_bul[kk],n_bul[kk],c_bul[kk],chi2[kk],Model[kk],status[kk]	

	else:
		way0 = way
		fd = open(file_out_ini,'a')
		print >>fd, '%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s' % (number,m0d0,magDisk0,h0,z00,elld0,cd0,rtr0,meb0,magBul0,reb0,n0,ellb0,cb0,way0)
		fd.close()
		return k_iter-1,xc_d,yc_d,m0_d,h,q_d,z0,rtr,xc_bul,yc_bul,me_bul,re_bul,q_bul,n_bul,c_bul,chi2,Model,status







def second_iter(number,Model,k_iter,nx,ny,xc,yc,pix2sec,m0,EXPTIME,sky_level,file_out_txt,file_out_pdf,filter_name,colour_name,colour_value,Aext,z,D,galaxy_name,find_psf,m0d0,h0,z00,elld0,rmax,meb0,reb0,n0,ellb0,cb0,chi2):

		if rmax>10000.:	rmax = 99999.
		status = GALFIT.dec_galfit(find_psf,m0-2.5*log10(EXPTIME),pix2sec,nx,ny,xc,yc,m0d0,99999.,h0,z00,elld0,-1.,rmax,meb0,99999.,reb0,n0,ellb0,cb0,sky_level,PAd=90.,PAb=90.)


		if status==0:
			xc_d,yc_d,m0_d,h,q_d,z0,rtr,xc_bul,yc_bul,me_bul,re_bul,q_bul,n_bul,c_bul,chi2 = res.results(galaxy_name,Model,'0',file_out_txt,file_out_pdf,filter_name,colour_name,colour_value,Aext,z,pix2sec,m0,D,status,tables='tog',gal_number=0,out=1)


			shutil.move(file_galfit_outimage,'./pics/%i/out_%i.fits' % (number,k_iter+1))
			shutil.move('subcomps.fits','./pics/%i/sub_%i.fits' % (number,k_iter+1))
			os.remove('galfit.01')

		return k_iter+1,xc_d,yc_d,m0_d,h,q_d,z0,rtr,xc_bul,yc_bul,me_bul,re_bul,q_bul,n_bul,c_bul,chi2,Model,status

	
