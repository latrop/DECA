#!/usr/bin/python
# Results for the model with spirals
# -*- coding:  cp1251 -*-

import random as random_number
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
from os.path import exists
import fileinput
import pyfits
import re
#import cosmolopy.distance as cd
from scipy import special
import fpdf
#import aplpy
from fpdf import FPDF
tmp_out = sys.stdout


DECA_PATH = os.path.split(os.path.dirname(__file__))[0]

sys.path.append(DECA_PATH)
sys.path.append(DECA_PATH + '/ini_modules')
sys.path.append(DECA_PATH + '/analysis_modules')
sys.path.append(DECA_PATH + '/output_modules')
sys.path.append(DECA_PATH + '/prep_modules')


import setup
import cosmo
import distance as cd

omega_M = setup.omega_M
omega_L = setup.omega_L
omega_k = setup.omega_k
h = setup.h

s = "disc+bulge"





def read_fits(gal): 
# The function to read the fits file.
	hdulist = pyfits.open(gal)
	prihdr = hdulist[2].header
	scidata = hdulist[2].data
	return scidata

def image(gal,name,number,ratio):
	dx = 5.
	dy = 5./ratio
	try:
		fig = aplpy.FITSFigure(gal, hdu=number, figsize=(dx, dy))
		fig.show_grayscale(stretch='log',vmid=-10,vmin=-5, invert=True)
		fig.show_contour(gal, colors='red')
		#fig.tick_labels.set_font(size='small')
		fig.axis_labels.hide() # to hide axis labels
		fig.tick_labels.hide() # to hide tick labels
		#fig.add_label(0.2, 0.9, name, relative=True,size=28)
		if name=='galaxy':	fig.add_scalebar(0.002777)
		if name=='galaxy':	fig.scalebar.set_label(' 10" ')
		pic = name + '.png'
		savefig(pic,bbox_inches=0)
	except:
		print "AplPy module is not found!"

def filters(filter_name):
# from http://www.ucolick.org/~cnaw/sun.html
	if filter_name=='U':
		Msun = 5.61
	if filter_name=='B':
		Msun = 5.48
	if filter_name=='V':
		Msun = 4.83
	if filter_name=='R':
		Msun = 4.42
	if filter_name=='I':
		Msun = 4.08
	if filter_name=='J':
		Msun = 3.64
	if filter_name=='H':
		Msun = 3.32
	if filter_name=='K':
		Msun = 3.28
	if filter_name=='L':
		Msun = 3.25

	#SDSS-bands:
	if filter_name=='u':
		Msun = 6.80
	if filter_name=='g':
		Msun = 5.45
	if filter_name=='r':
		Msun = 4.76
	if filter_name=='i':
		Msun = 4.58
	if filter_name=='z':
		Msun = 4.51
	return Msun

def incl(q_vis,q_0=0.2):
	return np.arccos( np.sqrt(((q_vis)**2 - (q_0)**2) / (1.-(q_0)*2)) )

def Aint_tully(q_vis,a1=0.92,a2=1.63):
	return (a1 + b1*(log10(W)-2.5)) * log10(1./q)

def Aint(q):
	return 2.5*log10(1./q)

def D_L(z):
	cosmo = {'omega_M_0':omega_M, 'omega_lambda_0':omega_L, 'omega_k_0':omega_k, 'h':h}
	d_L = cd.luminosity_distance(z, **cosmo)	#in Mpc
	return d_L

def R(c):
	return math.pi*(c+2.) / (4.*special.beta(1./(c+2.),1.+1./(c+2.)))

def scale(par,D):
	return par*D*1000./206265. 

def magDisk_face_on(m0d,h,q,c):
	return -2.5*log10(2.*math.pi) + m0d - 5.*log10(h) + Aint(q) + 2.5*log10(R(c))

def magDisk_edge_on(S0d,h,z0,c):
	return -2.5*log10(2.*math.pi) + S0d - 2.5*log10(z0*h) + 2.5*log10(R(c))

def MagDisk_edge_on(S0d,h,z0,c,z,filter_name,colour_name,colour_value,Aext):
	return magDisk_edge_on(S0d,h,z0,c) - 5.*log10(D_L(z)) - Aext - cosmo.A_kcor(filter_name, z, colour_name, colour_value) - 25.

def magGaussian(m0,fwhm,q,c):
	return -2.5*log10(2.*math.pi) + m0 - 5.*log10(fwhm) + Aint(q) + 2.5*log10(R(c))

def MagGaussian(m0,fwhm,q,c,z,filter_name,colour_name,colour_value,Aext):
	return magGaussian(m0,fwhm,q,c) - 5.*log10(D_L(z)) - Aext - cosmo.cosmo.A_kcor(filter_name, z, colour_name, colour_value) - 25.

def MagDisk_face_on(m0d,h,q,c,z,filter_name,colour_name,colour_value,Aext):
	return magDisk_face_on(m0d,h,q,c) - 5.*log10(D_L(z)) - Aext - cosmo.A_kcor(filter_name, z, colour_name, colour_value) - 25.

def Luminosity(Magnitude,filter_name):
	return 10**((filters(filter_name)-Magnitude)/2.5)

def magBulge_f(reb,meb,n,q,c):
	nu = 1.9987*n-0.3267
	An = 2.5*log10(2.*math.pi*n/(nu**(2.*n)) * special.gamma(2.0*n)) + 1.0857*nu
	return - 5.*log10(reb) + meb - An - 2.5*log10(q) + 2.5*log10(R(c))

def MagBulge_f(reb,meb,n,q,c,z,filter_name,colour_name,colour_value,Aext):
	return magBulge_f(reb,meb,n,q,c) - 5.*log10(D_L(z)) - Aext - cosmo.A_kcor(filter_name, z, colour_name, colour_value) - 25.

def m0d_corr(m0d,q,Aext,z,a3=0.5):
	return m0d + a3*log10(1./q) - Aext - 2.5*log10( (1.+z)**3 )

def meb_corr(meb,q,Aext,z,a4=0.5):
	return meb + a4*log10(1./q) - Aext - 2.5*log10( (1.+z)**3 )

def MagDisk_face_on_corr(MagDisk,q,a5):
	return MagDisk -a5*log10(1./q)

def MagBulge_corr(MagBulge,q,a6):
	return MagBulge -a6*log10(1./q)


def BT_opt(m0_d,h_d,q_d,z0,c_d,me_bul,re_bul,n_bul,q_bul,c_bul,m0):
	if z0!=99999.:	magDisk = magDisk_edge_on(m0_d,h_d,z0,c_d)
	else:	magDisk = magDisk_face_on(m0_d,h_d,q_d,c_d)
	lumDisk = 10**(0.4*(m0-magDisk))

	magBul = magBulge_f(re_bul,me_bul,n_bul,q_bul,c_bul)
	lumBul = 10**(0.4*(m0-magBul))
	BT = lumBul / (lumBul+lumDisk)
	return BT



def sextr(file_out_sex,number,RA1,DEC1,galaxy_name,filter_name,FWHM,PA1,a_image,b_image,flux_radius,mag_auto,ellip,kron_r,pix2sec):
	q0 = 0.2
	q = 1.-ellip
	a_image = kron_r*a_image*pix2sec
	b_image = kron_r*b_image*pix2sec
	flux_radius = flux_radius*pix2sec
	incl = degrees(arccos(sqrt( (q*q-q0*q0)/(1.-q0*q0) )))
	if incl+3.<=90. and incl+3.>=0.:	incl = incl+3.
	if os.path.exists(file_out_sex)==True:
		f = open(file_out_sex, "a") 
		sys.stdout = f
		print "%i\t%.3f\t%.3f\t%s\t%s\t%.3f\t%.1f\t%.2f\t%.2f\t%.3f\t%.2f\t%.2f\t%.1f\t" % (number,RA1,DEC1,galaxy_name,filter_name,FWHM,PA1,a_image,b_image,q,flux_radius,mag_auto,incl)
		sys.stdout = tmp_out
		f.close()
	else:
		f = open(file_out_sex, "w") 
		sys.stdout = f
		print "#\tRA\tDEC\tname\tF\tfwhm\tPA\tsma\tsmb\tq\tre\tmag\ti"
		print "%i\t%.3f\t%.3f\t%s\t%s\t%.3f\t%.1f\t%.2f\t%.2f\t%.3f\t%.2f\t%.2f\t%.1f\t" % (number,RA1,DEC1,galaxy_name,filter_name,FWHM,PA1,a_image,b_image,q,flux_radius,mag_auto,incl)
		sys.stdout = tmp_out
		f.close()

def results(galaxy_name,model,out_image,file_out_txt,file_out_pdf,filter_name,colour_name,colour_value,Aext,z,pix2sec,m0,D,status,tables='tog',gal_number=0,out=1,ima_number=0,fix=1):
	if out_image!='0':
		file_galfit_outimage = out_image
	else:
		file_galfit_outimage = 'galfit_out.fits'

	if out==0:
		hdulist1 = pyfits.open('subcomps.fits')
		if model=='exp+ser' or model=='edge+ser' or model=='edge' or model=='exp':
			scidata_disk = hdulist1[2].data
		if model=='exp+ser' or model=='edge+ser':
			scidata_bulge = hdulist1[3].data
		if model=='ser':
			scidata_bulge = hdulist1[2].data			
		
	hdulist = pyfits.open(file_galfit_outimage)
	prihdr = hdulist[2].header
	m0_d=[99999.,99999.]
	me_bul=[99999.,99999.]
	z0=[99999.,99999.]
	z_0_k = [99999.,99999.]	
	q_d=[99999.,99999.]
	q_bul=[99999.,99999.]
	lumDisk = 0.
	lumBul = 0.
	MagBul = 99999.
	M_gal = 99999.
	re_bul_k = [99999.,99999.]
	rtr = [99999.,99999.]
	c_bul = [99999.,99999.]

	magTot1 = 99999.
	MagTot1 = 99999.
	MagDisk1 = 99999.
	MagBul1 = 99999.
	LumTot1 = 99999.
	LumDisk1 = 0.
	LumBul1 = 0.
	lumDisk1 = 0.
	lumBul1 = 0.
	BT1 = 99999.
	BD1 = 99999.

	r_in = [99999.,99999.]
	r_out = [99999.,99999.]
	r_angle = [99999.,99999.]
	incl_angle = [99999.,99999.]
	spa = [99999.,99999.]
	pa = [99999.,99999.]


	'''
	xc_d=[9999.,9999.]
	yc_d=[9999.,9999.]
	xc_bul=[9999.,9999.]
	yc_bul=[9999.,9999.]
	m0_d=[9999.,9999.]
	h_d=[9999.,9999.]
	ell_d=[9999.,9999.]
	re_bul=[9999.,9999.]
	ell_bul=[9999.,9999.]
	n_bul=[9999.,9999.]
	'''

	if type(D)==str and type(z)!=str:
		D = D_L(z)
		D_A = D / (1.+z)**2
		Scale='kpc'
	elif type(D)!=str and type(z)!=str:
		D_A = D / (1.+z)**2
		Scale='kpc'
	elif type(D)!=str and type(z)==str:
		Scale='kpc'
		z = D*73./300000.
		D_A = D 
	else:
		z = 'NO'
		Scale='arcsec'

	if model=="exp" or model=="exp+ser":
			m0_d = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_MAG']))
			h_d = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_RS']))
			PA_d =  map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_PA']))
			q_d = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_AR']))
			xc_d = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_XC']))
			yc_d = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_YC']))
			r_in = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_RIN']))
			r_out = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_ROUT']))
			r_angle = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_RANG']))
			incl_angle = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_INCL']))
			spa = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_SPA']))
			pa = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_PA']))

			if '1_C0' in prihdr:
				c_d = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_C0']))
			else:
				c_d = [0.,0.]
			h_d[0] = h_d[0]*pix2sec
			h_d[1] = h_d[1]*pix2sec
			magDisk = magDisk_face_on(m0_d[0],h_d[0],q_d[0],c_d[0])
			
			lumDisk = 10**(0.4*(m0-magDisk))

			if out==0:			
				# REAL LUMINOSITY:
				lumDisk1 = sum(scidata_disk)
				magDisk1 = m0-2.5*log10(lumDisk1)
				if z!='NO':	MagDisk1 = magDisk1 - 5.*log10(D_L(z)) - Aext - cosmo.A_kcor(filter_name, z, colour_name, colour_value) - 25.
				else:	MagDisk1 = 99999.


			if z!='NO':
				MagDisk = MagDisk_face_on(m0_d[0],h_d[0],q_d[0],c_d[0],z,filter_name,colour_name,colour_value,Aext)
				LumDisk = Luminosity(MagDisk,filter_name)
				h_d_k = scale(h_d[0],D)
				h_d_k_e = scale(h_d[1],D)
			else:
				MagDisk = 99999.
				LumDisk = 99999.
				h_d_k = 99999.
				h_d_k_e = 99999.




	if model=="exp+ser" and fix==1:
		def magBulge_f(reb,meb,n,q,c):
			nu = 1.9987*n-0.3267
			An = 2.5*log10(2.*math.pi*n/(nu**(2.*n)) * special.gamma(2.0*n)) + 1.0857*nu
			return - 5.*log10(reb) + meb - An - 2.5*log10(q) + 2.5*log10(R(c))
		def MagBulge_f(reb,meb,n,q,c,z,filter_name,colour_name,colour_value,Aext):
			return magBulge_f(reb,meb,n,q,c) - 5.*log10(D_L(z)) - Aext - cosmo.A_kcor(filter_name, z, colour_name, colour_value) - 25.
		me_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_MU_E']))
		re_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_RE']))
		n_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_N']))
		PA_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_PA']))
		if setup.fix_ell==1:
			q_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_AR']))
		else:
			q_bul = [1.,0.]
		xc_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_XC']))
		yc_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_YC']))
		if '2_C0' in prihdr:
			c_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_C0']))
		else:
			c_bul = [0.,0.]
		re_bul[0] = pix2sec*re_bul[0]
		re_bul[1] = pix2sec*re_bul[1]

		magBul = magBulge_f(re_bul[0],me_bul[0],n_bul[0],q_bul[0],c_bul[0])
		lumBul = 10**(0.4*(m0-magBul))

		if out==0:			
			# REAL LUMINOSITY:
			lumBul1 = sum(scidata_bulge)
			magBul1 = m0-2.5*log10(lumBul1)
			if z!='NO':	MagBul1 = magBul1 - 5.*log10(D_L(z)) - Aext - cosmo.A_kcor(filter_name, z, colour_name, colour_value) - 25.
			else:	MagBul1 = 99999.

		if z!='NO':
			MagBul = MagBulge_f(re_bul[0],me_bul[0],n_bul[0],q_bul[0],c_bul[0],z,filter_name,colour_name,colour_value,Aext)
			LumBul = Luminosity(MagBul,filter_name)
			re_bul_k = [scale(re_bul[0],D),scale(re_bul[1],D)]
		else:
			MagBul = 99999.
			LumBul = 99999.
			re_bul_k = [99999.,99999.]
		#z0_k=99999.
		#z0_k_e=99999.



	chi2 = prihdr['CHI2NU']

	#BT = lumBul/(lumDisk+lumBul)
	#print m0
	m_gal = m0 - 2.5*log10(lumBul+lumDisk)
	#print magDisk
	#print magBul
	if type(z)!=str:	M_gal = m_gal - 5.*log10(D_L(z)*10**6) + 5.
	BT = lumBul / (lumBul+lumDisk)



	if out==0:
		if MagBul1!=99999.:	LumBul1 = Luminosity(MagBul1,filter_name)
		else:	LumBul1 = 0.	
		if MagDisk1!=99999.:	LumDisk1 = Luminosity(MagDisk1,filter_name)
		else:	LumDisk1 = 0.
		LumGal1 = LumBul1 + LumDisk1

		magGal1 = m0-2.5*log10(lumBul1+lumDisk1)

		BT1 = LumBul1 / (LumBul1+LumDisk1)
		BD1 = LumBul1 / (LumDisk1)

		MagGal1 = magGal1 - 5.*log10(D_L(z)*10**6) + 5.


	if m0_d[0]==99999.:
		xc_d=[99999.,99999.]
		yc_d=[99999.,99999.]
		h_d=[99999.,99999.]
		h_d_k = 99999.
		h_d_k_e = 99999.
		z0=[99999.,99999.]
		z_0_k=[99999.,99999.]
		PA_d=[99999.,99999.]
		q_d=[99999.,99999.]
		MagDisk=99999.
		magDisk=99999.
		BT=1.
		BT1=1.
		BD1=99999.
	if me_bul[0]==99999.:
		xc_bul=[99999.,99999.]
		yc_bul=[99999.,99999.]
		re_bul=[99999.,99999.]
		re_bul_k=[99999.,99999.]
		n_bul=[99999.,99999.]
		PA_bul=[99999.,99999.]
		q_bul=[99999.,99999.]
		MagBulge=99999.
		magBulge=99999.
		BT=0.
		BT1=0.
		BD1=0.



	if out==1 and fix==0:                        # here the disk scale was corrected (Sergey). q_d was added (Sergey)
		return xc_d[0],yc_d[0],m0_d[0],h_d[0]/pix2sec,r_in[0],r_out[0],r_angle[0],incl_angle[0],spa[0],pa[0], q_d[0]
	if out==1 and fix==1:
		return xc_d[0],yc_d[0],m0_d[0],h_d[0],r_in[0],r_out[0],r_angle[0],incl_angle[0],spa[0],pa[0],rtr[0],xc_bul[0],yc_bul[0],me_bul[0],re_bul[0],q_bul[0],n_bul[0],c_bul[0],chi2








	if tables=="tog" and out==0:
		f = open(file_out_txt, "a") 
		sys.stdout = f
		print "%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%i\t%i" % (gal_number,m0_d[0],m0_d[1],h_d[0],h_d[1],h_d_k,h_d_k_e,PA_d[0],PA_d[1],z0[0],z0[1],z_0_k[0],z_0_k[1],rtr[0],rtr[1],MagDisk,MagDisk1,LumDisk1,me_bul[0],me_bul[1],re_bul[0],re_bul[1],re_bul_k[0],re_bul_k[1],n_bul[0],n_bul[1],PA_bul[0],PA_bul[1],q_bul[0],q_bul[1],c_bul[0],c_bul[1],MagBul,MagBul1,LumBul1,BT,BT1,BD1,MagGal1,LumGal1,magGal1,chi2,model,status,ima_number)
		sys.stdout = tmp_out
		f.close()


	if (file_out_pdf!='NONE' or file_out_pdf!='NO') and status==0:
		galfit_out = 'galfit_out.fits'
		scidata=read_fits(galfit_out)
		[ny,nx] = scidata.shape
		l=200.
		ratio=float(nx)/float(ny)
		#ratio = 3
		dimy = l/(3.*ratio)
		dimx=dimy*ratio


		try:
			import ds9
			ds9.ds9_image(galfit_out)
			jpeg = 0
		except:
			import aplpy
			image(galfit_out,'galaxy',1,ratio)
			image(galfit_out,'model',2,ratio)
			image(galfit_out,'residual',3,ratio)
			jpeg = 1

		pdf=FPDF()
		pdf.add_page()
		pdf.set_font('Times','B',18)
		pdf.set_text_color(220,50,50)
		pdf.cell(170,10,galaxy_name,0,1,'C')
		pdf.set_text_color(0,0,0)
		pdf.set_font('Times','B',14)
  
  
		if model=='exp+ser':
				text = 'Model: %s\n Disk: m0d=%.2f +/- %.2f mag arcsec-2, h=%.2f +/- %.2f arcsec, h_k=%.2f +/- %.2f kpc, e=%.2f +/- %.2f, PA=%.2f +/- %.2f deg\n Bulge: meb=%.2f +/- %.2f mag arcsec-2, reb=%.2f +/- %.2f arcsec, reb_k=%.2f +/- %.2f kpc, n=%.2f +/- %.2f, e=%.2f +/- %.2f, PA=%.2f +/- %.2f deg, c0=%.1f\n Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, chi2=%.2f' % (model,m0_d[0],m0_d[1],h_d[0],h_d[1],h_d_k,h_d_k_e,1.-q_d[0],q_d[1],PA_d[0],PA_d[1],me_bul[0],me_bul[1],re_bul[0],re_bul[1],re_bul_k[0],re_bul_k[1],n_bul[0],n_bul[1],1.-q_bul[0],q_bul[1],PA_bul[0],PA_bul[1],c_bul[0],magGal1,MagGal1,BT1,chi2)


		if model=='exp':
				text = 'Model: %s\n Disk: m0d=%.2f +/- %.2f mag arcsec-2, h=%.2f +/- %.2f arcsec, h_k=%.2f +/- %.2f kpc, e=%.2f +/- %.2f, PA=%.2f +/- %.2f deg\n  Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, chi2=%.2f' % (model,m0_d[0],m0_d[1],h_d[0],h_d[1],h_d_k,h_d_k_e,1.-q_d[0],q_d[1],PA_d[0],PA_d[1],magGal1,MagGal1,BT1,chi2)

		if model=='edge+ser':
				text = 'Model: %s\n Disk: m0d=%.2f +/- %.2f mag arcsec-2, h=%.2f +/- %.2f arcsec, h_k=%.2f +/- %.2f kpc, z0=%.2f +/- %.2f arcsec, z0_k=%.2f +/- %.2f kpc\n Bulge: meb=%.2f +/- %.2f mag arcsec-2, reb=%.2f +/- %.2f arcsec, reb_k=%.2f +/- %.2f kpc, n=%.2f +/- %.2f, e=%.2f +/- %.2f, c0=%.1f\n Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, Rtr=%.2f arcsec, chi2=%.2f' % (model,m0_d[0],m0_d[1],h_d[0],h_d[1],h_d_k,h_d_k_e,z0[0],z0[1],z_0_k[0],z_0_k[1],me_bul[0],me_bul[1],re_bul[0],re_bul[1],re_bul_k[0],re_bul_k[1],n_bul[0],n_bul[1],1.-q_bul[0],q_bul[1],c_bul[0],magGal1,MagGal1,BT1,rtr[0],chi2)

		if model=='edge':
				text = 'Model: %s\n Disk: m0d=%.2f +/- %.2f mag arcsec-2, h=%.2f +/- %.2f arcsec, h_k=%.2f +/- %.2f kpc, z0=%.2f +/- %.2f arcsec, z0_k=%.2f +/- %.2f kpc\n Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, Rtr=%.2f arcsec, chi2=%.2f' % (model,m0_d[0],m0_d[1],h_d[0],h_d[1],h_d_k,h_d_k_e,z0[0],z0[1],z_0_k[0],z_0_k[1],magGal1,MagGal1,BT1,rtr[0],chi2)

		if model=='ser':
				text = 'Model: %s\n Bulge: meb=%.2f +/- %.2f mag arcsec-2, reb=%.2f +/- %.2f arcsec, reb_k=%.2f +/- %.2f kpc, n=%.2f +/- %.2f, e=%.2f +/- %.2f, PA=%.2f +/- %.2f deg, c0=%.1f\n Galaxy: m=%.2f mag, M=%.2f mag, B/T=%.2f, chi2=%.2f' % (model,me_bul[0],me_bul[1],re_bul[0],re_bul[1],re_bul_k[0],re_bul_k[1],n_bul[0],n_bul[1],1.-q_bul[0],q_bul[1],PA_bul[0],PA_bul[1],c_bul[0],magGal1,MagGal1,BT1,chi2)
		
  		pdf.write(5,text)

		if dimy>dimx:	dimy = dimx
		if model=='edge' or model=='edge+ser':
 	 		pdf.image('galaxy.png',x=5,y=51.,h=dimy,w=dimx)
	  		pdf.image('model.png',x=75,y=51.,h=dimy,w=dimx)
			pdf.image('residual.png',x=145,y=51.,h=dimy,w=dimx)
			try:
				pdf.image('ell.png',x=5,y=100.,w=65.,h=160.)
			except:
				print 'There is no ell.png!'
			pdf.image('edge_prof.png',x=75,y=115,w=130.)


		if model=='exp' or model=='exp+ser' or model=='ser':
 	 		pdf.image('galaxy.png',x=5,y=51.,h=dimy,w=dimx)
	  		pdf.image('model.png',x=75,y=51.,h=dimy,w=dimx)
			pdf.image('residual.png',x=145,y=51.,h=dimy,w=dimx)
			try:
				pdf.image('ell.png',x=5,y=112.,w=80.,h=190.)
			except:
				print 'There is no ell.png!'
			#pdf.image('major_prof.png',x=70,y=110.,w=70.)
			try:
				pdf.image('azim_prof.png',x=100,y=115.,w=90.)
			except:
				print 'There is no azim_prof.png!'
			try:
				pdf.image('radius_prof.png',x=100,y=205,w=90.)
			except:
				print 'There is no radius_prof.png!'
		#pdf.cell(100,100,text,0,1,'B')

		pdf.output(file_out_pdf,'F')





