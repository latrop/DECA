# THE SCRIPT TO RUN SEXTRACTOR, READ ITS OUTPUT FILE AND CREATE MASK FILES
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
import pyfits

DECA_PATH = os.path.split(os.path.dirname(__file__))[0]

sys.path.append(DECA_PATH)
sys.path.append(DECA_PATH + '/ini_modules')
sys.path.append(DECA_PATH + '/analysis_modules')
sys.path.append(DECA_PATH + '/output_modules')
sys.path.append(DECA_PATH + '/prep_modules')

import setup

tmp_out = sys.stdout

#*** FROM SETUP FILE ***
radius_find = setup.radius_find
semimajor_find = setup.semimajor_find

'''
minarea = setup.DETECT_MINAREA
thresh = setup.DETECT_THRESH
filt_name = setup.FILTER_NAME
deb_nthr = setup.DEBLEND_NTHRESH
deb_minc = setup.DEBLEND_MINCONT 
cl_par = setup.CLEAN_PARAM 
ph_auto = str(setup.PHOT_AUTOPARAMS[0]) + "," + str(setup.PHOT_AUTOPARAMS[1])
back_size = setup.BACK_SIZE
filt_size = setup.BACK_FILTERSIZE
back_type = setup.BACKPHOTO_TYPE
mem = setup.MEMORY_PIXSTACK
'''

def run_sextr(survey,file_in,m0,GAIN,pix2sec,fwhm,DECA_PATH,case=0):
		if case==0:
			minarea = setup.DETECT_MINAREA
			thresh = setup.DETECT_THRESH
			filt_name = setup.FILTER_NAME
			deb_nthr = setup.DEBLEND_NTHRESH
			deb_minc = setup.DEBLEND_MINCONT 
			cl_par = setup.CLEAN_PARAM 
			ph_auto = str(setup.PHOT_AUTOPARAMS[0]) + "," + str(setup.PHOT_AUTOPARAMS[1])
			back_size = setup.BACK_SIZE
			filt_size = setup.BACK_FILTERSIZE
			back_type = setup.BACKPHOTO_TYPE
			mem = setup.MEMORY_PIXSTACK
		elif case==1:
			# Slightly diffirent conditions from case=0
			minarea = setup.DETECT_MINAREA
			thresh = setup.DETECT_THRESH
			filt_name = setup.FILTER_NAME
			deb_nthr = setup.DEBLEND_NTHRESH
			deb_minc = setup.DEBLEND_MINCONT 
			cl_par = setup.CLEAN_PARAM 
			ph_auto = str(setup.PHOT_AUTOPARAMS[0]) + "," + str(setup.PHOT_AUTOPARAMS[1])
			back_size = setup.BACK_SIZE
			filt_size = setup.BACK_FILTERSIZE
			back_type = setup.BACKPHOTO_TYPE
			mem = setup.MEMORY_PIXSTACK
		elif survey=='SDSS_dr7':
			# Slightly diffirent conditions from case=0
			minarea = setup.DETECT_MINAREA
			thresh = setup.DETECT_THRESH
			filt_name = setup.FILTER_NAME
			deb_nthr = setup.DEBLEND_NTHRESH
			deb_minc = setup.DEBLEND_MINCONT 
			cl_par = setup.CLEAN_PARAM 
			ph_auto = str(setup.PHOT_AUTOPARAMS[0]) + "," + str(setup.PHOT_AUTOPARAMS[1])
			back_size = setup.BACK_SIZE
			filt_size = setup.BACK_FILTERSIZE
			back_type = setup.BACKPHOTO_TYPE
			mem = setup.MEMORY_PIXSTACK


		# FUNCTION TO RUN SEXTRACTOR WITH DEFAULT INITIAL PARAMETERS
		path1 = DECA_PATH + '/modules/sextractor/default.conv'
		shutil.copy(path1,r"default.conv")
		path2 = DECA_PATH + '/modules/sextractor/default.psf'
		shutil.copy(path2,r"default.psf")
		path3 = DECA_PATH + '/modules/sextractor/default.sex'
		shutil.copy(path3,r"default.sex")
		path4 = DECA_PATH + '/modules/sextractor/default.param'
		shutil.copy(path4,r"default.param")
		path5 = DECA_PATH + '/modules/sextractor/default.nnw'
		shutil.copy(path5,r"default.nnw")
		pars = " -MAG_ZEROPOINT " + str(m0) + " -GAIN " + str(GAIN) + " -PIXEL_SCALE " + str(pix2sec) + " -SEEING_FWHM " + str(fwhm) + " -DETECT_MINAREA " + str(minarea) + " -DETECT_THRESH " + str(thresh) + " -FILTER_NAME " + str(filt_name) + " -DEBLEND_NTHRESH " + str(deb_nthr) + " -DEBLEND_MINCONT " + str(deb_minc) + " -CLEAN_PARAM " + str(cl_par) + " -PHOT_AUTOPARAMS " + str(ph_auto) + " -BACK_SIZE " + str(back_size) + " -BACK_FILTERSIZE " + str(filt_size) + " -BACKPHOTO_TYPE " + str(back_type) + " -MEMORY_PIXSTACK " + str(mem)
		s = "sex " + file_in + " -c" + " default.sex " + pars
		print s
		subprocess.call(s, shell=True)



def bad_pix_mask_stars(file_segm,file_badpix,xc,yc,cl_star,rmax):
		# FUNCTION TO CREATE MASK FILE FOR GALFIT
		hdulist = pyfits.open(file_segm, do_not_scale_image_data=True)
		img = hdulist[0].data
		(dimy,dimx) = img.shape

		f = open(file_badpix, "w") 
		sys.stdout = f

		'''
		for i in range(dimy):
			for k in range(dimx):
				repeat = 0
				for m in range(len(cl_star)):
					if img[i,k]>0 and cl_star[m]>0.9 and img[int(yc[m]),int(xc[m])]==img[i,k] and repeat==0:
						print "%i %i" % (k,i)
						repeat = repeat + 1
		'''

		for m in range(len(cl_star)):
			if cl_star[m]>0.5:
				for i in range(int(floor(yc[m]-rmax[m])),int(ceil(yc[m]+rmax[m])),1):
					for k in range(int(floor(xc[m]-rmax[m])),int(ceil(xc[m]+rmax[m])),1):
						if i<dimy and i>0 and k<dimx and k>0 and img[int(yc[m]),int(xc[m])]==img[i,k]:
							print "%i %i" % (k+1,i+1)




		sys.stdout = tmp_out
		f.close()


def bad_pix_mask(file_segm,file_badpix,xc,yc):
		# FUNCTION TO CREATE MASK FILE FOR GALFIT
		hdulist = pyfits.open(file_segm, do_not_scale_image_data=True)
		img = hdulist[0].data
		#print img

		(dimy,dimx) = img.shape
		#exit()
	
		f = open(file_badpix, "w") 
		sys.stdout = f
		if type(xc)==list:
			for i in range(dimy):
				for k in range(dimx):
					for m in range(len(xc)):
						if img[i,k]>0 and img[i,k]!=img[yc[m],xc[m]]:
							print "%i %i" % (k+1,i+1)
		
		else:
			for i in range(dimy):
				for k in range(dimx):
					if img[i,k]>0 and img[i,k]!=img[yc,xc]:
						print "%i %i" % (k+1,i+1)
		

		sys.stdout = tmp_out
		f.close()	
		
def gal_clean_ima(file_segm,file_gal,file_gal_clean,xc,yc,sky_level):
		# FUNCTION TO CREATE THE IMAGE WITHOUT CONTAMINENTS (THE MASK IS APPLIED)
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
				if img1[i,k]>0 and img1[i,k]!=img1[yc,xc] or img3[i,k]<=-sky_level+1.:
					img3[i,k] = 0.
		
		hdulist3.flush()

def gal_clean_ima_stars(file_segm,file_gal,file_gal_clean,xc,yc,cl_star,rmax):
		# FUNCTION TO CREATE THE IMAGE WITHOUT CONTAMINENTS (THE MASK IS APPLIED)
		hdulist1 = pyfits.open(file_segm, do_not_scale_image_data=True)
		img1 = hdulist1[0].data
		hdulist2 = pyfits.open(file_gal, do_not_scale_image_data=True)
		img2 = hdulist2[0].data

		(dimy,dimx) = img1.shape
		shutil.copy(file_gal,file_gal_clean) 
		hdulist3 = pyfits.open(file_gal_clean, do_not_scale_image_data=True, mode='update')
		img3 = hdulist3[0].data

		for m in range(len(cl_star)):
			if cl_star[m]>0.5:
				for i in range(int(floor(yc[m]-rmax[m])),int(ceil(yc[m]+rmax[m])),1):
					for k in range(int(floor(xc[m]-rmax[m])),int(ceil(xc[m]+rmax[m])),1):
						if i<dimy and i>0 and k<dimx and k>0 and img1[int(yc[m]),int(xc[m])]==img1[i,k]:
							img3[i,k] = 0.

	
		hdulist3.flush()



'''
def gal_clean_ima1(survey,file_segm,file_gal,file_gal_clean,xc,yc):
		# FUNCTION TO CREATE THE IMAGE WITHOUT CONTAMINENTS (THE MASK IS APPLIED)
		if survey=='UKIDSS':
			hdulist1 = pyfits.open(file_segm, do_not_scale_image_data=True)
			img1 = hdulist1[0].data
			hdulist2 = pyfits.open(file_gal, do_not_scale_image_data=True)
			img2 = hdulist2[1].data

			(dimy,dimx) = img1.shape
			from pyraf import iraf
			#shutil.copy(file_gal,file_gal_clean) 
			iraf.imcopy(file_gal+'[1]',file_gal_clean)
			hdulist3 = pyfits.open(file_gal_clean, do_not_scale_image_data=True, mode='update')
			img3 = hdulist3[0].data

			for i in range(dimy):
				for k in range(dimx):
					if img1[i,k]>0 and img1[i,k]!=img1[yc,xc]:
						img3[i,k] = 0.
		
			hdulist3.flush()
		else:
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
'''
def read_sextr_find(RA,DEC,file_out,coords=1):
	# FUNCTION TO FIND THE OBJECT AND RETURN ITS PARAMETERS FROM THE SEXTRACTOR CATALOGUE
	number,x_image,y_image,x_world,y_world,XPEAK_IMAGE,YPEAK_IMAGE,XPEAK_WORLD,YPEAK_WORLD,flux_radius25,flux_radius50,flux_radius75,flux_radius99,mag_auto,a_image,b_image,theta_image,ellipticity,kron_radius,backgr,class_star = loadtxt(file_out, usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], unpack=True, skiprows = 18)



	#if type(x_image)!='numpy.float64':	#dtype for 1 object!
	if type(x_image)==np.ndarray:	#dtype for 1 object!
		#print x_world


		if fabs(max(x_world)-fabs(min(x_world)))>1. and fabs(RA-fabs(max(x_world)))>1.:		#### !!!!!!!!!!!!!!!! EQUAL to 1 DEGREE !!!
			radius_find1 = radius_find
			print 'I am here coords=1'
		elif fabs(max(x_world)-fabs(min(x_world)))<1. and fabs(RA-fabs(max(x_world)))<1.:	
			radius_find1 = radius_find/3600.	
			print 'I am here coords=2'
		elif fabs(max(x_world)-fabs(min(x_world)))<1. and fabs(RA-fabs(max(x_world)))>1.:	
			radius_find1 = radius_find
			x_world = x_image
			y_world = y_image
			print 'I am here coords=3'



		n_iter = len(x_image)

		number_found = 0
		PA1 = [0.,0.]; xc1 = [0.,0.]; yc1 = [0.,0.]; RA_c1 = [0.,0.]; DEC_c1 = [0.,0.]; a1 = [0.,0.]; b1 = [0.,0.]
		flu_rad1 = [0.,0.]; ma_auto1 = [0.,0.]; ell1 = [0.,0.]; kron_r1 = [0.,0.]; C311 = [0.,0.]; rmax1 = [0.,0.]


		
		if setup.FIND_MAJOR_OBJECT == 0:
				a_image = list(a_image)
				iii = a_image.index(max(a_image))

				print '**0**'
				PA = - theta_image[iii]
				xc = x_image[iii]
				yc = y_image[iii]
				RA_c = x_world[iii]
				DEC_c = y_world[iii]
				a = float(a_image[iii])
				b = float(b_image[iii])
				flu_rad = flux_radius50[iii]
				ma_auto = mag_auto[iii]
				ell = ellipticity[iii]
				kron_r = kron_radius[iii]
				C31 = flux_radius75[iii]/flux_radius25[iii]
				rmax = flux_radius99[iii] 			



		if setup.FIND_MAJOR_OBJECT == 1:
			for k in range(n_iter):
				if fabs(RA-x_world[k])< radius_find1 and fabs(DEC-y_world[k])< radius_find1 and a_image[k]*kron_radius[k] > semimajor_find and class_star[k]<0.1:
					PA1.append(- theta_image[k])
					xc1.append(x_image[k])
					yc1.append(y_image[k])
					RA_c1.append(x_world[k])
					DEC_c1.append(y_world[k])
					a1.append(float(a_image[k]))
					b1.append(float(b_image[k]))
					flu_rad1.append(flux_radius50[k])
					ma_auto1.append(mag_auto[k])
					ell1.append(ellipticity[k])
					kron_r1.append(kron_radius[k])
					C311.append(flux_radius75[k]/flux_radius25[k]) 
					rmax1.append(flux_radius99[k])
					number_found = number_found + 1
					print '**1**'



			if number_found==0:
					PA = - theta_image[0]
					xc = x_image[0]
					yc = y_image[0]
					RA_c = x_world[0]
					DEC_c = y_world[0]
					a = float(a_image[0])
					b = float(b_image[0])
					flu_rad = flux_radius50[0]
					ma_auto = mag_auto[0]
					ell = ellipticity[0]
					kron_r = kron_radius[0]
					C31 = flux_radius75[0]/flux_radius25[0]
					rmax = flux_radius99[0] 
					print '**2**'

			elif number_found==1:
					PA = PA1[2]
					xc = xc1[2]
					yc = yc1[2]
					RA_c = RA_c1[2]
					DEC_c = DEC_c1[2]
					a = a1[2]
					b = b1[2]
					flu_rad = flu_rad1[2]
					ma_auto = ma_auto1[2]
					ell = ell1[2]
					kron_r = kron_r1[2]
					C31 = C311[2]
					rmax = rmax1[2]
					print '**3**'

			elif number_found>1:
					iii = a1.index(max(a1))
					PA = PA1[iii]
					xc = xc1[iii]
					yc = yc1[iii]
					RA_c = RA_c1[iii]
					DEC_c = DEC_c1[iii]
					a = a1[iii]
					b = b1[iii]
					flu_rad = flu_rad1[iii]
					ma_auto = ma_auto1[iii]
					ell = ell1[iii]
					kron_r = kron_r1[iii]
					C31 = C311[iii]
					rmax = rmax1[iii]
					print '**4**'
	else:
					PA = - theta_image
					xc = x_image
					yc = y_image
					RA_c = x_world
					DEC_c = y_world
					a = float(a_image)
					b = float(b_image)
					flu_rad = flux_radius50
					ma_auto = mag_auto
					ell = ellipticity
					kron_r = kron_radius
					C31 = flux_radius75/flux_radius25
					rmax = flux_radius99
					print '**5**'


				
	return xc,yc,RA_c,DEC_c,PA,a,b,flu_rad,ma_auto,ell,kron_r,C31,rmax



	
	
