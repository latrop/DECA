#!/usr/bin/python

############## DEComposition Analysis ############## 
##############     Version 1.0.1      ##############          
# Written by Mosenkov Aleksandr as a part of my PhD work. St. Petersburg State University, RUSSIA, 2013.
# It is absolutely for free according to the GNU license but please remember to give proper acknowledgement if you're using this code.
# You can also edit it if you want and email me about that. Let me know if you've found any mistakes in the code so I could fix them.
# My email:	mosenkovAV@gmail.com
####################################################


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
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
import time
import pyfits
import re


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

# We are using this directory where you've created input files
configdir = '.'
sys.path.append(configdir)
DECA_PATH = os.path.dirname(__file__)


#*** DECA MODULES ***
from prep_modules import sextr
from prep_modules import changer
from prep_modules import psf
from prep_modules import image 
from analysis_modules import ell
from output_modules import out
from output_modules import res
from prep_modules import sky_est
import setup
from analysis_modules import soph
from ini_modules import fitting
from ini_modules import fitting_opt_model
from ini_modules import fitting_edge
from ini_modules import fitting_incl
from prep_modules import survey_descr

tmp_out = sys.stdout

if str(sys.argv[1])=='-h':
	# Help for DECA
	print '\n\nUsing:'
	print ' \tpython [path_to_deca]/deca.py [number_of_galaxy] [mode]'
	print 'Arguments:'
	print '\t[number_of_galaxy]	1 --- For one image described in the input file deca_input.py'
	print '				[N] --- For the image [N] from the sample of images given in the input file deca_input.dat\n'
	print '\t[mode]			0 --- Standart running'
	print '				1 --- Only Sextractor running'
	print '				20 --- 1 + Psf determination and building psf image'
	print '				21 --- The same but without building psf image'  
	print '				30 --- 20 + Creation of the final image for Galfit'
	print '				4 --- 30 + Initial guesses'
	print '				5 --- 4 + Sophisticated analysis of edge-on galaxies\n\n'
	exit()



# SALUTATION
print bcolors.OKBLUE+ '########################################################\n' + bcolors.ENDC
print bcolors.OKBLUE+ 'dddddddddd.   eeeeeeeeeeee    cccccc           a  '   + bcolors.ENDC     
print bcolors.OKBLUE+ ' DDD     DDd   EEE       E  cCC    CCc        AAA  '   + bcolors.ENDC     
print bcolors.OKBLUE+ ' DDD      DDD  EEE         CCC               AAAAA  '   + bcolors.ENDC    
print bcolors.OKBLUE+ ' DDD      DDD  EEEeeeeE    CCC              A   AAA  '   + bcolors.ENDC   
print bcolors.OKBLUE+ ' DDD      DDD  EEE         CCC             AAaaaAAAA  '   + bcolors.ENDC 
print bcolors.OKBLUE+ ' DDD     dDD   EEE       e  CCc    ccc    A       AAA  '   + bcolors.ENDC
print bcolors.OKBLUE+ 'dDDDddddDD    eEEEeeeeeeeE   CCccccCC   aAAa     aAAAAa  '   + bcolors.ENDC
print bcolors.OKBLUE+ '\n############# WRITTEN BY MOSENKOV ALEKSANDR #############\n\n' + bcolors.ENDC


#*** Input files ***
if os.path.exists('setup.py')==False:
	print bcolors.OKBLUE+ 'The default setup file is used!' + bcolors.ENDC
else:
	print bcolors.OKBLUE+ 'The user defined setup file is used!' + bcolors.ENDC	

deca_input = 'deca_input.py'
deca_input_sample = 'deca_input.dat'
deca_initials = 'initials.py'
deca_initials_sample = 'initials.dat'

#*** Output files ***
file_ellipse = 'ellipse.png'
file_galfit_outimage = 'galfit_out.fits'
file_dustmod_output = 'dustmod.txt'
file_log = 'out.log'

#*** Other files ***
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


time_begin = time.time() 

# *** FROM SETUP.PY ***
M_psf = setup.M_psf
convolution = setup.convolution
filter_names = setup.filter_names
models = setup.models
file_out_txt = setup.file_out_txt
file_out_ini = setup.file_out_ini
file_out_pdf = setup.file_out_pdf
file_out_sex = setup.file_out_sex
file_ell_txt = setup.file_ell_txt
file_backgr_txt = setup.file_backgr_txt
file_out_sph = setup.file_out_sph
file_out_dust = setup.file_out_dust
ellip_cr = setup.ellip_cr
file_psf_txt = setup.file_psf_txt



#*** INPUT FILES SPECIFIED BY USER ***

mode=0
if os.path.exists('login.cl')==False:
	sys.exit(bcolors.FAIL+ 'WARNING! File login.cl does not exist! Run mkiraf!'+ bcolors.ENDC)

if os.path.exists(deca_input)==True:
	#*** PARAMETERS FROM THE INPUT FILE deca_input ***
	import deca_input as inp
	print bcolors.OKBLUE+ 'The file deca_input.py is used as the input file!' + bcolors.ENDC

 	NUMB = int(sys.argv[1])
	mode = int(sys.argv[2]) 
	
	image_name = inp.image_name
	weight_name = inp.weight_name
	psf_name = inp.psf_name
	survey = inp.survey
	m0 = inp.m0
	xsize = inp.xsize
	ysize = inp.ysize
	fwhm = inp.fwhm

	GAIN = inp.GAIN
	#read_out_noise = inp.read_out_noise
	NCOMBINE = inp.NCOMBINE
	EXPTIME = inp.EXPTIME
	sky_level = inp.sky


	if m0=='NO' or m0=='NONE':
		airmass = inp.airmass	
		aa = inp.aa	
		kk = inp.kk

	RA = inp.coords[0]
	DEC = inp.coords[1]
	if RA=='NO' or RA=='NONE' or RA==0 or DEC=='NO' or DEC=='NONE' or DEC==0:
		sys.exit(bcolors.FAIL + 'The coordinate(s) are not found! Specify \'coords\'!'+ bcolors.ENDC)

	
	filter_name = inp.filter_name
	if filter_name =='NONE' or filter_name == 'NO' or filter_name in filter_names==False:
		sys.exit(bcolors.FAIL + 'Such filter name is not found! Specify \'filter_name\'!'+ bcolors.ENDC)	 
	galaxy_name = inp.galaxy_name
	if galaxy_name == 'NONE':
		galaxy_name = '(%.4f_%.4f)' % (RA,DEC)

	colour_name = inp.colour_name
	colour_value = inp.colour_value
	z = inp.z
	D = inp.D
	if type(colour_name)=='str' or type(colour_value)=='str' or type(z)=='str':
		print bcolors.FAIL+ 'Can\'t find k-coorection and some other cosmological parameters!'+ bcolors.ENDC
	if type(D)=='str' or D==0 or D==0.:
		print bcolors.FAIL+ 'Can\'t find distance to the object. All output values will be observable!'	+ bcolors.ENDC

	Aext = inp.Aext 
	if type(Aext)=='str' or Aext==0 or Aext==0.:
		print bcolors.FAIL+ 'External extinction can\'t be found! All output values will be uncorrected for that!'+ bcolors.ENDC

	image_full = inp.image 
	if image=='NONE':
		sys.exit(bcolors.FAIL+ 'Specify \'image\' with \'field\' or \'full\'!'+ bcolors.ENDC)
	find_sky = inp.find_sky
	find_psf = inp.find_psf	
	model = inp.model
	find_model = inp.find_model

	if model=='NONE' or model in models==False:
		sys.exit(bcolors.FAIL+ 'Such model is not found! Check \'model\'!'+ bcolors.ENDC)	 		
	numb_iter = inp.numb_iter
	sph = inp.sph	 			
	if type(numb_iter)=='str' or numb_iter>3 or numb_iter<1:
		sys.exit(bcolors.FAIL+ 'Problem with \'numb_iter\'! Should be in [0,1,2,3]!'+ bcolors.ENDC)
	ini_use = inp.ini_use

	if ini_use!='YES' and ini_use!='NO' and ini_use!='yes' and ini_use!='no':
		sys.exit(bcolors.FAIL+ 'Problem with \'ini_use\'! Should be \'YES\' or \'NO\'!'+ bcolors.ENDC)

	if ini_use=='YES' or ini_use=='yes':
		if os.path.exists('initials.py')==False:
			sys.exit(bcolors.FAIL+ 'The input file initials.py is not found!'+ bcolors.ENDC)
		else:
			print bcolors.OKBLUE+ 'The program is using your initial parameters from the file initials.py!' + bcolors.ENDC
			import initials

		if model=='exp' or model=='exp+ser' or model=='edge' or model=='edge+ser':
			m0d0 = initials.m0_d
			magDisk0 = initials.magDisk
			h0 = initials.h_d
			z00 = initials.z0
			elld0 = initials.ell_d
			cd0 = initials.c_d
			rtr0 = initials.r_tr
			if m0d0<=0. and magDisk0<=0. or h0<=0. or elld0<0.:
				sys.exit(bcolors.FAIL+ 'Problem with initial(model) parameters! Should be >0.!'+ bcolors.ENDC)			

		if model=='ser' or model=='exp+ser' or model=='edge+ser':
			meb0 = initials.me_bul
			magBul0 = initials.magBul
			reb0 = initials.re_bul
			n0 = initials.n_bul
			ellb0 = initials.ell_bul
			cb0 = initials.c_bul
			if meb0==0 and magBul0==0 or reb0==0 or ellb0==1:
				sys.exit(bcolors.FAIL+ 'Problem with initial(model) parameters! Should be >0.!'+ bcolors.ENDC)

	n_it = 1
	number = 1
	k = number

if os.path.exists(deca_input_sample)==True:
	#*** PARAMETERS FROM THE INPUT FILE deca_input.dat as a list of objects ***
	with open(deca_input_sample, 'r') as f:
		lines = f.readlines()
		num_lines = len([l for l in lines if l.strip(' \n') != ''])
	n_it = num_lines
	NUMB = int(sys.argv[1])
	mode = int(sys.argv[2])

	

if os.path.exists(deca_initials_sample)==True:
	print bcolors.OKBLUE+ 'The program is using your initial parameters from the file initials.dat!' + bcolors.ENDC
	#	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14
	#	#	filter	S0d	mag_d	h_d	z0_d	ell_d	c_d	r_tr	me_bul	mag_bul	re_bul	n_bul	ell_bul	c_bul	

	m0d0s,magDisk0s,h0s,z00s,elld0s,cd0s,rtr0s,meb0s,magBul0s,reb0s,n0s,ellb0s,cb0s = loadtxt('initials.dat', usecols=[2,3,4,5,6,7,8,9,10,11,12,13,14],dtype=float, unpack=True,skiprows=0)
	numbs = loadtxt('initials.dat', usecols=[0],dtype=int, unpack=True,skiprows=0)
	filts = loadtxt('initials.dat', usecols=[1],dtype=str, unpack=True,skiprows=0)



if os.path.exists(deca_input_sample)==True:
	f = open(deca_input_sample,'r')


# Creation of the main ouput txt-files:
#1)
if not os.path.exists("%s" % (file_out_txt)):
	ffff = open(file_out_txt,'w')
	print >>ffff, "#\tm0_d\tm0_d_e\th_d\th_d_e\th_d_k\th_d_k_e\tPA_d\tPA_d_e\tz0\tz0_e\tz0_k\tz0_k_e\trtr\trtr_e\tM_d\tM_d1\tL_d1\tme_bul\tme_bul_e\tre_bul\tre_bul_e\tre_bul_k\tre_bul_k_e\tn_bul\tn_bul_e\tPA_bul\tPA_bul_e\tq_bul\tq_bul_e\tc_bul\tc_bul_e\tM_b\tM_b1\tL_b1\tBT\tBT1\tBD1\tM_g1\tL_g1\tm_g1\tchi2\tmodel\tstatus\tpic_number"
	print >>ffff, '\tmas\tmas\tarcsec\tarcsec\tkpc\tkpc\tdeg\tdeg\tarcsec\tarcsec\tkpc\tkpc\tarcsec\tarcsec\tmag\tmag\tL_sun\tmas\tmas\tarcsec\tarcsec\tkpc\tkpc\t\tdeg\tdeg\t\t\t\t\tmag\tmag\tL_sun\t\t\tmag\tL_sun\tmag'
	ffff.close()	


#2)
if not os.path.exists("%s" % (file_out_sex)):
	fff = open(file_out_sex,'w')
	print >>fff, '#\tRA_c\tDEC_c\txc\tyc\tPA\ta_image\tb_image\tflux_radius\tmag_auto\tellip\tkron_r\tC31\trmax'
	print >>fff, '\tdeg\tdeg\tpix\tpix\tdeg\tarcsec\tarcsec\tarcsec\tmag\t\tsma(smb)\t\tarcsec'
	fff.close()

#3)
if not os.path.exists("%s" % (file_backgr_txt)):
		ff_b = open(file_backgr_txt,'w')
		print >>ff_b, '#\tsky\tvariance\tsky_way\tsky_inp'
		print >>ff_b, '\tDN\tDN^2\t\tDN'
		ff_b.close()

#4)
if not os.path.exists("%s" % (file_psf_txt)):
		ff_psf = open(file_psf_txt,'w')
		print >>ff_psf, '#\tFWHM\tNstars\tpsf_way\tfwhm_inp\tm0'
		print >>ff_psf, '\tarcsec\t\t\tarcsec\tmag arcsec^-2'
		ff_psf.close()

#5)
if not os.path.exists("%s" % (file_ell_txt)):
		ff_ell = open(file_ell_txt,'w')
		print >>ff_ell, '#\tRmax_bar'
		print >>ff_ell, '\tarcsec'
		ff_ell.close()


if n_it>1:
	lines=f.readlines()
	s = re.split('\t',lines[NUMB-1])
	
	number = int(s[0])
	galaxy_name = s[1]
	RA = float(s[2]) 
	DEC = float(s[3]) 
	if galaxy_name == 'NONE':
		galaxy_name = '(%.4f_%.4f)' % (RA,DEC)
	filter_name = s[4]
	colour_name = s[5]
	colour_value = float(s[6])
	z = s[7] #float(z)
	D = s[8]
	if z!='NO' and z!='NONE':
		z = float(z)
	if D!='NO' and D!='NONE':
		D = float(D)
	Aext = float(s[9])
	survey = s[10]
	m0 = float(s[11])
	xsize = float(s[12])
	ysize = float(s[13])
	fwhm = float(s[14])
	GAIN = float(s[15])
	read_out_noise = float(s[16])
	NCOMBINE = int(s[17])
	EXPTIME = float(s[18])
	sky_level = s[19]
	find_sky = int(s[20])
	find_psf = int(s[21])
	image_full = s[22]
	model = s[23]
	find_model = s[24]
	numb_iter = int(s[25])
	sph = s[26]
	ini_use = s[27]
	image_name = s[28]
	weight_name = s[29]
	#psf_name = s[30]
	p_split = re.compile(r"\n")
	psf_name = str(p_split.split(s[30])[0])

	k = NUMB


if n_it>0:

	#############################################################################################################################
	#############################################################################################################################
	#############################################################################################################################
	def main():
		global EXPTIME, GAIN, NCOMBINE,pix2sec,m0,sky_level,weight_name,fwhm,model,find_model,Aext,file_out_pdf,k,read_out_noise

		# Checking the existence of files:
		if image_name.find('.fits')>=1 or image_name.find('.fit')>=1 or image_name.find('.fits.gz')>=1 or image_name.find('.fit.gz')>=1 or image_name.find('.fits.bz2')>=1 or image_name.find('.fit.bz2')>=1:
			file_image = image_name
		else:
			sys.exit(bcolors.FAIL+ 'The format of the input file is incorrect! The image should have .fits extension and the list of objects should have .list!'+ bcolors.ENDC)	

		if psf_name.find('.fits')>=1 or psf_name.find('.fit')>=1:
			psf_image = psf_name
			print psf_image, file_image
			if os.path.exists(psf_image)==False:
				sys.exit(bcolors.FAIL+ 'File %s does not exist!' % (psf_image) + bcolors.ENDC)
		else:
			psf_image = psf_name

		if weight_name.find('.fits')>=1 or weight_name.find('.fit')>=1:
			weight_image = weight_name
			if os.path.exists(weight_image)==False:
				sys.exit(bcolors.FAIL+ 'File %s does not exist!'+ bcolors.ENDC % (weight_image))

		if psf_name=='NONE' or psf_name=='NO':
			psf_image = psf_name


		if os.path.exists(file_image)==False:
			sys.exit(bcolors.FAIL+ 'File %s does not exist!'+ bcolors.ENDC % (file_image))

		ff = open('out.log', 'a') # Out to the out.log

		print bcolors.OKGREEN+ '\n************ OBJECT %i ************' % (k) + bcolors.ENDC
		# STEP 1: Copying the image to the current directory
		print bcolors.OKBLUE+ '1. Copying the input image ...' + bcolors.ENDC
		print >> ff, '**** The object number: %i ****' % (number)
		print >> ff, '1. Copying the input image %s to the current directory.' % (file_image)
		print file_image

		if file_image.find('.fits.gz')>=1 or file_image.find('.fit.gz')>=1:
			print ' Unzipping ...'
			print >> ff, ' Unzipping ...'
			shutil.copy(file_image,file_field+'.gz')
			subprocess.call("gzip -d %s" % (file_field+'.gz'), shell=True)
		elif file_image.find('.fits.bz2')>=1 or file_image.find('.fit.bz2')>=1:
			print ' Unzipping ...'
			print >> ff, ' Unzipping ...'
			shutil.copy(file_image,file_field+'.bz2')
			subprocess.call("bzip2 -d %s" % (file_field+'.bz2'), shell=True)
		else:
			shutil.copy(file_image,file_field)
		
		#exit()

		# CCD parameters of the surveys: SDSS_dr7, SDSS_dr8, UKIDSS, 2MASS, EFIGI, WISE, Spitzer 
		if survey=='NO' or survey=='NONE':
			pix2sec = xsize

			if pix2sec=='NO' or pix2sec=='NONE':
				sys.exit(bcolors.FAIL+ 'The pixel scale is not found! Check \'xsize\' and \'ysize\'!'+ bcolors.ENDC)
			if m0=='NO' or m0=='NONE':
				sys.exit(bcolors.FAIL+ 'The zeropint m0 cannot be found!')
		else:
			if survey == 'EFIGI':
				nx,ny,pix2sec,GAIN,read_out_noise,sky_level1,fwhm1,m0,Aext1,EXPTIME,NCOMBINE = survey_descr.EFIGI(file_field,aa,kk,airmass,GAIN)
			elif survey == 'SDSS_dr7':
				nx,ny,pix2sec,GAIN,read_out_noise,sky_level1,fwhm1,m0,Aext1,EXPTIME,NCOMBINE = survey_descr.SDSS_dr7(file_field,aa,kk,airmass)
			elif survey == 'SDSS_dr8':
				nx,ny,pix2sec,GAIN,read_out_noise,sky_level1,fwhm1,m0,Aext1,EXPTIME,NCOMBINE = survey_descr.SDSS_dr8(file_field)
			elif survey == '2MASS':
				nx,ny,pix2sec,GAIN,read_out_noise,sky_level1,fwhm1,m0,Aext1,EXPTIME,NCOMBINE = survey_descr.MASS2(file_field)
			elif survey == 'UKIDSS':
				nx,ny,pix2sec,GAIN,read_out_noise,sky_level1,fwhm1,m0,Aext1,EXPTIME,NCOMBINE = survey_descr.UKIDSS(file_field)
			elif survey == 'WISE':
				nx,ny,pix2sec,GAIN,read_out_noise,sky_level1,fwhm1,m0,Aext1,EXPTIME,NCOMBINE = survey_descr.WISE(file_field)
			elif survey == 'Spitzer':
				nx,ny,pix2sec,GAIN,read_out_noise,sky_level1,fwhm1,m0,Aext1,EXPTIME,NCOMBINE = survey_descr.Spitzer(file_field)

			if Aext1!='not_found':
				Aext = Aext1

			if sky_level1!='not_found':
				sky_level = sky_level1

			if fwhm1!='not_found':
				fwhm = fwhm1


		print 'pix2sec= %.3f arcsec/pix' % (pix2sec)
		print 'm0 = %.3f mag arcsec^-2' % (m0)
		print >>ff, 'pix2sec= %.3f arcsec/pix' % (pix2sec)
		print >>ff, 'm0 = %.3f mag arcsec^-2' % (m0)




		#*****************************************************************************************
		#******************************* THE PROGRAM ITSELF **************************************



		#********************************PREPARATION*********************************
		'''
		print bcolors.OKGREEN+ '\n************ OBJECT %i ************' % (k) + bcolors.ENDC
		# STEP 1: Copying the image to the current directory
		print bcolors.OKBLUE+ '1. Copying the input image ...' + bcolors.ENDC
		print >> ff, '**** The object number: %i ****' % (number)
		print >> ff, '1. Copying the input image %s to the current directory.' % (file_image)
		

		if file_image.find('.fits.gz')>=1 or file_field.find('.fit.gz')>=1:
			shutil.copy(file_image,file_field+'.gz')
			subprocess.call("gzip -d %s" % (file_field+'.gz'), shell=True)
		elif file_image.find('.fits.bz2')>=1 or file_field.find('.fit.bz2')>=1:
			shutil.copy(file_image,file_field+'.bz2')
			subprocess.call("bzip2 -d %s" % (file_field+'.bz2'), shell=True)
		else:
			shutil.copy(file_image,file_field)
		
		'''
		if survey=='NO' or survey=='NONE':
			hdulist = pyfits.open(file_field)
			prihdr = hdulist[0].header
			nx = prihdr['NAXIS1']
			ny = prihdr['NAXIS2']


		# STEP 2-0: Image preparation  (Surveys only)
		if survey=='SDSS_dr7':
			print bcolors.OKBLUE+ '2. Subtracting soft bias ...' + bcolors.ENDC
			print >> ff, '2. Subtracting soft bias ...'
			try:
				from pyraf import iraf
				iraf.imarith(file_field,'-',1000.,file_field)	# Subtracting 1000 DN of soft bias
			except:
				f_res = open("resid.cl", "w") 
				sys.stdout = f_res
				print "# Preparation of SDSS dr7 frames"
				print "images"
				print "imutil"
				print "imarith %s - 1000.0 %s" % (file_field,file_field)
				print "logout"	
				sys.stdout = tmp_out
				f_res.close()
				os.chmod(r"resid.cl",0777)
				subprocess.call("cl < resid.cl -o", shell=True)
				os.remove('resid.cl')
		if survey=='SDSS_dr8':
			#darkVariance = (read_out_noise/GAIN)**2
			print bcolors.OKBLUE+ '2. Converting nMgy to standart DN ...' + bcolors.ENDC
			print >> ff, '2. Converting nMgy to standart DN ...'
			image.prep_ima(survey,file_field,file_field1, GAIN)
			os.rename(file_field1,file_field)


		# STEP 3-1: psf.fits copying if exists
		if psf_image.find('.fits')>=1 or psf_image.find('.fit')>=1:
			print bcolors.OKBLUE+ '3. Copying the input psf image ...' + bcolors.ENDC
			print >> ff, '3. Copying the input psf image ...'
			shutil.copy(psf_image,file_psf)

		# STEP 3-2: weight.fits copying if exists
		if weight_name.find('.fits')>=1 or weight_name.find('.fit')>=1:
			print bcolors.OKBLUE+ '  Copying the input sigma image ...' + bcolors.ENDC
			print >> ff, '  Copying the input sigma image ...'
			shutil.copy(weight_name,file_weight)
		else:
			file_weight = 'none'
		#print survey	
		#exit()

		if survey=='SDSS_dr8':
			print file_field2
			sky_est.imcopy_func('field.fits[0]',file_field2,1,1,nx,ny)
		elif survey=='UKIDSS':
			sky_est.imcopy_func('field.fits[1]',file_field2,1,1,nx,ny)
			hdulist1 = pyfits.open(file_field2, do_not_scale_image_data=True, mode='update')
			prihdr1 = hdulist1[0].header
			prihdr1.update('CTYPE1','RA---TAN')
			prihdr1.update('CTYPE2','DEC--TAN')
			hdulist1.flush()			
		else:
			#try:
			#	sky_est.imcopy_func('field[0].fits',file_field2,1,1,nx,ny)	# [0] is added
			#	sky_est.imcopy_func('field[0].fits','field3.fits',1,1,nx,ny)	# ADDED
			#	shutil.move('field3.fits',file_field)
			#except:
			try:
				sky_est.imcopy_func('field.fits[0]',file_field2,1,1,nx,ny)	
				sky_est.imcopy_func('field.fits[0]','field3.fits',1,1,nx,ny)	
			except:
				sky_est.imcopy_func('field.fits',file_field2,1,1,nx,ny)	
				sky_est.imcopy_func('field.fits','field3.fits',1,1,nx,ny)	
							

		#exit()
		# STEP 4: First time SExtractor runing: for the field image
		print bcolors.OKBLUE+ '4. SExtractor runing for the first time ...' + bcolors.ENDC
		print >> ff, '4. SExtractor runing for the first time ...'
		if NCOMBINE>1:
			sextr.run_sextr(survey,file_field,m0,GAIN*NCOMBINE,pix2sec,fwhm,DECA_PATH)
		else:
			sextr.run_sextr(survey,file_field,m0,GAIN,pix2sec,fwhm,DECA_PATH)

		print 'We are searching for an object: %.6f,%.6f' % (RA,DEC) 
		xc,yc,RA_c,DEC_c,PA,a_image,b_image,flux_radius,mag_auto,ellip,kron_r,C31,rmax = sextr.read_sextr_find(RA,DEC,file_sex_field)
		RA1=RA_c; DEC1=DEC_c; PA1 = PA
		print 'The coordinates of the object from SExtractor: (RA,DEC) =  %.4f,%.4f' % (RA1,DEC1)
		print 'Positional angle from SExtractor: PA =  %.1f' % (PA)
		print >>ff, '  The coordinates of the object: (RA,DEC) =  %.4f,%.4f' % (RA1,DEC1)
		print >>ff, '  Positional angle: PA =  %.1f' % (PA)

		#%%%%%%%%		
		#find_sky=0
		#%%%%%%%%

		#exit()



		# NOT TESTING
		way_backr = setup.sky_find
		sky_level0 = sky_level
		# STEP 5: THE BACKGROUND ESTIMATION
		if find_sky==1 or type(sky_level)=='str':
			#way_backr = setup.sky_find 
			print bcolors.OKBLUE+ '5. The background estimation ...' + bcolors.ENDC
			print >> ff, '5. The background estimation ...'
			#if survey=='2MASS':
			#	way_backr = 4
			
			if max([nx/2.,ny/2.])>2.*kron_r*a_image and math.pi*(kron_r*a_image)**2<nx*ny/3.:
				sky_level,variance = sky_est.sky_estim(survey,xc,yc,kron_r*a_image,ellip,PA,way_backr)
			else:
				way_backr_altern = 4
				way_backr = way_backr_altern
				sky_level,variance = sky_est.sky_estim(survey,xc,yc,kron_r*a_image,ellip,PA,way_backr_altern)
				print 'There is not enough sky backround to estimate it properly! The measured sky_level and variance may appear overesimated!' 
			print 'SKY_LEVEL = %.3f' % (sky_level)
			print >> ff, '  SKY_LEVEL = %.3f' % (sky_level)
			print 'Variance = %.3f' % (variance)
			print >> ff, '  Variance = %.3f' % (variance)
		elif find_sky!=1 and type(sky_level)!='str':
			variance = (read_out_noise/GAIN)**2
			#sky_level4,variance = sky_est.sky_estim(survey,xc,yc,kron_r*a_image,ellip,PA,way_backr)
		else:
			sys.exit(bcolors.FAIL+ 'The sky_level is not given and find_sky is equal to 0 (should be 0 for a given sky_level)!')

		#print type(sky_level,variance,k
		ff_b = open(file_backgr_txt,'a')
		print >>ff_b, '%i\t%.3f\t%.3f\t%i\t%s' % (k,float(sky_level),variance,way_backr,str(sky_level0))
		ff_b.close()

		if mode==1:
			exit()




		'''
		# THE BACKGROUND ESTIMATION: DIFFIRENT WAYS
		sky_level1 = []
		variance1 = []
		for way_backr in range(4,6,1):
			if find_sky==1 or type(sky_level)=='str':
				print bcolors.OKBLUE+ 'The initial background estimation ...' + bcolors.ENDC
				print >> ff, 'The initial background estimation ...'
				sky_level2,variance2 = sky_est.sky_estim(survey,xc,yc,kron_r*a_image,ellip,PA,way_backr)
				print 'SKY_LEVEL = %.3f' % (sky_level2)
				print >> ff, 'SKY_LEVEL = %.3f' % (sky_level2)
				print 'Variance = %.3f' % (variance2)
				print >> ff, 'Variance = %.3f' % (variance2)
				sky_level1.append(sky_level2)
				variance1.append(variance2)
			else:
				variance = (read_out_noise/GAIN)**2


		ff_b = open(file_backgr_txt,'a')

		#print >>ff_b, '%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' % (k,sky_level1[0],variance1[0],sky_level1[1],variance1[1],sky_level1[2],variance1[2],sky_level1[3],variance1[3],sky_level1[4],variance1[4],sky_level1[5],variance1[5],sky_level,(read_out_noise/GAIN)**2)

		print >>ff_b, '%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' % (k,float(sky_level1[0]),variance1[0],float(sky_level1[1]),variance1[1],float(sky_level),(read_out_noise/GAIN)**2)

		ff_b.close()
		
		if mode==1:
			exit()
		'''

		

		#%%%%%%%%
		#find_psf = 0
		#%%%%%%%%

		#print sky_level
		#exit()
		# Subtraction of the sky level
		try:
			from pyraf import iraf
			iraf.imarith(file_field2,'-',sky_level,file_field2)
			#iraf.imarith(file_field,'-',sky_level,file_field)
		except:
			#print file_field2,sky_level,file_field2
			#exit()

			f_res = open("resid.cl", "w") 
			sys.stdout = f_res
			print "Script: esdid"
			print "images"
			print "imutil"
			print "imarith %s - %.3f %s" % (file_field2,float(sky_level),file_field2)
			#print "imarith %s - %.3f %s" % (file_field,sky_level,file_field)
			print "logout"	
			sys.stdout = tmp_out
			f_res.close()
			os.chmod(r"resid.cl",0777)
			subprocess.call("cl < resid.cl -o", shell=True)
			os.remove('resid.cl')


		# STEP 6: PSF model
		if find_psf==0:
			Nstars = 1
			FWHM = fwhm
			print 'Given PSF image will be used!'		
		if find_psf>0:
			print bcolors.OKBLUE+ '6. psf.fits creating...' + bcolors.ENDC
			print >> ff, '6. psf.fits creating...'		
		if find_psf == 1:
			Nstars = 1
			psf.mod_psf(fwhm/pix2sec,M_psf,0.,0.,m0 - 2.5*log10(EXPTIME),pix2sec,pix2sec)
			FWHM = fwhm
		if find_psf == 2:
			try:			
				FWHM,Nstars = psf.PSF_model_sph(survey,file_field2,m0,EXPTIME,pix2sec,pix2sec,GAIN,sky_level)
			except:
				psf.mod_psf(fwhm/pix2sec,M_psf,0.,0.,m0 - 2.5*log10(EXPTIME),pix2sec,pix2sec)
				FWHM = fwhm
		if find_psf == 3:
			try:
				FWHM,Nstars = psf.PSF_model_best(survey,file_field2,m0,EXPTIME,pix2sec,pix2sec,GAIN,sky_level)		
			except:
				print 'Not enough stars with good S/N!'
				psf.mod_psf(fwhm/pix2sec,M_psf,0.,0.,m0 - 2.5*log10(EXPTIME),pix2sec,pix2sec)
				FWHM = fwhm
				Nstars = 0
		if find_psf == 4:
				try:
					FWHM,Nstars = psf.PSF_star_best(survey,file_field2,m0,pix2sec,pix2sec,GAIN,sky_level,EXPTIME)
				except:
					psf.mod_psf(fwhm/pix2sec,M_psf,0.,0.,m0 - 2.5*log10(EXPTIME),pix2sec,pix2sec)
					FWHM = fwhm
			
		if find_psf == 5: #and (survey=='SDSS_dr7' or survey=='SDSS_dr8'):
				Nstars = 1
				#try:
				psf.PSF_extract(xc,yc,filter_name)

				FWHM = fwhm
				print 'psField file from SDSS was used to extract the psf image!'
				#except:
				#psf.mod_psf(fwhm/pix2sec,M_psf,0.,0.,m0 - 2.5*log10(EXPTIME),pix2sec,pix2sec)
				#FWHM = fwhm	
		if find_psf == 10:
			Nstars = 0
			FWHM = fwhm
			print 'There is no convolution!!!'	

		print '  FWHM = %.3f [arcsec]' % (FWHM) 
		print >> ff, '  FWHM = %.3f [arcsec]' % (FWHM) 
		#exit()

		ff_psf = open(file_psf_txt,'a')
		print >>ff_psf, '%i\t%.3f\t%i\t%i\t%.3f\t%.3f' % (k,FWHM,Nstars,find_psf,fwhm,m0)
		ff_psf.close()
		
		if mode==20 or mode==21:
			if mode==20:
				if not os.path.exists("./pics/%i" % (k)):
					os.makedirs("./pics/%i" % (k))
				path ='./pics/' + str(k) + '/'
				shutil.move('psf.fits',path+'psf.fits')
				exit()
			else:
				exit()

		if survey=='UKIDSS':
			hdulist = pyfits.open(file_field, do_not_scale_image_data=True, mode='update')
			prihdr = hdulist[1].header
			prihdr['CTYPE1'] = 'RA'
			prihdr['CTYPE2'] = 'DEC'
			prihdr.update('SYSTEM','FK5')
			hdulist.flush()

		#%%%%%%%%
		#image_full='field'
		#%%%%%%%%
		#exit()

		# STEP 7: Field image rotating and the galaxy image extracting
		if image_full=='full':
			os.rename(file_field,file_gal)
		if image_full=='field':
			if model=='edge' or model=='edge+ser':
				# Rotation of the image about (xc,yc), angle=PA
				print bcolors.OKBLUE+ '7. Preparation of the image ...' + bcolors.ENDC
				print >> ff, '7.  Preparation of image ...'
				print >> ff, '  The rotation angle for IRAF: PA=%.2f deg.' % (-PA)	
				changer.rotate(survey,file_field,file_field1,-PA,xc,yc,nx,ny)
				os.remove(file_field)
				os.rename(file_field1,file_field)
				print RA_c,DEC_c
				#exit()
				#if survey!='UKIDSS':
				#	xc,yc=changer.coord_convert(file_field,RA_c,DEC_c)
				#print xc,yc
			 	# Extraction of the object image
				
				print bcolors.OKBLUE+ 'Extracting the image ...' + bcolors.ENDC
				print >> ff, '  Extracting the image ...'	
				changer.extract1(survey,file_field,file_gal,nx/2,ny/2,1.5*rmax,1.5*rmax*(1.-ellip),0.)
				os.remove(file_field)
				#exit()
			else:
				# STEP: Extraction of the object image
				# Rotation of the image about (xc,yc), angle=PA
				print bcolors.OKBLUE+ '7. Preparation of the image ...' + bcolors.ENDC
				print >> ff, '7.  Preparation of image ...'
				print >> ff, '  The rotation angle  for IRAF: PA=%.2f degrees.' % (-PA)	
				try:
					changer.rotate(survey,file_field+'[0]',file_field1,-PA,xc,yc,nx,ny)
				except:
					changer.rotate(survey,file_field,file_field1,PA,xc,yc,nx,ny)					
				os.remove(file_field)
				os.rename(file_field1,file_field)
				if survey=='UKIDSS':
					sextr.run_sextr(survey,file_field,m0,GAIN,pix2sec,fwhm,DECA_PATH)
					xc,yc,RA_c,DEC_c,PA,a_image,b_image,flux_radius,mag_auto,ellip,kron_r,C31,rmax = sextr.read_sextr_find(RA,DEC,file_sex_field)
				#else:
				#	xc,yc=changer.coord_convert(file_field,RA_c,DEC_c)

				print bcolors.OKBLUE+ 'Extracting the image ...' + bcolors.ENDC
				print >> ff, '  Extracting the image ...'	
				extr_coeff = setup.extr_coeff
				changer.extract1(survey,file_field,file_gal,nx/2,ny/2,1.5*rmax,1.5*rmax*(1.-ellip),PA)
				os.remove(file_field)

				hdulist = pyfits.open(file_gal)
				prihdr = hdulist[0].header
				nxx = prihdr['NAXIS1']
				nyy = prihdr['NAXIS2']
				#print nxx/2.,nyy/2.,2.*rmax
				#nxx=398
				#nyy=318
				#print int(nx/2-1.5*rmax),int(ny/2-1.5*rmax*(1.-ellip)),int(nx/2+1.5*rmax),int(ny/2+1.5*rmax*(1.-ellip))
				#exit()
				
				try:
					'''
					if nxx/2.>3.*rmax:
						sky_est.imcopy_func(file_gal,'galaxy1.fits',int(ceil(nxx/2-3.*rmax)),int(floor(nyy/2-3.*rmax*(1.-ellip))),int(ceil(nxx/2+3.*rmax)),int(floor(nyy/2+3.*rmax*(1.-ellip))))
						shutil.copy('galaxy1.fits',file_gal)
						os.remove('galaxy1.fits')
						print '1!'
					else:
						print '2!'
					'''
					AA = 2
					if nxx/2.>AA*rmax+2:
						sky_est.imcopy_func(file_gal,'galaxy1.fits',int(ceil(nxx/2-AA*rmax)),int(floor(nyy/2-AA*rmax*(1.-ellip))),int(ceil(nxx/2+AA*rmax)),int(floor(nyy/2+AA*rmax*(1.-ellip))))
						shutil.copy('galaxy1.fits',file_gal)
						os.remove('galaxy1.fits')
						print '1!'
					else:
						print '2!'
				except:
					print 'The cutted image is ok!'
				
				#exit()




			os.remove(file_sex_field)
			os.remove(file_segm)
			#os.remove(file_badpix)
			os.remove('aper.fits')
			os.remove('objects.fits')
			os.remove('fon.fits')
			#os.remove('modelPSF.txt')
		#iraf.imarith(file_gal,'-',sky_level,file_gal)
		
		#exit()
		# STEP 8: Second time SExtractor runing: for galaxy image
		print bcolors.OKBLUE+ '8. Second time SExtractor runing ...' + bcolors.ENDC
		print >> ff, '8. Second time SExtractor runing ...'
		if NCOMBINE>1:
			sextr.run_sextr(survey,file_gal,m0,GAIN*NCOMBINE,pix2sec,fwhm,DECA_PATH)
		else:
			sextr.run_sextr(survey,file_gal,m0,GAIN,pix2sec,fwhm,DECA_PATH)
		hdulist = pyfits.open(file_gal)
		scidata = hdulist[0].data
		ny,nx = np.shape(scidata)
		xc,yc,RA_c,DEC_c,PA,a_image,b_image,flux_radius,mag_auto,ellip,kron_r,C31,rmax = sextr.read_sextr_find((nx/2),(ny/2),file_sex_field,coords=2)	#earlier was RA,DEC
		# all scales are in pixels and SB are in mag/arcsec^2
		rmax = fabs(rmax)
		

		fff = open(file_out_sex,'a')
		print >> fff, '%i\t%.6f\t%.6f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' % (k,RA_c,DEC_c,xc,yc,PA,a_image*pix2sec,b_image*pix2sec,flux_radius*pix2sec,mag_auto,ellip,kron_r,C31,rmax*pix2sec)
		fff.close()

		#print mag_auto, m0,xc,yc,RA_c,DEC_c
		
		

		#exit()


		# STEP 9: Bad pixel mask and cleaning galaxy image
		print bcolors.OKBLUE+ '9. Bad pixel mask creating ...' + bcolors.ENDC
		print >> ff, '9. Bad pixel mask creating ...'
		if setup.mask_stars==1:
			sextr.bad_pix_mask(file_segm,file_badpix,xc,yc)
		else:
			xcc,ycc,cl_star,rmax_star = loadtxt(file_sex_field, usecols=[1,2,20,12], unpack=True, skiprows = 18)
			sextr.bad_pix_mask_stars(file_segm,file_badpix,xcc,ycc,cl_star,rmax_star)
		print bcolors.OKBLUE+ 'Cleaning the galaxy image from contaminents ...' + bcolors.ENDC
		print >> ff, '  Cleaning the galaxy image from contaminents ...'



		# Subtraction of the sky level
		try:
			from pyraf import iraf
			iraf.imarith(file_gal,'-',sky_level,file_gal)	# Subtracting 1000 DN of soft bias
		except:
			f_res = open("resid.cl", "w") 
			sys.stdout = f_res
			print "Script: resid"
			print "images"
			print "imutil"
			print "imarith %s - %.3f %s" % (file_gal,float(sky_level),file_gal)
			print "logout"	
			sys.stdout = tmp_out
			f_res.close()
			os.chmod(r"resid.cl",0777)
			subprocess.call("cl < resid.cl -o", shell=True)
			os.remove('resid.cl')
		

		if setup.mask_stars==1:
			sextr.gal_clean_ima(file_segm,file_gal,file_gal_clean,xc,yc,float(sky_level))
		else:
			
			sextr.gal_clean_ima_stars(file_segm,file_gal,file_gal_clean,xcc,ycc,cl_star,rmax_star)

		if mode==30:
			if not os.path.exists("./pics/%i" % (k)):
				os.makedirs("./pics/%i" % (k))
			path ='./pics/' + str(k) + '/'
			shutil.move(file_gal_clean,path+file_gal_clean)
			exit()
		elif mode==31:
			exit()
		#exit()

		#print xc,yc
		if mode!=4:
			#STEP 10: Ellipse fitting
			#rmax=a_image*kron_r
			try:
				#if model=='exp' or model=='exp+ser' or model=='ser' or model=='edge' or model=='edge+ser':
				step = 1.0
				minsma = 1.0
				maxsma = rmax + rmax/2.
				print bcolors.OKBLUE+ '10. Ellipse fitting ...' + bcolors.ENDC
				print >> ff, '10. Ellipse fitting ...'
				ell.main_ell(file_gal_clean,xc,yc,step,minsma,maxsma,m0,pix2sec,FWHM,rmax,sqrt(variance))
				rmaxBAR = ell.rmax_bar(rmax)


				ff_ell = open(file_ell_txt,'a')
				print >>ff_ell, '%i\t%.3f' % (k,rmaxBAR*pix2sec)
				ff_ell.close()				

			except:
				print '10. Ellipse fitting failed!!!'
				print >> ff, '10. Ellipse fitting failed!'

		#exit()


		##############################&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&







		SKY_LEVEL = 0.
		NOISE = variance	# NOISE=(std.dev)^2

		#exit()
		
		# STEP 11: Updating the header with GAIN, NCOMBINE and EXPTIME keywords. M0 (zeropoint) is found.
		print bcolors.OKBLUE+ '11. Updating the image descriptor: ' + bcolors.ENDC
		print >> ff, '11. Updating the image descriptor: GAIN=%.3f, EXPTIME=%.3f, NCOMBINE=%i' % (GAIN,EXPTIME,NCOMBINE)
		image.add_to_header(file_gal_clean,EXPTIME,m0,GAIN,NCOMBINE,pix2sec,NOISE)


		#%%%%%%%%%%%%%
		#find_model='opt'
		#find_model='NO'
		#find_model='edge'
		#%%%%%%%%%%%%%

		
		if setup.fit_sky==1:
			sky_level = 99999.
		else:
			sky_level = 0.

		#print 'sky',sky_level,find_model
		#exit()
		# STEP 12: INITIAL GUESSES AND FIRST RUN OF GALFIT
		print bcolors.OKBLUE+ '12. Initial guesses finding ...' + bcolors.ENDC
		print >> ff, '12. Initial guesses finding  and the GALFIT runing for the first time ...'
		if find_model=='opt':
			numb,xc_d,yc_d,m0_d,h,q_d,z0,rtr,xc_bul,yc_bul,me_bul,re_bul,q_bul,n_bul,c_bul,chi2,model,status = fitting_opt_model.main(ff,number,xc,yc,kron_r,flux_radius,mag_auto,a_image,b_image,PA,sqrt(NOISE),pix2sec,m0,FWHM,EXPTIME,ini_use,find_model,model,ellip,C31,sky_level,file_out_txt,file_out_pdf,filter_name,colour_name,colour_value,Aext,z,D,galaxy_name,find_psf,rmax,mode)
		elif find_model=='edge':
			numb,xc_d,yc_d,m0_d,h,q_d,z0,rtr,xc_bul,yc_bul,me_bul,re_bul,q_bul,n_bul,c_bul,chi2,model,status = fitting_edge.main(ff,number,xc,yc,kron_r,flux_radius,mag_auto,a_image,b_image,PA,sqrt(NOISE),pix2sec,m0,FWHM,EXPTIME,ini_use,find_model,model,ellip,C31,sky_level,file_out_txt,file_out_pdf,filter_name,colour_name,colour_value,Aext,z,D,galaxy_name,find_psf,rmax,mode)
		elif find_model=='incl':
			numb,xc_d,yc_d,m0_d,h,q_d,z0,rtr,xc_bul,yc_bul,me_bul,re_bul,q_bul,n_bul,c_bul,chi2,model,status = fitting_incl.main(ff,number,xc,yc,kron_r,flux_radius,mag_auto,a_image,b_image,PA,sqrt(NOISE),pix2sec,m0,FWHM,EXPTIME,ini_use,find_model,model,ellip,C31,sky_level,file_out_txt,file_out_pdf,filter_name,colour_name,colour_value,Aext,z,D,galaxy_name,find_psf,rmax,mode)
		else:
			numb,xc_d,yc_d,m0_d,h,q_d,z0,rtr,xc_bul,yc_bul,me_bul,re_bul,q_bul,n_bul,c_bul,chi2,model,status = fitting.main(ff,number,xc,yc,kron_r,flux_radius,mag_auto,a_image,b_image,PA,sqrt(NOISE),pix2sec,m0,FWHM,EXPTIME,ini_use,find_model,model,ellip,C31,sky_level,file_out_txt,file_out_pdf,filter_name,colour_name,colour_value,Aext,z,D,galaxy_name,find_psf,rmax,mode)

		
		#numb_iter = 1
		if int(status)==0:
			try:
				if numb_iter>=2 and find_model=='edge':
					sub_out = './pics/'+str(number)+'/out_'+str(numb)+'.fits'		
					shutil.copy(sub_out,'out.fits')
					#import mask_contam
					#mask_contam.mask('out.fits',xc_d,yc_d,rmax,rmax*(1.-ellip))	
				
					#exit()
					# STEP 14: GALFIT: SECOND TIME RUNING
					print bcolors.OKBLUE+ 'Galfit runing for the second time ...' + bcolors.ENDC
					print >> ff, 'Galfit runing for the second time ...'

					numb,xc_d,yc_d,m0_d,h,q_d,z0,rtr,xc_bul,yc_bul,me_bul,re_bul,q_bul,n_bul,c_bul,chi2,model,status = fitting_edge.second_iter(number,model,numb,nx,ny,xc,yc,pix2sec,m0,EXPTIME,sky_level,file_out_txt,file_out_pdf,filter_name,colour_name,colour_value,Aext,z,D,galaxy_name,find_psf,m0_d,h,z0,99999.,rtr,me_bul,re_bul,n_bul,1.-q_bul,c_bul,chi2)
			except:
				print 'The Second time running crashed!'

			#try:
			if numb_iter>=2 and find_model=='incl':
					print 'Galfit runing for the second time ...'
					sub_out = './pics/'+str(number)+'/out_'+str(numb)+'.fits'		
					shutil.copy(sub_out,'out.fits')
					#import mask_contam
					#mask_contam.mask('out.fits',xc_d,yc_d,rmax,rmax*(1.-ellip))	
				
					#exit()
					# STEP 14: GALFIT: SECOND TIME RUNING
					print bcolors.OKBLUE+ 'Galfit runing for the second time ...' + bcolors.ENDC
					print >> ff, 'Galfit runing for the second time ...'

					numb,xc_d,yc_d,m0_d,h,q_d,z0,rtr,xc_bul,yc_bul,me_bul,re_bul,q_bul,n_bul,c_bul,chi2,model,status = fitting_incl.second_iter(number,model,numb,nx,ny,xc,yc,pix2sec,m0,EXPTIME,sky_level,file_out_txt,file_out_pdf,filter_name,colour_name,colour_value,Aext,z,D,galaxy_name,find_psf,m0_d,h,z0,99999.,rtr,me_bul,re_bul,n_bul,1.-q_bul,c_bul,chi2,FWHM)
			#except:
			#	print 'The Second time running crashed!'				
			#exit()
			PA = PA+90.
			# STEP 13: Plots
			if numb_iter<=2:
				image_out = './pics/'+str(number)+'/out_'+str(numb)+'.fits'
				sub_out = './pics/'+str(number)+'/sub_'+str(numb)+'.fits'		
				shutil.copy(image_out,file_galfit_outimage)
				shutil.copy(sub_out,'subcomps.fits')		

			m_crit = -2.5*log10(sqrt(NOISE)) + m0 + 5.*log10(pix2sec)
			print bcolors.OKBLUE+ 'Picture creating ...' + bcolors.ENDC
			print >> ff, 'Picture creating ...'
			if model=='exp' or model=='exp+ser' or model=='ser':
				try:
					step = 1.0
					minsma = 1.0
					maxsma = rmax + rmax/2.
					out.profiles_galfit_extract(xc,yc,PA,model,z0=0.,h=0.)
					out.major_prof(m0,pix2sec,rmax,m_crit,model)
					out.profiles_galfit_extract(xc,yc,PA-90.,model,z0=0.,h=0.)
					out.minor_prof(m0,pix2sec,rmax*(1.-ellip),m_crit,model)
					out.azim_prof(xc,yc,step,minsma,maxsma,m0,pix2sec)
					out.radius_prof(xc,yc,m0,pix2sec,rmax)
				
					if model=='exp+ser':
						h1 = h
						out.isophotes(number,file_gal_clean,xc,yc,m0,pix2sec,m_crit,rmax,rmax*(1.-ellip),h=h1,z0=0.,re=0.)
					if model=='exp':
						h1 = h
						m_crit = m0_d + 4.
						out.isophotes(number,file_gal_clean,xc,yc,m0,pix2sec,m_crit,rmax,rmax*(1.-ellip),h=h1,z0=0.,re=0.)
					if model=='ser':
						re1 = re_bul
						out.isophotes(number,file_gal_clean,xc,yc,m0,pix2sec,m_crit,rmax,rmax*(1.-ellip),h=0,z0=0.,re=re1)
				except:
					print bcolors.FAIL+ 'THERE IS A PROBLEM WITH ISOPHOTE PLOTTING!'+ bcolors.ENDC

			if (model=='edge' or model=='edge+ser') and setup.warps_cut==0 and sph=='NO':
				sub_out = './pics/'+str(number)+'/out_'+str(numb)+'.fits'
				hdulist = pyfits.open(sub_out)
				prihdr = hdulist[1].header
				nx1 = prihdr['NAXIS1']
				ny1 = prihdr['NAXIS2']	
				os.remove('galaxy_clean.fits')
				sky_est.imcopy_func(sub_out+'[1]','galaxy_clean.fits',1,1,nx1,ny1)
				#print 'ok!'
			#exit()

			if model=='edge+ser':
				h1 = h	
				z01 = z0
				out.profiles_galfit_extract(xc,yc,PA,model,z0=z01,h=h1)
				out.edge_prof(file_gal_clean,xc,yc,m0,pix2sec,rmax*pix2sec*1.5,m_crit,rmax*ellip*pix2sec*1.5,z0,h)

			if model=='edge':
				h1 = h	
				z01 = z0
				m_crit = m0_d + 4.
				out.profiles_galfit_extract(xc,yc,PA,model,z0=z01,h=h1)
				out.edge_prof(file_gal_clean,xc,yc,m0,pix2sec,rmax*pix2sec*2.,m_crit,rmax*ellip*pix2sec*2.,z0,h)

			# STEP 14: Output files
			file_out_pdf = './pics/'+str(number)+'/'+file_out_pdf
			print bcolors.OKBLUE+ 'Creating the output files ...' + bcolors.ENDC
			print >> ff, 'Creating the output files ...'
			status = int(status)
			if n_it==1:
				res.results(galaxy_name,model,image_out,file_out_txt,file_out_pdf,filter_name,colour_name,colour_value,Aext,z,pix2sec,m0,D,status,tables='tog',gal_number=k,out=0,ima_number=numb)
			else:
				res.results(galaxy_name,model,file_galfit_outimage,file_out_txt,file_out_pdf,filter_name,colour_name,colour_value,Aext,z,pix2sec,m0,D,status,tables='tog',gal_number=k,out=0,ima_number=numb)
			#print rmax

			path_new = './pics/'+str(number)+'/'
			def mover(file_from):
				try:
					shutil.move(file_from,path_new+file_from)
				except:
					print 'There is no %s!' % (file_from)

			if setup.delete_contam==0:
				sub_out = './pics/'+str(number)+'/out_'+str(numb)+'.fits'		
				shutil.copy(sub_out,'out.fits')
				from prep_modules import del_contam
				del_contam.fill_cont()
				mover('galaxy_wt_cont.fits')

			mover('ell.png')
			mover('galaxy.png')
			mover('major_prof.png')
			mover('minor_prof.png')
			mover('model.png')
			mover('radius_prof.png')
			mover('residual.png')
			mover('galaxy.jpeg')
			mover('model.jpeg')
			mover('residual.jpeg')
			#if sph=='NO':	mover('galaxy_clean.fits')
			mover('psf.fits')
			mover('models.txt')
			mover('out.log')
			mover('azim_prof.png')
			mover('ini.txt')
			mover('pitch.png')
			mover('pitch.dat')
			file_out_pdf = setup.file_out_pdf


		
			print '\n********************* RESULTS *********************'
			if model=='exp' or model=='exp+ser':
				print 'DISK:'
				print 'm0d=%.2f (mag arcsec-2), h=%.2f (arcsec)' % (m0_d,h)
				print >> ff, 'm0d=%.2f (mag arcsec-2), h=%.2f (arcsec)' % (m0_d,h)
			if model=='edge' or model=='edge+ser':
				print 'EDGE-ON DISK:'
				print 'm0d=%.2f (mag arcsec-2), h=%.2f (arcsec), z0=%.2f (arcsec)' % (m0_d,h,z0)	
				print >> ff, 'm0d=%.2f (mag arcsec-2), h=%.2f (arcsec), z0=%.2f (arcsec)' % (m0_d,h,z0)			
			if model=='edge+ser' or model=='exp+ser' or model=='ser':
				print 'BULGE:'
				print 'meb=%.2f (mag arcsec-2), reb=%.2f (arcsec), n=%.2f' % (me_bul,re_bul,n_bul)
				print >> ff, 'meb=%.2f (mag arcsec-2), reb=%.2f (arcsec), n=%.2f' % (me_bul,re_bul,n_bul)
			print '\n****************************************************'
			if sph=='YES':
				h = h/pix2sec
				z0 = z0/pix2sec
				re_bul = re_bul/pix2sec
				rmax = rmax

		else:
			if mode!=5:
				print bcolors.FAIL+ 'GALFIT could not find the good model fot this object!!!'+ bcolors.ENDC
				ffff = open(file_out_txt,'a')
				print >>ffff, "%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%i\t%i" % (number,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,99999.,model,1,0)
				ffff.close()

		#print 'rmax=%.3f' % (rmax)
		#exit()
		if sph=='YES' or mode==5:
			print '/pics/%i/galaxy_wt_cont.fits' % (number)
			if os.path.exists('pics/%i/galaxy_wt_cont.fits' % (number))==True:
				shutil.copy('pics/%i/galaxy_wt_cont.fits' % (number),'galaxy_clean.fits')
				#print 'ok!'
			#exit()

			status = 1
			# Sophisticated analysis of galaxies
			print bcolors.OKBLUE+ 'Sophisticated analysis of the edge-on galaxy:' + bcolors.ENDC
			
			soph.main(ff,file_gal_clean,model,number,numb,me_bul,re_bul,n_bul,NOISE,m0,pix2sec,FWHM,xc,yc,m0_d,h,z0,b_image*kron_r,rmax,status)
			path_new = './pics/'+str(number)+'/'
			def mover(file_from):
				try:
					shutil.move(file_from,path_new+file_from)
				except:
					print 'There is no %s!' % (file_from)

			mover('n_r_sph.png')
			mover('yc_r_sph.png')
			mover('yc_r.png')
			mover('z0_r.png')
			mover('z0_r_sph.png')
			mover('gal_med.png')
			mover('edge_sph_prof.png')
			mover('galaxy_clean.fits')
			'''
			# 1. What image to work with:
			#image_out = './pics/'+str(number)+'/out_'+str(numb)+'.fits'
			sub_out = './pics/'+str(number)+'/sub_'+str(numb)+'.fits'		
			#shutil.copy(image_out,file_galfit_outimage)
			shutil.copy(sub_out,'subcomps.fits')






			print bcolors.OKBLUE+ 'Sophisticated analysis ...' + bcolors.ENDC
			# STEP 15: Sophisticated Analysis
			if model=='ser' or model=='exp+ser' or model=='edge+ser':
				#1. Bulge shape picture
				median_filter(gal,'galaxy_median.fits')
				plots.image('galaxy_median.fits','galaxy_median',0,float(ny)/float(nx))

			if model=='edge' or model=='edge+ser':
				#2. Edge-on disk
				xc,yc,m0d,h,z0,xc_s,yc_s,m0d_s,h_s,z0_s,n_s,xx,yy = edge_soph.main(ff,file_gal_clean_new,xc,yc,rmin,rmax,zmin1,zmin2,zmax,z0,m0,pix2sec)
				#3. Warps
				omega,omega_l,omega_r,alfa_s,r_last_right,r_last_left,y_last_right,y_last_left,psi,A_right,C_right,A_left,C_left = warp_pars(xx,yy)
			'''	

		if numb_iter==3:
			# DUST attenuation
			s=2
		

		time_finish = time.time() 
		duration = time_finish - time_begin
		print '*********** NUMBER %s is done ***********' % (str(number))
		print bcolors.OKBLUE+ 'THE END' + bcolors.ENDC
		print 'time=%.3fsec' % (duration)
		print '****************************************'
		print '========================================\n\n\n'
		print >> ff, 'time=%.3fsec' % (duration)	# The duration of this DECA running 
		ff.close()
		cleaner = 'sh ' + DECA_PATH + '/cleaner.sh'
		subprocess.call(cleaner, shell=True)	# To remove all old files
		subprocess.call('killall -9 vocl.e', shell=True)	# To remove all old files
		subprocess.call('killall -9 galfit', shell=True)	# To remove all old files
		subprocess.call('killall -9 python', shell=True)	# To remove all old files

		#subprocess.call("sh cleaner.sh", shell=True)	# To remove all old files

	main()








