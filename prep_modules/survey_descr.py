#!/usr/bin/python
import sys
import shutil
import pyfits
import numpy as np
import math

def UKIDSS(gal_image):
	hdulist = pyfits.open(gal_image)#, do_not_scale_image_data=True, mode='update')
	prihdr = hdulist[1].header
	nx = prihdr['NAXIS1']
	ny = prihdr['NAXIS2']
	pix2sec = prihdr['PIXLSIZE']
	GAIN = prihdr['GAIN']
	read_out_noise = prihdr['READNOIS']
	sky_level = prihdr['SKYLEVEL']
	sky_noise = prihdr['SKYNOISE']
	fwhm = prihdr['SEEING'] * pix2sec # now in arcsec
	A = prihdr['EXTINCT']
	prihdr0 = hdulist[0].header
	exptime = prihdr0['EXP_TIME']
	NCOMBINE = prihdr0['NEXP']

	m0 = prihdr['MAGZPT']+2.5*math.log10(exptime) #5.*math.log10(pix2sec)
	#del prihdr['CTYPE1']
	#del prihdr['CTYPE2']
	#hdulist.flush()
	print nx,ny,pix2sec,GAIN,read_out_noise,sky_level,fwhm,m0,A,exptime,NCOMBINE
	return nx,ny,pix2sec,GAIN,read_out_noise,sky_level,fwhm,m0,A,exptime,NCOMBINE


def MASS2(gal_image):
	#http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_2.html
	#http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_8a.html

	def gain_ron(date,hemi,f):
		import datetime, calendar, re
		def day_counter(d):
			year1 = 1997; month1 = 1; day1 = 1
			year2 = int(d[0:4]); month2 = int(d[4:6]); day2 = int(d[6:8])
			date1 = datetime.date(year1, month1, day1) # month and day are 1-base
			date2 = datetime.date(year2, month2, day2)
			days_in_first_year = (datetime.date(year1,12,31)-date1).days
			days_in_last_year = (date2 - datetime.date(year2, 1, 1)).days
			if year1 != year2:
			    n_days = days_in_first_year
			    try:
				    for year in range(year1+1, year2):
					n_days = n_days + 365 + (1*calendar.isleap(year))
			    except:
			    	    no=1
			    n_days = n_days + days_in_last_year 
			else:
			    n_days = abs((datetime.date(year1,12,31)-date1).days - (datetime.date(year1,12,31)-date2).days) 

			return n_days

		#print day_counter(str(20030603))
	

		if int(date[0:2])<=99 and int(date[0:2])>13:
			date = '19' + date
		else:
			date = '20' + date
		#print date,f,hemi
		#exit()
		f =str(f)
		hemi = str(hemi)
		'''
		print date,hemi,f
		date = '19970603'; hemi = 'n'; f = 'k'
		print date,hemi,f
		exit()
		'''

		if day_counter(date)>=day_counter('19970603') and day_counter(date)<=day_counter('19971231') and hemi=='n' and f=='k':
			 gron = [8.0,55.]
		if day_counter(date)>=day_counter('19980101') and day_counter(date)<=day_counter('19981023') and hemi=='n' and f=='k':
			 gron = [8.0,57.]
		if day_counter(date)>=day_counter('19981024') and day_counter(date)<=day_counter('19990913') and hemi=='n' and f=='k':
			 gron = [6.5,54.]
		if day_counter(date)>=day_counter('19990914') and day_counter(date)<=day_counter('20030101') and hemi=='n' and f=='k':
			 gron = [6.5,53.]

		if day_counter(date)>=day_counter('19970603') and day_counter(date)<=day_counter('19971231') and hemi=='n' and f=='h':
			 gron = [10.,40.]
		if day_counter(date)>=day_counter('19980101') and day_counter(date)<=day_counter('19981023') and hemi=='n' and f=='h':
			 gron = [10.,40.]
		if day_counter(date)>=day_counter('19981024') and day_counter(date)<=day_counter('19990913') and hemi=='n' and f=='h':
			 gron = [8.,50.]
		if day_counter(date)>=day_counter('19990914') and day_counter(date)<=day_counter('20030101') and hemi=='n' and f=='h':
			 gron = [6.,41.]

		if day_counter(date)>=day_counter('19970603') and day_counter(date)<=day_counter('19971231') and hemi=='n' and f=='j':
			 gron = [10.,38.]
		if day_counter(date)>=day_counter('19980101') and day_counter(date)<=day_counter('19981023') and hemi=='n' and f=='j':
			 gron = [10.,42.]
		if day_counter(date)>=day_counter('19981024') and day_counter(date)<=day_counter('19990913') and hemi=='n' and f=='j':
			 gron = [10.,44.]
		if day_counter(date)>=day_counter('19990914') and day_counter(date)<=day_counter('20030101') and hemi=='n' and f=='j':
			 gron = [6.,39.]

		if day_counter(date)>=day_counter('19970101') and day_counter(date)<=day_counter('19980228') and hemi=='s' and f=='k':
			 gron = [9.9,45.]			# THERE IS NO DATA ABOUT THIS!!!
		if day_counter(date)>=day_counter('19970101') and day_counter(date)<=day_counter('19980228') and hemi=='s' and f=='h':
			 gron = [9.9,45.]			# THERE IS NO DATA ABOUT THIS!!!

		if day_counter(date)>=day_counter('19970101') and day_counter(date)<=day_counter('19980228') and hemi=='s' and f=='j':
			 gron = [9.9,45.]			# THERE IS NO DATA ABOUT THIS!!!
		if day_counter(date)>=day_counter('19980301') and day_counter(date)<=day_counter('19990226') and hemi=='s' and f=='k':
			 gron = [9.9,45.]
		if day_counter(date)>=day_counter('19990227') and day_counter(date)<=day_counter('20030101') and hemi=='s' and f=='k':
			 gron = [10.0,50.]
		if day_counter(date)>=day_counter('19980301') and day_counter(date)<=day_counter('19990226') and hemi=='s' and f=='h':
			 gron = [8.0,45.]
		if day_counter(date)>=day_counter('19990227') and day_counter(date)<=day_counter('20030201') and hemi=='s' and f=='h':
			 gron = [6.3,41.]
		if day_counter(date)>=day_counter('19980301') and day_counter(date)<=day_counter('19990226') and hemi=='s' and f=='j':
			 gron = [8.5,43.]
		if day_counter(date)>=day_counter('19990227') and day_counter(date)<=day_counter('20030201') and hemi=='s' and f=='j':
			 gron = [6.8,45.]

		return gron[0],gron[1]

	hdulist = pyfits.open(gal_image)
	prihdr = hdulist[0].header


	SH = float(prihdr['SEESH'])
	hemi = str(prihdr['SCANDIR'])
	date = str(prihdr['ORDATE'])
	f = str(prihdr['FILTER'])
	sky_level = float(prihdr['SKYVAL'])
	dev = float(prihdr['SKYSIG'])
	m0 = float(prihdr['MAGZP'])
	nx = int(prihdr['NAXIS1'])
	ny = int(prihdr['NAXIS2'])

	pix2secx = float(prihdr['CDELT1'])*3600.
	pix2secy = float(prihdr['CDELT2'])*3600.

	if np.around(math.fabs(pix2secx),decimals=4)==np.around(math.fabs(pix2secy),decimals=4):
		pix2sec = math.fabs(np.around(pix2secx,decimals=4))
	else:
		pix2sec = [pix2secx,pix2secy]

	exptime = 1.3	#1.3 or 7.8 ???
	NCOMBINE = 6

	#print date,hemi,f
	#exit()
	GAIN,read_out_noise = gain_ron(date,hemi,f)
	fwhm = 3.13*SH - 0.46
	#pix2sec = 1.
	A = 'not_found'
	#print nx,ny,pix2sec,GAIN,read_out_noise,sky_level,fwhm,m0,A,exptime,NCOMBINE
	#!!!!!!!!!!! read-out-noise is in e
	#GAIN = GAIN*NCOMBINE

	return nx,ny,pix2sec,GAIN,read_out_noise,sky_level,fwhm,m0,A,exptime,NCOMBINE
	

def WISE(gal_image):
	# http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec1_2.html#atlas
	# http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4c.html
	hdulist = pyfits.open(gal_image)
	prihdr = hdulist[0].header


	nx = int(prihdr['NAXIS1'])
	ny = int(prihdr['NAXIS2'])
	pix2secx =  float(prihdr['PXSCAL1'])
	pix2secy =  float(prihdr['PXSCAL2'])
	band = int(prihdr['BAND'])
	if band==1 or band==2:
		exptime = 7.7
	else:
		exptime = 8.8
	NCOMBINE = 1 

	if np.around(math.fabs(pix2secx),decimals=4)==np.around(math.fabs(pix2secy),decimals=4):
		pix2sec = math.fabs(np.around(pix2secx,decimals=4))
	else:
		pix2sec = [pix2secx,pix2secy]

	GAIN = 10.
	read_out_noise = 10.

	if band==1:
		fwhm = 6.08
	elif band==2:
		fwhm = 6.84
	elif band==3:
		fwhm=7.36
	elif band==4:
		fwhm==11.99
	A = 'not_found'
	m0 = prihdr['MAGZP']
	sky_level = 'not_found'
	return nx,ny,pix2sec,GAIN,read_out_noise,sky_level,fwhm,m0,A,exptime,NCOMBINE


def gain_dark_SDSS(camcol,band,run):
	if band=='u':
		if camcol==1:
			gaidark = [1.62,9.61]
		if camcol==2:
			if run<1100:
				gaidark = [1.595,12.6025]		
			else:
				gaidark = [1.825,12.6025]
		if camcol==3:
			gaidark = [1.59,8.7025]
		if camcol==4:
			gaidark = [1.6,12.6025]
		if camcol==5:
			gaidark = [1.47,9.3025]
		if camcol==6:
			gaidark = [2.17,7.0225]
	if band=='g':
		if camcol==1:
			gaidark = [3.32,15.6025]
		if camcol==2:
			gaidark = [3.855,1.44]		
		if camcol==3:
			gaidark = [3.845,1.3225]
		if camcol==4:
			gaidark = [3.995,1.96]
		if camcol==5:
			gaidark = [4.05,1.1025]
		if camcol==6:
			gaidark = [4.035,1.8225]
	if band=='r':
		if camcol==1:
			gaidark = [4.71,1.8225]
		if camcol==2:
			gaidark = [4.6,1.00]		
		if camcol==3:
			gaidark = [4.72,1.3225]
		if camcol==4:
			gaidark = [4.76,1.3225]
		if camcol==5:
			gaidark = [4.725,0.81]
		if camcol==6:
			gaidark = [4.895,0.9025]
	if band=='i':
		if camcol==1:
			gaidark = [5.165,7.84]
		if camcol==2:
			if run<1500:
				gaidark = [6.565,5.76]
			if run>1500:
				gaidark = [6.565,6.25]			
		if camcol==3:
			gaidark = [4.86,4.6225]
		if camcol==4:
			if run<1500:
				gaidark = [4.885,6.25]
			if run>1500:
				gaidark = [4.885,7.5625]	
		if camcol==5:
			gaidark = [4.64,7.84]
		if camcol==6:
			gaidark = [4.76,5.0625]

	if band=='z':
		if camcol==1:
			gaidark = [4.745,0.81]
		if camcol==2:
			gaidark = [5.155,1.0]		
		if camcol==3:
			gaidark = [4.885,1.0]
		if camcol==4:
			if run<1500:
				gaidark = [4.775,9.61]
			if run>1500:
				gaidark = [4.775,12.6025]	
		if camcol==5:
			if run<1500:
				gaidark = [3.48,1.8225]
			if run>1500:
				gaidark = [3.48,2.1025]	
		if camcol==6:
			gaidark = [4.69,1.21]
	return gaidark[0],gaidark[1]


def SDSS_dr7(gal_image,aa,kk,airmass):
	#http://www.sdss.org/dr7/algorithms/fluxcal.html
	pix2sec = 0.396
	exptime = 53.907456
	import surveys
	m0 = surveys.m0_dr7(exptime,aa,kk,airmass,pix2sec)


	hdulist = pyfits.open(gal_image)
	prihdr = hdulist[0].header
	nx = int(prihdr['NAXIS1'])
	ny = int(prihdr['NAXIS2'])
	run = int(prihdr['RUN'])
	band = str(prihdr['FILTER'])
	camcol = int(prihdr['CAMCOL'])
		
	GAIN,var = gain_dark_SDSS(camcol,band,run)
	read_out_noise = math.sqrt(var)*GAIN
	sky_level = 'not_found'
	NCOMBINE = 1
	A = 'not_found'
	fwhm = 'not_found'

	return nx,ny,pix2sec,GAIN,read_out_noise,sky_level,fwhm,m0,A,exptime,NCOMBINE


def SDSS_dr8(gal_image):
	#http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
	hdulist = pyfits.open(gal_image)
	prihdr = hdulist[0].header
	nx = int(prihdr['NAXIS1'])
	ny = int(prihdr['NAXIS2'])
	run = int(prihdr['RUN'])
	band = str(prihdr['FILTER'])
	camcol = int(prihdr['CAMCOL'])
	kkk = prihdr['NMGY']
	import surveys
	m0 = surveys.m0_dr8_dn(kkk)
		
	GAIN,var = gain_dark_SDSS(camcol,band,run)
	read_out_noise = math.sqrt(var)*GAIN	# should be <~5 e
	sky_level = 0.
	exptime = 53.907456
	NCOMBINE = 1
	A =  'not_found'
	fwhm = 'not_found'
	pix2sec = 0.396

	return nx,ny,pix2sec,GAIN,read_out_noise,sky_level,fwhm,m0,A,exptime,NCOMBINE

def EFIGI(gal_image,aa,kk,airmass,GAIN):
	read_out_noise = 4.
	hdulist = pyfits.open(gal_image)
	prihdr = hdulist[0].header

	exptime = 53.907456
	nx = int(prihdr['NAXIS1'])
	ny = int(prihdr['NAXIS2'])
	pix2sec =  math.fabs(prihdr['CD1_1'])*3600.
	import surveys
	m0 = surveys.m0_dr7(exptime,aa,kk,airmass,0.396)-5.*math.log10(pix2sec)+5.*math.log10(0.396)
	sky_level = 0.
	NCOMBINE = 1
	A =  'not_found'
	fwhm = 'not_found'
	return nx,ny,pix2sec,GAIN,read_out_noise,sky_level,fwhm,m0,A,exptime,NCOMBINE

def Spitzer(gal_image):
	# http://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/mopex/mopexusersguide/11/
	return nx,ny,pix2sec,GAIN,read_out_noise,sky_level,fwhm,m0,A,exptime,NCOMBINE

#MASS2('gal.fits')	
#header_extr('2435_938_40_13_3.fits')	
	
