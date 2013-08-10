################### The configuration file for DECA ###################

#********* About the image *********#
image_name = 's108_aJ_asky_980321n1060021.fits.gz' 	# Full path to the image
weight_name = 'NONE' 				# Full path to the weight image or 'NONE' if you haven't it
psf_name = 'NONE'			   	# Full path to the psf image or 'NONE' if you haven't it

survey = '2MASS' 				# 'NONE' or 'SDSS_dr7', 'SDSS_dr8', '2MASS', 'WISE', 'UKIDSS'
m0 = 28.0 					# Null point for Poghson formula: m0 = magZP + 2.5log(EXPTIME). 'NO' if you don't know it  yet (only for SDSS_dr7 and SDSS_dr8).
xsize = 1.					# Plate scale x: [arcsec/pix]
ysize = 1.					# Plate scale y: [arcsec/pix]
GAIN = 4.6					# [electrons/ADU] - see GALFIT manual
NCOMBINE = 1					# Number of combined images
EXPTIME = 53.907				# The exposure time: [sec]

fwhm = 2.3 					# [arcsec]. If 'NONE' the program will use the psf image (if it exists) or create it.
sky = 'NONE'					# Estimated sky background: [ADU]. If 'NONE' the program will found it by itself.

# For Poghson formula to find m0 (e.g. for SDSS_dr7). If m0 is given the lines below can be neglected.
# See e.g.http://www.sdss.org/dr7/algorithms/fluxcal.html
airmass = 1.1846848				# Airmass	
aa = -23.749599					# Zeropoint	
kk = 0.062895					# Extinction coefficient


#********** ABOUT THE OBJECT ********** #
coords = (205.5347,35.6542) 			# Write (<RA>,<DEC>), decimals are needed
filter_name = 'J'				# Photometric band
galaxy_name = 'NGC5273'				# This name is necessary for the output file	
colour_name = 'J - K'				# Known colour of this object. If it's unknown use 'NO'. Spaces between letters are necessary!
colour_value = 0.835				# The colour value: [mag]. If it's unknown use 'NO'.
z = 0.0042					# Redshift. If it's unknown use 'NO'.
D = 'NO'					# Distance: [Mpc]. If it's unknown or z is already specified then use 'NO'.
Aext = 0.007					# External extinction: [mag]. If it's unknown use 'NO'.

#********** DECOMPOSITION ANALYSIS ************#
find_sky = 1 					# Options: 0 - use this in GALFIT input, 1 - estimate or check the given sky level
find_psf = 0 					# Options: 0 - psf.fits is given, 1 - psf model on your input fwhm,
						# 2 - psf model on finding best psf parameters (slow),
						# 3 - psf model of the best PSF star, 4 - cut of the best PSF star from the image,
						# 5 - psField will be used to extract psf (only for SDSS),10 - NO CONVOLUTION 
image = 'field' 				# 'full' if you don't need to cut out the investigating object and 'field' if you want to cut out the object
model = ('exp+ser') 				# Available: 'exp', 'ser', 'gau', '2exp', 'edge', 'edge2'
find_model = 'incl'				# Let the program find the best model using minimum CHI2 ('YES' or 'NO')
numb_iter = 1 					# 0 - a model with parameters given below, 1 - only one iteration, 2 - two iterations, 3 - dust model (Baes code)
sph = 'NO'					# Sophisticated analysis of the edge-on galaxy ('YES' or 'NO')

#******** Initial parameters or parameters for the modeling (optional) *********#
ini_use = 'NO'					# 'YES' if you want to use your own initial parameters or 'NO' if not. In case of numb_iter = 0 (model building) this line is not in use.



