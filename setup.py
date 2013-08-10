# DEFAULT INPUT PARAMETERS FOR DECA:
# READ CAREFULLY THE FOLLOWING LINES:

#******************* OUTPUT FILES *******************
file_out_txt = 'deca_out.txt'	
file_out_ini = 'ini_out.txt'	
file_out_pdf = 'plots.pdf'	
file_out_sph = 'sph_out.txt'	
file_out_sex = 'sextr_out.txt'
file_out_sph = 'sph_out.txt'	
file_out_dust = 'dust_out.txt'	
file_backgr_txt = 'backgr.txt'
file_psf_txt = 'psf_estim.txt'
file_inimodels_txt = 'models_ini.txt'
file_out_edge = 'soph_edge.txt'
file_ell_txt = 'ellips.txt'


# INTERRUPTION TIME:
TIME = 1000			# Time of the galfit process in sec. After that time galfit will be interrupted!


#***SExtractor setup***:
# READ MANUAL TO SEXTRACTOR
DETECT_MINAREA = 4
DETECT_THRESH = 3.0
FILTER_NAME = 'default.conv'
DEBLEND_NTHRESH = 5
DEBLEND_MINCONT = 0.00005
CLEAN_PARAM = 2.0

PHOT_AUTOPARAMS = (2.0, 4.0)

BACK_SIZE = 300
BACK_FILTERSIZE = 5
BACKPHOTO_TYPE = 'LOCAL'

MEMORY_PIXSTACK = 3000000 

# For SExtractor, how to select the object from the Sextractor table:
FIND_MAJOR_OBJECT = 0		# 0 - yes, 1 - no
radius_find = 20. 		# in sec, 30, find a target within a box: 2radius_find X 2radius_find
semimajor_find = 20. 		# in sec, 30, find a target with a semimajor axis > semimajor_find

# Filters:
filter_names = ['U','B','V','R','I','J','H','K','u','g','r','i','z']

# Models:
models = ['exp','exp+ser','ser+exp','ser','edge','edge+ser']


# Masking the contaminents:
mask_stars = 1			# 1; 0 - mask only stars in the image, 1 - mask all the objects exepting the galaxy 
star_detect = 0.97		# 0.97; cl_star from Sextractor. star_detect==1 for stars with S/N=inf
delete_contam = 1		# Create an output image without contaminents (they would be placed by pixels from model). 0 - yes, 1 -no.
warps_cut = 1			# Cut out warps during decomposition (for edge-on galaxies with warped discs). 0 - yes, 1 - no.

deca_del_contam = 1		# If you want to let DECA mask contaminents using special algorithm (0), or if not (1)
add_contam_sextr = 0		# + add contaminents from Sextractor (0), or only mask objects found by DECA (1)
lim = 0				# if you want to mask objects with I>lim*rms that are out the faintest isophote of the galaxy, type 0 if you don't want to do that.
gap = 1				# aperture width (=[gap]  in pixels) for the masking contaminents algorithm

# Sky background
sky_find = 5 			# 5; Method to find sky background


# How to cut out the galaxy image
extr_coeff = 2.			# box: extr_coeff*R_Kron*a_image x extr_coeff*R_Kron*b_image 


#PSF: good star selecting
convolution = 'YES'		# 'YES' or "NO'
M_psf = 14.857			# NOT USED
SN_limit = 100			# Minimal Signal-to-noise ratio of the psf stars
star_cl = 0.99			# Selection of stars for psf fitting
gal_cl = 0.05			# cl_star from Sextractor. gal_cl==0. for diffused objects
gal_infl_coef = 2.		# Coefficient of influence of galaxies: kron_radius[k]*max(a_image[k],b_image[k])*gal_infl_coef - galaxy neighbourhood
bright_coeff = 6. 		# mag
kron_coef = 2. 			# To find good isolated PSF stars

#*********** For PSF modelling *************#
window = 'gauss'		# Function to create psf image: 'moffat' or 'gauss'
M_psf = 14.4325 
box_psf = 50			# 50*50
beta = 4.765			# For Moffat function


#******* Initials ********************#
rmin_bulge = 1.			# minimum radius of bulge (in fwhm units: rmin_bulge*fwhm = arcsec). We mask central circle rmin_bulge*fwhm.


#******* Cosmology ******************#
# Model:
omega_M = 0.3		
omega_L = 0.7		
omega_k = 0.0		
h=0.72


# For finding initial parameters. Way 4 
ellip_cr = 0.6			# Non edge-on galaxy: 0<ell<ellip_cr; Probable edge-on or highly inclined:1>ell>=ellip_cr
C1 = 4.2			# Scodeggio et al. (2002) A&A 384,812-825
C2 = 6.0
#Ccr = 4.5			# lg(B/D)~-0.5
Ccr = 3.5			# CI to separate early and late type galaxies
El_cr = 0.8			# If ell<El_cr then it is an edge-on galaxy


# Initial parameters:
way_ini = 2			# Main way to find initial parameters. 
				# 2 - using azimuthal profile;
				# 3 - using the major axis cut profile

# If you want to perform only one decomposition regardless the results:
one_dec = 'NO'			

# Detection of the cutout of the galaxy:
N_min = 10 			# Minimum number of pixels to detect the cutout of the disk 
coeff_disk_backgr = 3.		# coeff_disk_backgr * rms
coeff_backgr = 1.		# Level of disk intensity

# Important!
inner_trunc = 1			# Inner truncation for edge-on galaxies:	0 - YES, 1 - NO
dust_mask = 1			# Mask a thin lane of dust (for edge-on galaxies). 0 - YES, 1 - NO
fit_sky = 1			# 1, if you do not want to let 
fix_ell = 1			# Fix ellipticity of bulge to 0. 0 - YES, 1 - NO
fix_c = 0
lim_pars = 0			# Limit changes of parameters in certain borders. 0 - YES, 1 - NO
