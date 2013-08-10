#! /usr/bin/env python

import sys
from os import path
import re
from math import log, ceil, degrees, radians, pi, tan, atan, copysign
import math

from numpy import cos, sin, log, empty, zeros, zeros_like, float_, where, arange, median, std
from scipy.weave import converters
from scipy import weave
from scipy.integrate import simps
from scipy.ndimage import gaussian_filter
from scipy.optimize import fsolve

import pyfits

from pylab import *

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


#SCALE = 0.396 # arcsec per pixel

def read(file_galfit, model):
    # All output scales are in pixels, coordinates are in pixels, surface brightnesses are in mag arcsec^-2. PA is as described in GALFIT
    # If bulge is not found (model='exp') then the output parameters for bulge will be equal 99999.
    hdulist = pyfits.open(file_galfit)
    prihdr = hdulist[2].header

    me_bul=[99999.,99999.]
    re_bul=[99999.,99999.]
    n_bul=[99999.,99999.]

    if model=="exp" or model=="exp+ser":
        m0_d = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_MAG']))
        h_d = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_RS']))
        PA_d =  map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_PA']))
        q_d = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_AR']))
        xc_d = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_XC']))
        yc_d = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_YC']))
        if '1_C0' in prihdr:
            c_d = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_C0']))
        else:
            c_d = [0.,0.]


    if model=="exp+ser":
        me_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_MU_E']))
        re_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_RE']))
        n_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_N']))
        PA_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_PA']))
        try:
            q_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_AR']))
        except:
            q_bul = [1.,0.]
        xc_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_XC']))
        yc_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_YC']))
        if '2_C0' in prihdr:
            c_bul = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['2_C0']))
        else:
            c_bul = [0.,0.]
    return xc_d[0], yc_d[0], q_d[0], PA_d[0], m0_d[0], h_d[0], me_bul[0], re_bul[0],n_bul[0]

def bulge_equal_disk(mueb, re, n, mu0d, h, SCALE):
    """Function finds radius r for which the surface brightness of the disk is
    equal to the surface brightness of the bulge."""
    k = 2.5/log(10)
    nu_n = 2*n - 1/3. + 4./(405*n) + 46./(25515*n*n)
    mu0b = mueb - k*nu_n
    # if the bulge is not brighter than the disk even at the center of the
    # galaxy then return just effective radius of the bulge
    if mu0d < mu0b:
        return int(ceil(re))
    # Convert scales from pixels to arcsecs
    re = re * SCALE
    h = h * SCALE
    # function f represents the difference between brightnesses of
    # the buglge and the disk
    f = lambda r: mu0b + k*nu_n*(abs(r)/re)**(1.0/n) - mu0d - k*abs(r)/h
    # find the zero of this function
    try:
        r0 = fsolve(func = f, x0 = 1)[0]
    except RuntimeWarning: # during the solving something goes wrong
        r0 = re        
    return int(ceil(r0/SCALE))

def pitch_fourier(inputImage, xCen, yCen, innerRadius, galSize, axisRatio, posang, step, nArms):
    """Finds thex fourier spectrum of an image"""
    #print innerRadius, galSize
    logfile = open("pitch_log.log", "a")
    logfile.write("%s   %s\n" % (str(innerRadius), str(galSize)))
    u_range = log(arange(innerRadius, galSize, 1))
    theta_range = arange(-pi, pi, 0.025)
    ut_map = empty((len(u_range), len(theta_range)), dtype=float)
    costra = cos(theta_range)
    sintra = sin(theta_range)
    code1 = r"""
    double x_inclined, y_inclined, x_projected, y_projected;
    double theta, expu;
    double cospa = cos(posang);
    double sinpa = sin(posang);
    int x_pixels, y_pixels;
    for (int i=0; i<u_range.length()[0]; i++){
        expu = exp(u_range(i));
        for (int j=0; j<sintra.length()[0]; j++){
            x_inclined = expu * costra(j) * axisRatio;
            y_inclined = expu * sintra(j);
            x_projected = x_inclined * cospa - y_inclined * sinpa;
            y_projected = x_inclined * sinpa + y_inclined * cospa;
            x_pixels = int(x_projected + xCen);
            y_pixels = int(y_projected + yCen);
            ut_map(i,j) = inputImage(x_pixels, y_pixels);
        }
    }
    """
    weave.inline(code1,
                 ['ut_map', 'posang', 'axisRatio', 'inputImage', 'u_range', 'costra', 'sintra', 'xCen', 'yCen'],
                 type_converters=converters.blitz,
                 compiler = 'gcc',
                 extra_compile_args = ['-O3'])
    rrr = simps(simps(ut_map, theta_range), u_range)
    D = abs(rrr)
    p_range = arange(-30, 30, step)
    code2 = r"""
    double x_inclined, y_inclined;
    double u, theta;
    double cospa = cos(posang);
    double sinpa = sin(posang);
    double pp = double(p);
    double mm = double(nArms);
    double D = 0.0;
    double expu, ppu;
    int x_pixels, y_pixels;
    for (int i=0; i<u_range.length()[0]; i++){
        u = u_range(i);
        expu = exp(u);
        ppu = pp * u;
        for (int j=0; j<theta_range.length()[0]; j++){
            theta = theta_range(j);
            x_inclined = expu * costra(j) * axisRatio;
            y_inclined = expu * sintra(j);
            x_pixels = int(x_inclined * cospa - y_inclined * sinpa + xCen);
            y_pixels = int(x_inclined * sinpa + y_inclined * cospa + yCen);
            D = -mm*theta - ppu;
            realPart(i,j) = inputImage(x_pixels, y_pixels) * cos(D);
            imagPart(i,j) = inputImage(x_pixels, y_pixels) * sin(D);
        }
    }
    """
    four_image_abs = []
    four_image = []
    ut_map = empty((len(u_range), len(theta_range)), dtype=complex)
    realPart = zeros_like(ut_map, dtype='float32')
    imagPart = zeros_like(ut_map, dtype='float32')
    for p in p_range:
        weave.inline(code2,
                     ['realPart',
                      'imagPart',
                      'posang',
                      'axisRatio',
                      'inputImage',
                      'p',
                      'nArms',
                      'u_range',
                      'xCen',
                      'yCen',
                      'costra',
                      'sintra',
                      'theta_range'],
                     type_converters=converters.blitz,
                     compiler = 'gcc',
                     extra_compile_args = ['-O2'])
        ut_map = realPart + 1j*imagPart
        rrr = simps(simps(ut_map, theta_range), u_range)
        four_image.append(rrr / D)
        four_image_abs.append([p, abs(rrr) / D])
    return four_image_abs, four_image, D

def get_pitch(fourier, m):
    """finds pitch angle by the fourier spectrum and the number of arms"""
    peaks = []
    for i in xrange(1, len(fourier) - 1):
        if (fourier[i][1] > fourier[i - 1][1]) and (fourier[i][1] > fourier[i + 1][1]):
            peaks.append(fourier[i])
    maximums = sorted(peaks, key=lambda x:x[1])
    maximums.reverse()
    critPitch = 90
    lim = m / abs(tan(radians(critPitch)))
    if len(maximums) == 1:
        if abs(maximums[0][0]) > lim:
            pmax = maximums[0][0]
            return pmax, degrees(atan(m / abs(pmax)))
        else:
            return -1, -1
    else:
        for i in maximums:
            if abs(i[0]) > lim:
                pmax = i[0]
                return  pmax, degrees(atan(m / abs(pmax)))
    return -1, -1



def getSpiralsModel(fitsFileName, model,SCALE):
    print bcolors.OKBLUE+ 'Spiral arms parameters. Fourier analysis.' + bcolors.ENDC
    import del_contam
    del_contam.fill_cont()
    fitsFileName = 'resid_wt_cont.fits'
    #exit()
    (galaxyCenterX, galaxyCenterY, diskAxisRatio, diskPositionAngle, diskCentralBrightness, diskExpScale,
     bulgeEffBrightness, bulgeEffRadius, bulgeSersicIndex) = read('out.fits', model)
    diskR25 = log(2.512**(25 - diskCentralBrightness)) * diskExpScale
    inputImage = pyfits.open(fitsFileName)
    galData = float_(inputImage[0].data)
    galData = galData.T
    galData[where(galData<0)] = 0.0    #FIXME: think how to run out of negtive values of intensity
    innerRadii = []
    pitchAngles = []
    maximums = []
    if model == "exp+ser":
        start = bulge_equal_disk(bulgeEffBrightness,
                                 bulgeEffRadius, 
                                 bulgeSersicIndex, 
                                 diskCentralBrightness, 
                                 diskExpScale, SCALE)
    else:
        start = 0.25 * diskExpScale
    biggestMaximumValue = 0.0
    optimalPitch = None
    optimalInnerRadius = 10
    start=10
    #print start,diskR25
    #exit()
    # Find optimal inner radius by checking value of the main maximum in the Fourier spectra
    for innerRadius in xrange(start, int(diskR25*0.33)):
        four, four_cmplx, normFactor = pitch_fourier(galData,
                                                     galaxyCenterX,
                                                     galaxyCenterY,
                                                     innerRadius,
                                                     diskR25,
                                                     diskAxisRatio,
                                                     radians(diskPositionAngle),
                                                     0.2, # step
                                                     2) # number of arms
        innerRadii.append(innerRadius)
        pitch = get_pitch(four, 2)[1]
        maximumValue = max([i[1] for i in four])
        maximums.append(maximumValue)
        if maximumValue > biggestMaximumValue:
            biggestMaximumValue = maximumValue
            optimalInnerRadius = innerRadius

    # Now compute the Fourier spectra with the optimal inner radius
    # and smaller step
    optimalFour, optimalFour_cmplx, normFactor = pitch_fourier(galData,
                                                               galaxyCenterX,
                                                               galaxyCenterY,
                                                               optimalInnerRadius,
                                                               diskR25,
                                                               diskAxisRatio,
                                                               radians(diskPositionAngle),
                                                               0.01, # step
                                                               2)
    # for given optimal fourier spectra compute pitch angle
    optimalPmax, optimalPitch = get_pitch(optimalFour, 2)

    angularRotation = degrees(log(diskR25/optimalInnerRadius)/tan(radians(optimalPitch)))
    angularRotation = copysign(angularRotation, -optimalPmax)
    Amax = max(optimalFour_cmplx, key=lambda x: abs(x))
    pa = degrees(atan(imag(Amax) / real(Amax)))

    # Compute inclination in degrees [see http://leda.univ-lyon1.fr/leda/param/incl.html for details]
    galaxyIncl = math.degrees(math.asin(math.sqrt((1-10**(-2*math.log10(1/diskAxisRatio)))/(1-10**(-2*0.5)))))

    '''
    # print spirals-related part of galfit config
    print "0) expdisk2                     # Object type"
    print "1)   %7s  %7s   0 0    # position x, y        [pixel]" % (galaxyCenterX, galaxyCenterY)
    print "3)   %7s             1      # total magnitude" % (diskCentralBrightness)
    print "4)   %7s             1      # disk scale-length    [Pixels]" % (diskExpScale)
    print "9)      0.65             1      # Axis ratio (b/a) -- width of spirals"
    print "10)  %7.2f             1      # Spirals position angle" % (pa)
    print "R0)      log                    # PA rotation func. (tanh, sqrt, log, linear, none)"
    print "R1)  %7i             1      # inner (bar radius)" % optimalInnerRadius
    print "R2)  %7i             1      # outer radius" % diskR25
    print "R3)  %7.1f             1      # Cumul. coord.rotation [degrees]" % (angularRotation)
    print "R4)        1             0      # winding scale"
    print "R9)  %7.1f             1      # inclination of entire galaxy" % (galaxyIncl)
    print "R10) %7.1f             1      # position angle (of galaxy)" % (diskPositionAngle-90)
    #print "C0)        0             1      # Fourier mode."
    print "Z)         0                    #  Skip this model in output image?  (yes=1, no=0)"
    '''
    # save obtained parameters to file
    outFile = open(path.join(path.split(fitsFileName)[0], "pitch.dat"), "w")
    outFile.truncate(0)
    # write header
    outFile.write("#pitch[degrees]  azimuth[degrees]   Pmax   A(Pmax)[]   innerRadius[pix]\n")
    outFile.write("%2.1f              %3.1f             %2.2f     %1.4f         %i\n" % (optimalPitch,
                                                                                      pa,
                                                                                      optimalPmax,
                                                                                      abs(Amax),
                                                                                      optimalInnerRadius))
    figure(figsize=(15/2.54, 13/2.54))
    ax1 = subplot(211)
    xlabel("$r_{inner}$", size=10)
    ylabel("$|A(p_{max}, r_{inner})|$", size=10)
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    plot(innerRadii, maximums)
    plot(optimalInnerRadius, biggestMaximumValue, "ro")

    ax2 = subplot(212)
    xlabel("$p$", size=10)
    ylabel("$|A(p)|$", size=10)
    for tick in ax2.xaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    plot([i[0] for i in optimalFour], [i[1] for i in optimalFour])

    outFigurePNG = path.join(path.split(fitsFileName)[0], "pitch.png")
    outFigureEPS = path.join(path.split(fitsFileName)[0], "pitch.eps")
    savefig(outFigurePNG)
    savefig(outFigureEPS)
    return galaxyCenterX, galaxyCenterY,diskCentralBrightness,diskExpScale,pa,optimalInnerRadius,diskR25,angularRotation,galaxyIncl,diskPositionAngle


