# Libraries
from math import acos, ceil, cos, exp, pi, sin, sqrt
import numpy as np
import variables as var


# wavelengths every 5nm
lam = np.arange(380., 781., 5)*10**-9
lamnm = np.arange(380, 781, 5)
# The CIE colour matching function for 380 - 780 nm in 5 nm intervals
cmf = np.loadtxt('cie-cmf.txt', usecols=(1,2,3))
# blackbody radiation
sun = (2*pi*var.h*var.c**2)/(lam**5*(np.exp((var.h*var.c)/(var.k*var.T*lam))-1))*10**-9
# irradiance on top of atmosphere
I = sun*(var.Rs**2/var.distance**2)
# Rayleigh scattering coefficient (at sea level)
rayleigh_coeff = ((8*pi**3)*(var.n**2-1)**2)/(3*var.N*lam**4)
# Mie scattering coefficient (m^-1)
mie_coeff = 21*10**-6
# Ozone cross section (in cm^2/molecule) to absorption coefficient
ozone_coeff = np.loadtxt('ozone.txt', usecols=(1))*10**-4*var.N


def density_rayleigh(height):
    return exp(-height/var.Hr)

def density_mie(height):
    return exp(-height/var.Hm)

def density_ozone(height):
    height /= 1000
    if height<9.52:
        p = 0.1
    elif height>=9.52 and height<37.33:
        p = 0.295*height-3.18
    elif height>=37.33:
        p = 106*exp(-0.0729*height)
    return p*10**-9

def phase_rayleigh(mu):
    return (3/(16*pi))*(1+mu**2)

def phase_mie(mu):
    return (3*(1-var.g**2)*(1+mu**2))/(8*pi*(2+var.g**2)*(1+var.g**2-2*var.g*mu)**(3/2))

def sphere_intersection(pos, rot, center, radius):
    a = rot[0]**2+rot[1]**2+rot[2]**2
    b = -2*(rot[0]*(center[0]-pos[0])+rot[1]*(center[1]-pos[1])+rot[2]*(center[2]-pos[2]))
    c = (center[0]-pos[0])**2+(center[1]-pos[1])**2+(center[2]-pos[2])**2-radius**2
    t = (-b+sqrt(b**2-4*a*c))/2*a
    return np.array([pos[0]+rot[0]*t, pos[1]+rot[1]*t, pos[2]+rot[2]*t])

def distance(P1, P2):
    return sqrt((P2[0]-P1[0])**2+(P2[1]-P1[1])**2+(P2[2]-P1[2])**2)

def normalize(lat, lon):
    return np.array([cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)])

def angle_vectors(rot1, rot2):
    return acos(rot1[0]*rot2[0]+rot1[1]*rot2[1]+rot1[2]*rot2[2])


def intensity(sun_rot, cam_pos, cam_rot, samples, samples_light):
    # intersection between camera and top of atmosphere
    B = sphere_intersection(cam_pos, cam_rot, var.earth_center, var.Ra)
    # distance from camera to top of atmosphere
    AB = distance(cam_pos, B)
    # length of each inscattering step
    segment = AB/samples
    # cosine of angle between camera and sun directions
    mu = cos(angle_vectors(cam_rot, sun_rot))
    # Rayleigh and Mie contribution
    sumR = 0
    sumM = 0
    #sumO = 0
    optical_depthR = 0
    optical_depthM = 0
    #optical_depthO = 0
    phaseR = phase_rayleigh(mu)
    phaseM = phase_mie(mu)

    # for each point along AB
    for i in range(samples):
        # distance between each sample point and A
        step = i*segment+segment/2
        # point for inscattering
        P = cam_pos+step*cam_rot
        # height of P from sea level
        height = distance(var.earth_center, P)-var.Re
        # optical depth
        pr = density_rayleigh(height)*segment
        pm = density_mie(height)*segment
        #po = density_ozone(height)*segment
        optical_depthR += pr
        optical_depthM += pm
        #optical_depthO += po
        # intersection between P and top of atmosphere
        C = sphere_intersection(P, sun_rot, var.earth_center, var.Ra)
        # distance from P to top of atmosphere
        PC = distance(P, C)
        # variables
        densityR = 0
        densityM = 0
        #densityO = 0
        # length of each outscattering step
        segment_light = PC/samples
        # adaptive samples
        adaptive_samples = ceil((segment_light/segment)*samples_light)
        if adaptive_samples>samples_light:
            adaptive_samples = samples_light
        # for each point along PA
        for j in range(adaptive_samples):
            # distance between each sample point and P
            step_light = j*segment_light+segment_light/2
            # point for outscattering
            Q = P+step_light*sun_rot
            # height of Q from sea level
            height_light = distance(var.earth_center, Q)-var.Re
            # optical depth
            densityR += density_rayleigh(height_light)
            densityM += density_mie(height_light)
            #densityO += density_ozone(height_light)
        # optical depth
        optical_depth_lightR = densityR*segment_light
        optical_depth_lightM = densityM*segment_light
        #optical_depth_lightO = densityO*segment_light
 
        #transmittance = rayleigh_coeff*(optical_depthR+optical_depth_lightR) + 1.11*mie_coeff*(optical_depthM+optical_depth_lightM) + ozone_coeff*(optical_depthO+optical_depth_lightO)
        transmittance = rayleigh_coeff*(optical_depthR+optical_depth_lightR) + 1.11*mie_coeff*(optical_depthM+optical_depth_lightM)
        attenuation = np.exp(-transmittance)
        sumR += attenuation*pr
        sumM += attenuation*pm
        #sumO += attenuation*po

    # total intensity at pixel
    #return sun*(sumR*rayleigh_coeff*phaseR + sumM*mie_coeff*phaseM + sumO*ozone_coeff)
    return sun*(sumR*rayleigh_coeff*phaseR + sumM*mie_coeff*phaseM)


def spec_to_srgb(spec, linear, exposure):
    # spectrum to XYZ
    XYZ = (np.sum(spec[:, np.newaxis]*cmf, axis=0)*5*10**-9)*683
    # XYZ to sRGB linear
    sRGBlinear = (var.D65 @ XYZ)
    if linear:
        return sRGBlinear
    # sRGB linear to non-linear sRGB gamma corrected
    else:
        sRGBlinear *= exposure
        sRGB = [0, 0, 0]
        for i in range(3):
            if sRGBlinear[i]>0.0031308:
                sRGB[i] = 1.055*sRGBlinear[i]**(1/2.4)-0.055
            else:
                sRGB[i] = 12.92*sRGBlinear[i]
        return sRGB


'''
samples = 10

def transmittanceR(S):
    segment = S/samples
    densityR = 0
    optical_depthR = 0
    # integral from Pa to Pb
    for sample in range(samples):
        step = sample*segment+segment/2
        P = Pa+step*dir
        height = distance(P, var.earth_center)-var.Re
        densityR += density_rayleigh(height)
    optical_depthR = densityR*segment
    return rayleigh_coeff*optical_depthR


def intensity_s(Pa, Pb, cam_rot, sun_rot):
    dist = distance(Pa, Pb)
    segment = dist/samples
    tR = 0
    tM = 0
    tR1 = 0
    tM1 = 0
    densityR = 0
    densityM = 0
    # integral from Pa to Pb
    for sample in range(samples):
        step = sample*segment+segment/2
        P = Pa+step*cam_rot
        height = distance(P, var.earth_center)-var.Re
        densityR += density_rayleigh(height)
        densityM += density_mie(height)
        tR += transmittanceR(P, Pc)
        tR1 += transmittanceR(PaP)
        tM += transmittanceM(PPc)
        tM1 += transmittanceM(PaP)

    integralR = densityR*np.exp(-tR-tR1)*segment
    integralM = densityM*np.exp(-tM-tM1)*segment
    Ir = sun*phase_rayleigh(mu)*(rayleigh_coeff/(4*math.pi))*integralR
    Im = sun*phase_mie(mu)*(var.mie_coeff/(4*math.pi))*integralM
    # total intensity
    Is = Ir+Im
'''