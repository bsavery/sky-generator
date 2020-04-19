# Libraries
import math
import numpy as np
import variables as var
import colour


# wavelengths every 5nm
lam = np.arange(380., 781., 5)*10**-9
lamnm = np.arange(380, 781, 5)
# blackbody radiation
sun = (2*math.pi*var.h*var.c**2)/(lam**5*(np.exp((var.h*var.c)/(var.k*var.T*lam))-1))*10**-9
# irradiance on top of atmosphere
I = sun*(var.Rs**2/var.distance**2)
# Rayleigh coefficient (at sea level)
rayleigh_coeff = (8*math.pi**3*(var.n**2-1)**2)/3*(1/var.N)*(1/lam**4)
# define color data
cmfs = colour.STANDARD_OBSERVERS_CMFS['CIE 1931 2 Degree Standard Observer']
illuminant = colour.ILLUMINANTS_SDS['D65']
# The CIE colour matching function for 380 - 780 nm in 5 nm intervals
cmf = np.loadtxt('cie-cmf.txt', usecols=(1,2,3))


def density_ratio(height, scale):
    return math.exp(-height/scale)

def phase_rayleigh(mu):
    return 3/(16*math.pi)*(1+(mu)**2)

def phase_mie(mu):
    return (3/(8*math.pi))*(((1-var.g**2)*(1+mu**2))/((2+var.g**2)*(1+var.g**2-2*var.g*mu)**(3/2)))

def sphere_intersection(pos, rot, center, radius):
    a = rot[0]**2+rot[1]**2+rot[2]**2
    b = -2*(rot[0]*(center[0]-pos[0])+rot[1]*(center[1]-pos[1])+rot[2]*(center[2]-pos[2]))
    c = (center[0]-pos[0])**2+(center[1]-pos[1])**2+(center[2]-pos[2])**2-radius**2
    t = (-b+math.sqrt(b**2-4*a*c))/2*a
    return np.array([pos[0]+rot[0]*t, pos[1]+rot[1]*t, pos[2]+rot[2]*t])

def distance_points(P1, P2):
    return math.sqrt((P2[0]-P1[0])**2+(P2[1]-P1[1])**2+(P2[2]-P1[2])**2)

def normalize(lat, lon):
    return np.array([math.cos(lat)*math.cos(lon), math.cos(lat)*math.sin(lon), math.sin(lat)])

def angle_vectors(rot1, rot2):
    return math.acos(rot1[0]*rot2[0]+rot1[1]*rot2[1]+rot1[2]*rot2[2])


def get_rgb(sun_rot, cam_pos, cam_rot, earth_center, samples):
    # intersection between camera and top of atmosphere
    B = sphere_intersection(cam_pos, cam_rot, earth_center, var.Ra)
    # distance from camera to top of atmosphere
    AB = distance_points(cam_pos, B)
    # length of each inscattering step
    segment = AB/samples
    # cosine of angle between camera and sun directions
    mu = math.cos(angle_vectors(cam_rot, sun_rot))
    # Rayleigh and Mie contribution
    sumR = 0
    sumM = 0
    optical_depthR = 0
    optical_depthM = 0
    phaseR = phase_rayleigh(mu)
    phaseM = phase_mie(mu)

    # for each point along AB
    for i in range(samples):
        # distance between each sample point and A
        distA = i*segment+segment/2
        # point for inscattering
        P = cam_pos+distA*cam_rot
        # height of P from sea level
        height = distance_points(earth_center, P)-var.Re
        # optical depth for light
        hr = math.exp(-height/var.Hr)*segment
        hm = math.exp(-height/var.Hm)*segment
        optical_depthR += hr
        optical_depthM += hm
        # intersection between P and top of atmosphere
        C = sphere_intersection(P, sun_rot, earth_center, var.Ra)
        # distance from P to top of atmosphere
        segment_light = distance_points(P, C)
        # variables
        optical_depth_lightR = 0
        optical_depth_lightM = 0
        # length of each outscattering step
        ds = segment_light/samples

        # for each point along PA
        for j in range(samples):
            # distance between each sample point and P
            distP = j*ds+ds/2
            # point for outscattering
            Q = P+distP*sun_rot
            # height of Q from sea level
            height_light = distance_points(earth_center, Q)-var.Re
            # optical depth
            optical_depth_lightR += math.exp(-height_light/var.Hr)*ds
            optical_depth_lightM += math.exp(-height_light/var.Hm)*ds
 
        tau = rayleigh_coeff*(optical_depthR+optical_depth_lightR)+var.mie_coeff*1.1*(optical_depthM+optical_depth_lightM)
        attenuation = np.exp(-tau)
        sumR += attenuation*hr
        sumM += attenuation*hm

    # total intensity at pixel
    I = sun*(sumR*rayleigh_coeff*phaseR+sumM*var.mie_coeff*phaseM)
    # convert to srgb
    rgb = spec_to_srgb(I)

    return rgb


def spec_to_srgb(spec):
    # spectrum to XYZ
    spec /= 600000
    XYZ = np.sum(spec[:, np.newaxis]*cmf, axis=0)*5
    # XYZ to RGB linear
    D65 = np.array([[3.2404542, -1.5371385, -0.4985314],
                [-0.9692660, 1.8760108, 0.0415560],
                [0.0556434, -0.2040259, 1.0572252]])
    E = np.array([[2.3706743, -0.9000405, -0.4706338],
                [-0.5138850, 1.4253036, 0.0885814],
                [0.0052982, -0.0146949, 1.0093968]])
    RGBlinear = D65 @ XYZ
    # RGB linear to sRGB gamma corrected
    sRGB = [0, 0, 0]
    for i in range(3):
        if RGBlinear[i]>0.0031308:
            sRGB[i] = 1.055*RGBlinear[i]**(1/2.4)-0.055
        else:
            sRGB[i] = 12.92*RGBlinear[i]
    for i in range(3):
        sRGB[i] *= 255

    return sRGB



def spec_to_srgb1(spec):
    data = {}
    for A, B in zip(lamnm, spec):
        data[A] = B/5500
    sd = colour.SpectralDistribution(data)
    # Calculating the sample spectral distribution *CIE XYZ* tristimulus values.
    XYZ = colour.sd_to_XYZ(sd, cmfs, illuminant)
    # The output domain of *colour.sd_to_XYZ* is [0, 100] and the input
    # domain of *colour.XYZ_to_sRGB* is [0, 1]. It needs to be accounted for,
    # thus the input *CIE XYZ* tristimulus values are scaled.
    return colour.XYZ_to_sRGB(XYZ/100)*255