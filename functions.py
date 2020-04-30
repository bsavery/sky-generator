# Libraries
from math import ceil, cos, exp, pi, radians, sin, sqrt, pow, radians
import numpy as np
import constants as con


# Functions
def density_rayleigh(height):
    return exp(-height/con.Hr)

def density_mie(height):
    return exp(-height/con.Hm)

def density_ozone(height):
    height /= 1000
    if height<9.52:
        p = 0.1
    elif height>=9.52 and height<37.33:
        p = 0.295*height-3.18
    elif height>=37.33:
        p = 106*exp(-0.0729*height)
    return (p/5)*10**-9

def phase_rayleigh(mu):
    return 3/(16*pi)*(1+mu**2)

def phase_mie(mu):
    return (3*(1-con.g**2)*(1+mu**2))/(8*pi*(2+con.g**2)*(1+con.g**2-2*con.g*mu)**1.5)

def atmosphere_intersection(pos, rot):
    a = pow(rot[0], 2)+pow(rot[1], 2)+pow(rot[2], 2)
    b = -2*(rot[0]*(-pos[0])+rot[1]*(-pos[1])+rot[2]*(-pos[2]))
    c = pow(-pos[0], 2)+pow(-pos[1], 2)+pow(-pos[2], 2)-con.Ra**2
    t = (-b+sqrt(pow(b, 2)-4*a*c))/(2*a)
    return np.array([pos[0]+rot[0]*t, pos[1]+rot[1]*t, pos[2]+rot[2]*t])

def normalize(lat, lon):
    return np.array([cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)])

def surface_intersection(pos, rot):
    if rot[2]>=0:
        return False
    x = -pos[0]
    y = -pos[1]
    z = -pos[2]
    i = rot[0]
    j = rot[1]
    k = rot[2]
    t = (x*i+y*j+z*k)/(pow(i, 2)+pow(j, 2)+pow(k, 2))
    D = pow(x, 2)-2*x*i*t+pow(i*t,2)+pow(y, 2)-2*y*j*t+pow(j*t, 2)+pow(z, 2)-2*z*k*t+pow(k*t, 2)
    if D<=con.Re**2:
        return True
    else:
        return False

def spec_to_srgb(spec, linear, exposure):
    # spectrum to XYZ
    XYZ = (np.sum(spec[:, np.newaxis]*con.cmf, axis=0)*5*10**-9)*683
    # XYZ to sRGB linear
    sRGBlinear = (con.D65 @ XYZ)
    if linear:
        sRGBlinear *= 4
        return sRGBlinear
    else:
        # sRGB linear to non-linear sRGB gamma corrected
        sRGBlinear *= exposure
        sRGB = [0, 0, 0]
        for i in range(3):
            if sRGBlinear[i]>0.0031308:
                sRGB[i] = 1.055*sRGBlinear[i]**(1/2.4)-0.055
            else:
                sRGB[i] = 12.92*sRGBlinear[i]
        return sRGB