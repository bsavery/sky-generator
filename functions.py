# Libraries
from math import cos, exp, pi, sin, sqrt, pow
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

def distance(P1, P2):
    return sqrt((P2[0]-P1[0])**2+(P2[1]-P1[1])**2+(P2[2]-P1[2])**2)

def normalize(lat, lon):
    return np.array([cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)])

def angle_vectors(rot1, rot2):
    return rot1[0]*rot2[0]+rot1[1]*rot2[1]+rot1[2]*rot2[2]

def atmosphere_intersection(pos, dir):
    a = pow(dir[0], 2)+pow(dir[1], 2)+pow(dir[2], 2)
    b = -2*(dir[0]*(-pos[0])+dir[1]*(-pos[1])+dir[2]*(-pos[2]))
    c = pow(-pos[0], 2)+pow(-pos[1], 2)+pow(-pos[2], 2)-con.Ra**2
    t = (-b+sqrt(pow(b, 2)-4*a*c))/(2*a)
    return np.array([pos[0]+dir[0]*t, pos[1]+dir[1]*t, pos[2]+dir[2]*t])

def surface_intersection(pos, dir):
    if dir[2]>=0:
        return False
    t = (-pos[0]*dir[0]-pos[1]*dir[1]-pos[2]*dir[2])/(pow(dir[0], 2)+pow(dir[1], 2)+pow(dir[2], 2))
    D = pow(-pos[0], 2)-2*-pos[0]*dir[0]*t+pow(dir[0]*t,2)+pow(-pos[1], 2)-2*(-pos[1])*dir[1]*t+pow(dir[1]*t, 2)+pow(-pos[2], 2)-2*(-pos[2])*dir[2]*t+pow(dir[2]*t, 2)
    if D<=con.Re**2:
        return True
    else:
        return False

def spec_to_xyz(spec):
    # integrate color matching function
    return (np.sum(spec[:, np.newaxis]*con.cmf, axis=0)*5*10**-9)*683

def xyz_to_rgb(xyz, linear, exposure):
    # XYZ to sRGB linear
    sRGBlinear = (con.D65 @ xyz)
    if linear:
        return sRGBlinear
    else:
        # adjust exposure
        sRGBlinear = 1-np.exp(-exposure*sRGBlinear)
        # sRGB linear to non-linear sRGB gamma corrected
        sRGBlinear *= exposure
        sRGB = [0, 0, 0]
        for i in range(3):
            if sRGBlinear[i]>0.0031308:
                sRGB[i] = 1.055*sRGBlinear[i]**(1/2.4)-0.055
            else:
                sRGB[i] = 12.92*sRGBlinear[i]
        return sRGB