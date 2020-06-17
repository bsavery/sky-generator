# Libraries
from math import cos, exp, pi, sin, sqrt, pow
import numpy as np
import constants as con


# Functions
def density_rayleigh(height):
    return exp(-height / con.Hr)

def density_mie(height):
    return exp(-height / con.Hm)

def density_ozone(height):
    if height < 10000 or height >= 40000:
        return 0
    elif height >= 10000 and height < 25000:
        return 1 / 15000 * height - 2 / 3
    elif height >= 25000 and height < 40000:
        return -(1 / 15000 * height - 8 / 3)

def phase_rayleigh(mu):
    return 3 / (16 * pi) * (1 + mu * mu)

def phase_mie(mu):
    return (3 * (1 - con.g * con.g) * (1 + mu * mu)) / (8 * pi * (2 + con.g * con.g) * pow((1 + con.g * con.g - 2 * con.g * mu), 1.5))

def distance(a, b):
    difference = a - b
    return sqrt(np.dot(difference, difference))

def geographical_to_direction(lat, lon):
    return np.array([cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat)])

def atmosphere_intersection(pos, dir):
    b = -2 * np.dot(dir, -pos)
    c = np.sum(pos * pos) - con.Ra * con.Ra
    t = (-b + sqrt(b * b - 4 * c)) / 2
    return pos + dir * t

def surface_intersection(pos, dir):
    if dir[2] >= 0:
        return False
    t = np.dot(dir, -pos) / np.sum(dir * dir)
    D = pos[0] * pos[0] - 2 * -pos[0] * dir[0] * t + pow(dir[0] * t, 2) + pos[1] * pos[1] - 2 * (-pos[1]) * dir[1] * t + pow(dir[1] * t, 2) + pos[2] * pos[2] - 2 * (-pos[2]) * dir[2] * t + pow(dir[2] * t, 2)
    if D <= con.Re*con.Re:
        return True
    else:
        return False

def spec_to_xyz(spec):
    # integrate color matching function
    return (np.sum(spec[:, np.newaxis] * con.cmf, axis=0) * 20 * 10**-9) * 683

def xyz_to_rgb(xyz, exposure):
    # XYZ to sRGB linear
    sRGBlinear = (con.Illuminant_D65 @ xyz)
    # adjust exposure
    sRGBlinear = 1 - np.exp(-exposure * sRGBlinear)
    # sRGB linear to non-linear sRGB gamma corrected
    sRGBlinear *= exposure
    sRGB = [0, 0, 0]
    for i in range(3):
        if sRGBlinear[i] > 0.0031308:
            sRGB[i] = 1.055 * sRGBlinear[i]**(1 / 2.4) - 0.055
        else:
            sRGB[i] = 12.92 * sRGBlinear[i]
    return sRGB