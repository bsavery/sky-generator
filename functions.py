# Libraries
import os
import sys
from math import cos, exp, pi, sin, sqrt, pow
import numpy as np
import constants as con


# Functions
def n_threads():
    if sys.platform == 'win32':
        return (int)(os.environ['NUMBER_OF_PROCESSORS'])
    else:
        return (int)(os.popen('grep -c cores /proc/cpuinfo').read())

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

def spec_to_xyz(spectrum):
    # xyz tristimulus values
    return (np.sum(spectrum[:, np.newaxis] * con.cmf, axis=0)) * con.wavelengths_step * con.max_luminous_efficacy

def xyz_to_rgb(xyz, exposure):
    # XYZ to sRGB linear
    sRGBlinear = np.dot(con.Illuminant_D65, xyz) * 120000
    # apply exposure
    sRGB_exposed = sRGBlinear * pow(2, exposure)
    # avoid infinite values
    for i in range(3):
        if sRGB_exposed[i] <= 0:
            sRGB_exposed[i] = 1e-5
    # apply filmic log encoding
    sRGB_log = (np.log2(sRGB_exposed / 0.18) + 10) / 16.5
    # avoid negative values
    for i in range(3):
        if sRGB_log[i] <= 0:
            sRGB_log[i] = 1e-5
    # apply look contrast
    sRGB = [0.0, 0.0, 0.0]
    for i in range(3):
        index = int(sRGB_log[i] * 4095)
        sRGB[i] = con.contrast_high[index]

    return sRGB