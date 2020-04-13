# Libraries
import math
import numpy as np


def xyz_from_xy(x, y):
    return np.array((x, y, 1-x-y))

class ColourSystem:
    # The CIE colour matching function for 380 - 780 nm in 5 nm intervals
    cmf = np.loadtxt('cie-cmf.txt', usecols=(1,2,3))

    def __init__(self, red, green, blue, white):
        # Chromaticities
        self.red, self.green, self.blue = red, green, blue
        self.white = white
        # The chromaticity matrix (rgb -> xyz) and its inverse
        self.M = np.vstack((self.red, self.green, self.blue)).T 
        self.MI = np.linalg.inv(self.M)
        # White scaling array
        self.wscale = self.MI.dot(self.white)
        # xyz -> rgb transformation matrix
        self.T = self.MI / self.wscale[:, np.newaxis]

    def xyz_to_rgb(self, xyz):
        rgb = self.T.dot(xyz)
        if np.any(rgb < 0):
            w = - np.min(rgb)
            rgb += w
        if not np.all(rgb==0):
            rgb /= np.max(rgb)
        rgb *= 255
        rgb_tuple = tuple(rgb.astype(int))
        return rgb_tuple

    def spec_to_xyz(self, spec):
        XYZ = np.sum(spec[:, np.newaxis] * self.cmf, axis=0)
        den = np.sum(XYZ)
        if den == 0.:
            return XYZ
        return XYZ / den

    def spec_to_rgb(self, spec):
        xyz = self.spec_to_xyz(spec)
        return self.xyz_to_rgb(xyz)

illuminant_D65 = xyz_from_xy(0.3127, 0.3291)
cs_srgb = ColourSystem(red=xyz_from_xy(0.64, 0.33),
                       green=xyz_from_xy(0.30, 0.60),
                       blue=xyz_from_xy(0.15, 0.06),
                       white=illuminant_D65)