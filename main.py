# Libraries
import math
import numpy as np
from PIL import Image
import spectrum as spec
import functions as fun
import variables as var


# Geometry
height = 0                                                  # position of camera from sea level (m) (max: 60km)
cam_pos = np.array([0, 0, var.Re+height], dtype=np.int64)   # position of camera
earth_center = np.array([0, 0, 0])                          # center of Earth
samples = 10                                                # number of divisions for the rays
# Sun rotation
sun_lat = math.radians(90)
sun_lon = math.radians(0)
# normalize sun rotation
sun_rot = fun.normalize(sun_lat, sun_lon)


# Image
pixelsx = 256
pixelsy = 64

img = Image.new('RGB', (pixelsx, pixelsy), "black")
pixels = img.load()


for i in range(img.size[0]):
    for j in range(img.size[1]):
        # camera rotation
        cam_lat = math.radians(j/pixelsy*90)
        cam_lon = math.radians(i/pixelsx*360-180)
        # normalize camera rotation
        cam_rot = fun.normalize(cam_lat, cam_lon)
        # angle between camera and sun directions
        angle = math.acos((cam_rot[0]*sun_rot[0]+cam_rot[1]*sun_rot[1]+cam_rot[2]*sun_rot[2])/(fun.module(cam_rot)*fun.module(sun_rot)))
        # intersection between camera and top of atmosphere
        B = fun.sphere_intersection(cam_pos, cam_rot, earth_center, var.Ra)
        # distance from camera to top of atmosphere
        AB = fun.distance_points(cam_pos, B)
        # length of each step
        AP = AB/samples

        pP = 0
        sumT = 0
        for k in range(samples):
            # distance between each sample point and A
            d = k*AP+AP/2
            # point for inscattering
            P = cam_pos+d*cam_rot
            # intersection between P and top of atmosphere
            C = fun.sphere_intersection(P, sun_rot, earth_center, var.Ra)
            # distance from P to top of atmosphere
            PC = fun.distance_points(P, C)
            # length of each step
            ds = PC/samples
            pQ = 0
            for l in range(samples):
                # distance between each sample point and P
                d = l*ds+ds/2
                # point for outscattering
                Q = P+d*sun_rot
                # height of Q from sea level
                hQ = fun.distance_points(earth_center, Q)-var.Re
                # density ratio for each Q point
                pQ += fun.density_ratio(hQ)
            # optical depth of CP
            Dcp = pQ*ds
            # height of P from sea level
            hP = fun.distance_points(earth_center, P)-var.Re
            # density ratio for each P point
            pP += fun.density_ratio(hP)
            # optical depth of PA
            Dpa = pP*AP
            # Tcp*Tpa
            Trans = np.exp(-fun.rayleigh_coeff_sea*Dcp*Dpa)
            sumT = Trans*pP
        # total intensity at pixel
        I = fun.sun*fun.rayleigh_coeff_sea*fun.phase_rayleigh(angle)*sumT*AP
        # convert to RGB
        rgb = spec.cs_srgb.spec_to_rgb(I)

        # print to the final pixels
        pixels[i,j] = rgb

img.show()
