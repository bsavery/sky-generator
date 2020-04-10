# Libraries
import math
import numpy as np
from PIL import Image
import spectrum as spec
import functions as fun
import variables as var


# Geometry
height = 0                          # position of camera from sea level (m) (max: 60km)
cam_pos = np.array([0, 0, var.Re+height], dtype=np.int64)
earth_center = np.array([0, 0, 0])
samples = 2
cam_lat = math.radians(0)
cam_lon = math.radians(0)
sun_lat = math.radians(40)
sun_lon = math.radians(0)


# normalize camera rotation
cam_rot = fun.normalize(cam_lat, cam_lon)
# normalize sun rotation
sun_rot = fun.normalize(sun_lat, sun_lon)
# intersection between camera and top of atmosphere
P = fun.sphere_intersection(cam_pos, cam_rot, earth_center, var.Ra)
# distance from camera to top of atmosphere
AB = fun.distance_points(cam_pos, P)
# length of each step
AP = AB/samples


for i in range(samples):
    # distance between each sample
    d = i*AP+AP/2
    # point for inscattering
    P1 = cam_pos+d*cam_rot

'''
APo = AB/steps
PoH = math.sin(latitude)*APo
AH = math.cos(latitude)*APo
# angle AHO
ang1 = math.atan(cam_pos/AH)
# angle HOA
ang2 = math.pi-math.pi/2-ang1
# angle OPoCo
alpha = ang2+math.pi/2+sun_rot[0]
print("alpha ", alpha)
OPo = math.sqrt(PoH**2+AH**2)
print("OPo ", OPo)
# angle OCoPo
beta = math.asin(cam_pos*math.sin(alpha)/OPo)
# angle PoOCo
gamma = math.pi-alpha-beta
# distance from point P to top of atmosphere C
PC = OPo*math.sin(gamma)/math.sin(beta)
print(PC)
'''

'''
# Image
pixelsx = 256
pixelsy = 64

img = Image.new('RGB', (pixelsx, pixelsy), "black")
pixels = img.load()
'''
#for i in range(img.size[0]):
#    for j in range(img.size[1]):
        # calculate Rayleigh scattering
        #rayleigh = fun.rayleigh_outscattering(lam, I, height, length)
        # covert spectrum to RGB
        #rgb = spec.cs_srgb.spec_to_rgb(rayleigh)
        # print to the final pixels
        #pixels[i,j] = rgb

#img.show()