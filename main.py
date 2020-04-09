# Libraries
import math
import numpy as np
from PIL import Image
import spectrum as spec
import functions as fun
import variables as var


# Geometry
height = 1                          # position of camera from sea level (m)
cam_pos = var.Re+height             # position of the camera
cam_rot = (0, 0)                    # rotation of the camera (from 0 to 89,9 degrees, from -180 to 180 degrees)

latitude = 89.999

# calculate distance camera - top atmosphere AB
if latitude == 0:
    AB = math.sqrt((var.Ra)**2-(cam_pos)**2)
if latitude > 0:
    # angle OAB
    alpha = math.radians(90+latitude)
    # angle OBA
    beta = math.asin(cam_pos*math.sin(alpha)/var.Ra)
    # angle AOB
    gamma = math.pi-alpha-beta
    AB = cam_pos*math.sin(gamma)/math.sin(beta)

print(AB)

# Image
pixelsx = 256
pixelsy = 64

img = Image.new('RGB', (pixelsx, pixelsy), "black")
pixels = img.load()

#for i in range(img.size[0]):
#    for j in range(img.size[1]):
        # calculate Rayleigh scattering
        #rayleigh = fun.rayleigh_outscattering(lam, I, height, length)
        # covert spectrum to RGB
        #rgb = spec.cs_srgb.spec_to_rgb(rayleigh)
        # print to the final pixels
        #pixels[i,j] = rgb

#img.show()