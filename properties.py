# Libraries
import constants as con
from math import radians
import functions as fun
import numpy as np


# Properties
# camera altitude from sea level (in m, max: 60km)
height = 10
# sun rotation (latitude and longitude in degrees)
sun_lat = 60
sun_lon = 0
# divisions of the rays (more divisions make more accurate results)
samples = 32
samples_light = 16
# number of processes (squared number must be near the number of logic processors of the CPU)
nprocess = 3
# image size (in pixels)
pixelsx = 128
pixelsy = int(pixelsx/2)
# render without black bottom
half = True
# save EXR image
linear = False
# save PNG image
save_img = False
exposure = 2
img_name = "sky"


# Definitions
# position of camera
cam_pos = np.array([0, 0, con.Re+height], dtype=np.int64)
# normalize sun rotation
sun_rot = fun.normalize(radians(sun_lat), 0)