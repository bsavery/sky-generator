# Libraries
import numpy as np


# Properties
# camera altitude from sea level (in m, max: 60km)
height = 10
# sun rotation (latitude and longitude in degrees)
sun_lat = 2
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
exposure = 300
img_name = "sky"