# Libraries
import math
import numpy as np
from PIL import Image
import spectrum
import functions as fun


# Wavelengths every 5 step
lam = np.arange(380., 781., 5)
lam = lam*10**-9


# Image
pixelsx = 256
pixelsy = 64

img = Image.new('RGB', (pixelsx, pixelsy), "black")
pixels = img.load()

for i in range(img.size[0]):
    for j in range(img.size[1]):
        AM = (i/pixelsx)*20+1
        # calculate blackbody emission
        sun = fun.blackbody(lam)
        # calculate irradiance on top of atmosphere
        I = fun.top_atmosphere(sun)
        # calculate Rayleigh scattering
        spec = fun.rayleigh(lam, AM, I)
        # covert spectrum to RGB
        rgb = spectrum.cs_srgb.spec_to_rgb(spec)
        # print to the final pixels
        pixels[i,j] = rgb

img.show()