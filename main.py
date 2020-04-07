# Libraries
import math
import numpy as np
from PIL import Image
import spectrum


# Globals
h = 6.62607004*10**-34  # Planck's constant
c = 299792458           # Speed of light (m/s)
T = 5778                # Sun's Temperature (k)
k = 1.38064852*10**-23  # Boltzmann constant
n = 1.0002926           # IOR Air
r = 695500000           # Radius of Sun (m)
Hr = 8434.5             # Atmospheric scale height (m)
No = 2.547305*10**25    # Molecular density
dep = 0.0279            # Depolarization factor
distance = 149.6*10**9  # Distance Earth-Sun (m)


# Wavelengths every 5 step
lam = np.arange(380., 781., 5)
lam = lam*10**-9

def blackbody(lam):
    # Planck's Law
    B = (2*math.pi*h*c**2)/(lam**5*(np.exp((h*c)/(k*T*lam))-1))*10**-9
    return B

def top_atmosphere(sun):
    # Irradiance on top of atmosphere
    I = sun*(r**2/distance**2)
    return I

def rayleigh(lam, AM, I):
    # Rayleigh Scattering
    Br = 24*math.pi**3*(Hr/No*lam**-4)*((n**2-1)/(n**2+2))**2*((6+3*dep)/(6-7*dep))
    Trs = np.exp(-Br*AM)
    Gbn = I*Trs
    return Gbn


# Image
pixelsx = 256
pixelsy = 64

img = Image.new('RGB', (pixelsx, pixelsy), "black")
pixels = img.load()

for i in range(img.size[0]):
    for j in range(img.size[1]):
        AM = (i/pixelsx)*20+1
        # calculate blackbody emission
        sun = blackbody(lam)
        # calculate irradiance on top of atmosphere
        I = top_atmosphere(sun)
        # calculate Rayleigh scattering
        spec = rayleigh(lam, AM, I)
        # covert spectrum to RGB
        rgb = spectrum.cs_srgb.spec_to_rgb(spec)
        # print to the final pixels
        pixels[i,j] = rgb

img.show()