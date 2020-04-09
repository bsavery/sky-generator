# Libraries
import math
import numpy as np
import variables as var

# Wavelengths every 5 nm
lam = np.arange(380., 781., 5)*10**-9
# blackbody radiation (TO VERIFY!!)
sun = (2*math.pi*var.h*var.c**2)/(lam**5*(np.exp((var.h*var.c)/(var.k*var.T*lam))-1))*10**-9
# irradiance on top of atmosphere
I = sun*(var.Rs**2/var.distance**2)
# Rayleigh coefficient (at sea level)
rayleigh_coeff_sea = (8*math.pi**3*(var.n**2-1)**2)/3*(1/var.N)*(1/lam**4)


def density_ratio(height):
    return math.exp(-height/var.Hr)

def phase_rayleigh(angle):
    # Rayleigh phase function
    return 3/(16*math.pi)*(1+(math.cos(angle))**2)

def rayleigh_outscattering(lam, I, height, length):
    # optical depth
    for i in range(steps):
        optical_depth += density_ratio(height)
    # transmittance
    Tr = np.exp(-B*optical_depth)
    return I*Tr