# Libraries
import math
import numpy as np
import variables as var


# wavelengths every 5nm
lam = np.arange(380., 781., 5)*10**-9
# blackbody radiation (TO VERIFY!!)
sun = (2*math.pi*var.h*var.c**2)/(lam**5*(np.exp((var.h*var.c)/(var.k*var.T*lam))-1))
# irradiance on top of atmosphere
I = sun*(var.Rs**2/var.distance**2)
# Rayleigh coefficient (at sea level)
rayleigh_coeff_sea = (8*math.pi**3*(var.n**2-1)**2)/3*(1/var.N)*(1/lam**4)


def density_ratio(height):
    return math.exp(-height/var.Hr)

def phase_rayleigh(angle):
    return 3/(16*math.pi)*(1+(math.cos(angle))**2)

def sphere_intersection(pos, rot, center, radius):
    a = rot[0]**2+rot[1]**2+rot[2]**2
    b = -2*(rot[0]*(center[0]-pos[0])+rot[1]*(center[1]-pos[1])+rot[2]*(center[2]-pos[2]))
    c = (center[0]-pos[0])**2+(center[1]-pos[1])**2+(center[2]-pos[2])**2-radius**2
    t = (-b+math.sqrt(b**2-4*a*c))/2*a
    return np.array([pos[0]+rot[0]*t, pos[1]+rot[1]*t, pos[2]+rot[2]*t])

def distance_points(P1, P2):
    return math.sqrt((P2[0]-P1[0])**2+(P2[1]-P1[1])**2+(P2[2]-P1[2])**2)

def normalize(lat, lon):
    return np.array([math.cos(lat)*math.cos(lon), math.cos(lat)*math.sin(lon), math.sin(lat)])

def module(P):
    return math.sqrt(P[0]**2+P[1]**2+P[2]**2)