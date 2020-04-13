# Libraries
import math
import numpy as np
import variables as var
import spectrum as spec


# wavelengths every 5nm
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


def get_rgb(sun_rot, cam_pos, cam_rot, earth_center, samples):
    # angle between camera and sun directions
    angle = math.acos((cam_rot[0]*sun_rot[0]+cam_rot[1]*sun_rot[1]+cam_rot[2]*sun_rot[2])/(module(cam_rot)*module(sun_rot)))
    # intersection between camera and top of atmosphere
    B = sphere_intersection(cam_pos, cam_rot, earth_center, var.Ra)
    # distance from camera to top of atmosphere
    AB = distance_points(cam_pos, B)
    # length of each inscattering step
    AP = AB/samples

    pP = 0
    sumT = 0
    for k in range(samples):
        # distance between each sample point and A
        d = k*AP+AP/2
        # point for inscattering
        P = cam_pos+d*cam_rot
        # intersection between P and top of atmosphere
        C = sphere_intersection(P, sun_rot, earth_center, var.Ra)
        # distance from P to top of atmosphere
        PC = distance_points(P, C)
        # length of each outscattering step
        ds = PC/samples
        pQ = 0
        for l in range(samples):
            # distance between each sample point and P
            d = l*ds+ds/2
            # point for outscattering
            Q = P+d*sun_rot
            # height of Q from sea level
            hQ = distance_points(earth_center, Q)-var.Re
            # density ratio for each Q point
            pQ += density_ratio(hQ)
        # optical depth of CP
        Dcp = pQ*ds
        # height of P from sea level
        hP = distance_points(earth_center, P)-var.Re
        # density ratio for each P point
        pP += density_ratio(hP)
        # optical depth of PA
        Dpa = pP*AP
        # Tcp*Tpa
        Trans = np.exp(-rayleigh_coeff_sea*(Dcp+Dpa))
        sumT += Trans*density_ratio(hP)
    # total intensity at pixel
    I = sun*rayleigh_coeff_sea*phase_rayleigh(angle)*sumT*AP
    # convert to RGB
    rgb = spec.cs_srgb.spec_to_rgb(I)
    
    return rgb