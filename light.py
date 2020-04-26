# Libraries
from math import ceil, cos, exp, pi, radians, sin, sqrt
import numpy as np
import constants as con
import properties as prop
import functions as fun


# Definitions
sun_rot = prop.sun_rot
samples = prop.samples
samples_light = prop.samples_light
lam = con.lam
lamnm = con.lamnm
cmf = con.cmf
sun = con.sun
I = con.I
rayleigh_coeff = con.rayleigh_coeff
mie_coeff = con.mie_coeff


# Pixel spectrum
def single_scattering(cam_pos, cam_rot):
    if not fun.surface_intersection(cam_pos, cam_rot):
        # intersection between camera and top of atmosphere
        B = fun.atmosphere_intersection(cam_pos, cam_rot)
        # distance from camera to top of atmosphere
        AB = np.linalg.norm(cam_pos-B)
        # length of each inscattering step
        segment = AB/samples
        # cosine of angle between camera and sun
        mu = np.dot(cam_rot, sun_rot)
        # Rayleigh and Mie contribution
        sumR = 0
        #sumM = 0
        #sumO = 0
        optical_depthR = 0
        #optical_depthM = 0
        #optical_depthO = 0
        phaseR = fun.phase_rayleigh(mu)
        #phaseM = phase_mie(mu)

        # for each point along AB
        for i in range(samples):
            # distance between each sample point and A
            step = i*segment+segment/2
            # point for inscattering
            P = cam_pos+step*cam_rot
            # height of P from sea level
            height = np.linalg.norm(con.earth_center-P)-con.Re
            # optical depth
            pr = fun.density_rayleigh(height)*segment
            #pm = density_mie(height)*segment
            #po = density_ozone(height)*segment
            optical_depthR += pr
            #optical_depthM += pm
            #optical_depthO += po
            if not fun.surface_intersection(P, sun_rot):
                # intersection between P and top of atmosphere
                C = fun.atmosphere_intersection(P, sun_rot)
                # distance from P to top of atmosphere
                PC = np.linalg.norm(P-C)
                # variables
                densityR = 0
                #densityM = 0
                #densityO = 0
                # length of each outscattering step
                segment_light = PC/samples
                # adaptive samples
                if segment_light>segment:
                    adaptive_samples = samples_light
                else:
                    adaptive_samples = ceil((segment_light/segment)*samples_light)
                # for each point along PA
                for j in range(adaptive_samples):
                    # distance between each sample point and P
                    step_light = j*segment_light+segment_light/2
                    # point for outscattering
                    Q = P+step_light*sun_rot
                    # height of Q from sea level
                    height_light = np.linalg.norm(con.earth_center-Q)-con.Re
                    # optical depth
                    densityR += fun.density_rayleigh(height_light)
                    #densityM += density_mie(height_light)
                    #densityO += density_ozone(height_light)
                # optical depth
                optical_depth_lightR = densityR*segment_light
                #optical_depth_lightM = densityM*segment_light
                #optical_depth_lightO = densityO*segment_light

                #transmittance = rayleigh_coeff*(optical_depthR+optical_depth_lightR) + mie_coeff*(optical_depthM+optical_depth_lightM) + ozone_coeff*(optical_depthO+optical_depth_lightO)
                #transmittance = rayleigh_coeff*(optical_depthR+optical_depth_lightR) + mie_coeff*(optical_depthM+optical_depth_lightM)
                transmittance = rayleigh_coeff*(optical_depthR+optical_depth_lightR)
                attenuation = np.exp(-transmittance)
                sumR += attenuation*pr
                #sumM += attenuation*pm
                #sumO += attenuation*po

        # total intensity at pixel
        #return sun*(sumR*rayleigh_coeff*phaseR + sumM*mie_coeff*phaseM + sumO*ozone_coeff)
        #return sun*(sumR*rayleigh_coeff*phaseR + sumM*mie_coeff*phaseM)
        return sun*(sumR*rayleigh_coeff*phaseR)
    else:
        return rayleigh_coeff*0


############################################################################
ray1 = fun.normalize(radians(45), radians(0))
ray2 = fun.normalize(radians(45), radians(180))
ray3 = fun.normalize(radians(-45), radians(180))
ray4 = fun.normalize(radians(-45), radians(0))
'''
ray1 = np.array([0.5, 0.5, 0.5])
ray2 = np.array([0.5, 0.5, -0.5])
ray3 = np.array([0.5, -0.5, 0.5])
ray4 = np.array([0.5, -0.5, -0.5])
ray5 = np.array([-0.5, 0.5, 0.5])
ray6 = np.array([-0.5, 0.5, -0.5])
ray7 = np.array([-0.5, -0.5, 0.5])
ray8 = np.array([-0.5, -0.5, -0.5])
'''

def gathered_light(point, dir):
    angle1 = np.dot(dir, ray1)
    angle2 = np.dot(dir, ray2)
    angle3 = np.dot(dir, ray3)
    angle4 = np.dot(dir, ray4)
    I1 = fun.phase_rayleigh(angle1)*single_scattering(point, ray1)
    I2 = fun.phase_rayleigh(angle2)*single_scattering(point, ray2)
    I3 = fun.phase_rayleigh(angle3)*single_scattering(point, ray3)
    I4 = fun.phase_rayleigh(angle4)*single_scattering(point, ray4)
    return (I1+I2+I3+I4)*(pi/2)


def multiple_scattering(cam_pos, cam_rot):
    # intersection between camera and top of atmosphere
    B = fun.atmosphere_intersection(cam_pos, cam_rot)
    # distance from camera to top of atmosphere
    AB = np.linalg.norm(cam_pos-B)
    # length of each inscattering step
    segment = AB/samples
    # Rayleigh and Mie contribution
    sumR = 0
    Tap = 0
    G = 0
    pr = 0
    # for each point along AB
    for i in range(samples):
        # distance between each sample point and A
        step = i*segment+segment/2
        # point for inscattering
        P = cam_pos+step*cam_rot
        # height of P from sea level
        height = np.linalg.norm(con.earth_center-P)-con.Re
        # optical depth
        G += gathered_light(P, cam_rot)
        pr += fun.density_rayleigh(height)
        Tap += np.exp(-rayleigh_coeff*pr)
        sumR += G*Tap*pr
    sumR *= segment
    total = (rayleigh_coeff/(2*pi))*sumR
    return total