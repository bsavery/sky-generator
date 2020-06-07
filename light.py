# Libraries
from math import ceil, cos, exp, pi, radians, sin
import numpy as np
import constants as con
import properties as prop
import functions as fun


# Definitions
# position of camera
cam_pos = np.array([0, 0, con.Re+prop.height], dtype=np.int64)
sun_dir = fun.normalize(radians(prop.sun_lat), 0)
samples = prop.samples
samples_light = prop.samples_light
irradiance = con.irradiance
rayleigh_coeff = con.rayleigh_coeff
mie_coeff = con.mie_coeff
ozone_coeff = con.ozone_coeff


def single_scattering(cam_dir):
    # intersection between camera and top of atmosphere
    B = fun.atmosphere_intersection(cam_pos, cam_dir)
    # distance from camera to top of atmosphere
    AB = fun.distance(cam_pos, B)
    # length of each inscattering step
    segment = AB/samples
    # cosine of angle between camera and sun
    mu = fun.angle_vectors(cam_dir, sun_dir)
    # Rayleigh and Mie contribution
    sumR = 0
    sumM = 0
    optical_depthR = 0
    optical_depthM = 0
    optical_depthO = 0
    phaseR = fun.phase_rayleigh(mu)
    phaseM = fun.phase_mie(mu)

    # for each point along AB
    for i in range(samples):
        # distance between each sample point and A
        step = i*segment+segment/2
        # point for inscattering
        P = cam_pos+step*cam_dir
        # height of P from sea level
        height = fun.distance(con.earth_center, P)-con.Re
        # optical depth
        pr = fun.density_rayleigh(height)*segment
        pm = fun.density_mie(height)*segment
        po = fun.density_ozone(height)*segment
        optical_depthR += pr
        optical_depthM += pm
        optical_depthO += po
        if not fun.surface_intersection(P, sun_dir):
            # intersection between P and top of atmosphere
            C = fun.atmosphere_intersection(P, sun_dir)
            # distance from P to top of atmosphere
            PC = fun.distance(P, C)
            # variables
            densityR = 0
            densityM = 0
            densityO = 0
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
                Q = P+step_light*sun_dir
                # height of Q from sea level
                height_light = fun.distance(con.earth_center, Q)-con.Re
                # optical depth
                densityR += fun.density_rayleigh(height_light)
                densityM += fun.density_mie(height_light)
                densityO += fun.density_ozone(height_light)
            # optical depth
            optical_depth_lightR = densityR*segment_light
            optical_depth_lightM = densityM*segment_light
            optical_depth_lightO = densityO*segment_light
            transmittance = rayleigh_coeff*(optical_depthR+optical_depth_lightR) + 1.11*mie_coeff*(optical_depthM+optical_depth_lightM) + ozone_coeff*(optical_depthO+optical_depth_lightO)
            attenuation = np.exp(-transmittance)
            sumR += attenuation*pr
            sumM += attenuation*pm

    # spectrum at pixel
    return irradiance*(sumR*rayleigh_coeff*phaseR + sumM*mie_coeff*phaseM)