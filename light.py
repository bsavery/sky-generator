# Libraries
from math import ceil, cos, degrees, exp, pi, radians, sin, sqrt
import numpy as np
import constants as con
import properties as prop
import functions as fun
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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
ozone_coeff = con.ozone_coeff


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
        sumM = 0
        sumO = 0
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
            P = cam_pos+step*cam_rot
            # height of P from sea level
            height = np.linalg.norm(con.earth_center-P)-con.Re
            # optical depth
            pr = fun.density_rayleigh(height)*segment
            pm = fun.density_mie(height)*segment
            po = fun.density_ozone(height)*segment
            optical_depthR += pr
            optical_depthM += pm
            optical_depthO += po
            if not fun.surface_intersection(P, sun_rot):
                # intersection between P and top of atmosphere
                C = fun.atmosphere_intersection(P, sun_rot)
                # distance from P to top of atmosphere
                PC = np.linalg.norm(P-C)
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
                    Q = P+step_light*sun_rot
                    # height of Q from sea level
                    height_light = np.linalg.norm(con.earth_center-Q)-con.Re
                    # optical depth
                    densityR += fun.density_rayleigh(height_light)
                    densityM += fun.density_mie(height_light)
                    densityO += fun.density_ozone(height_light)
                # optical depth
                optical_depth_lightR = densityR*segment_light
                optical_depth_lightM = densityM*segment_light
                optical_depth_lightO = densityO*segment_light

                transmittance = rayleigh_coeff*(optical_depthR+optical_depth_lightR) + mie_coeff*(optical_depthM+optical_depth_lightM) + ozone_coeff*(optical_depthO+optical_depth_lightO)
                attenuation = np.exp(-transmittance)
                sumR += attenuation*pr
                sumM += attenuation*pm
                sumO += attenuation*po

        # total intensity at pixel
        return sun*(sumR*rayleigh_coeff*phaseR + sumM*mie_coeff*phaseM - sumO*ozone_coeff)
    else:
        return rayleigh_coeff*0


def single_scattering_noozone(cam_pos, cam_rot):
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
        sumM = 0
        optical_depthR = 0
        optical_depthM = 0
        phaseR = fun.phase_rayleigh(mu)
        phaseM = fun.phase_mie(mu)

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
            pm = fun.density_mie(height)*segment
            optical_depthR += pr
            optical_depthM += pm
            if not fun.surface_intersection(P, sun_rot):
                # intersection between P and top of atmosphere
                C = fun.atmosphere_intersection(P, sun_rot)
                # distance from P to top of atmosphere
                PC = np.linalg.norm(P-C)
                # variables
                densityR = 0
                densityM = 0
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
                    densityM += fun.density_mie(height_light)
                # optical depth
                optical_depth_lightR = densityR*segment_light
                optical_depth_lightM = densityM*segment_light

                transmittance = rayleigh_coeff*(optical_depthR+optical_depth_lightR) + mie_coeff*(optical_depthM+optical_depth_lightM)
                attenuation = np.exp(-transmittance)
                sumR += attenuation*pr
                sumM += attenuation*pm

        # total intensity at pixel
        return sun*(sumR*rayleigh_coeff*phaseR + sumM*mie_coeff*phaseM)
    else:
        return rayleigh_coeff*0


# Pixel spectrum
def single_scattering_rayleigh(cam_pos, cam_rot):
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
        optical_depthR = 0
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
            optical_depthR += pr
            if not fun.surface_intersection(P, sun_rot):
                # intersection between P and top of atmosphere
                C = fun.atmosphere_intersection(P, sun_rot)
                # distance from P to top of atmosphere
                PC = np.linalg.norm(P-C)
                # variables
                densityR = 0
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
                # optical depth
                optical_depth_lightR = densityR*segment_light
                transmittance = np.exp(-rayleigh_coeff*(optical_depthR+optical_depth_lightR))
                sumR += transmittance*pr

        # total intensity at pixel
        return sun*(sumR*rayleigh_coeff*phaseR)
    else:
        return rayleigh_coeff*0


def single_scattering_mie(cam_pos, cam_rot):
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
        #sumR = 0
        sumM = 0
        #sumO = 0
        #optical_depthR = 0
        optical_depthM = 0
        #optical_depthO = 0
        #phaseR = fun.phase_rayleigh(mu)
        phaseM = fun.phase_mie(mu)

        # for each point along AB
        for i in range(samples):
            # distance between each sample point and A
            step = i*segment+segment/2
            # point for inscattering
            P = cam_pos+step*cam_rot
            # height of P from sea level
            height = np.linalg.norm(con.earth_center-P)-con.Re
            # optical depth
            #pr = fun.density_rayleigh(height)*segment
            pm = fun.density_mie(height)*segment
            #po = density_ozone(height)*segment
            #optical_depthR += pr
            optical_depthM += pm
            #optical_depthO += po
            if not fun.surface_intersection(P, sun_rot):
                # intersection between P and top of atmosphere
                C = fun.atmosphere_intersection(P, sun_rot)
                # distance from P to top of atmosphere
                PC = np.linalg.norm(P-C)
                # variables
                #densityR = 0
                densityM = 0
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
                    #densityR += fun.density_rayleigh(height_light)
                    densityM += fun.density_mie(height_light)
                    #densityO += density_ozone(height_light)
                # optical depth
                #optical_depth_lightR = densityR*segment_light
                optical_depth_lightM = densityM*segment_light
                #optical_depth_lightO = densityO*segment_light

                #transmittance = rayleigh_coeff*(optical_depthR+optical_depth_lightR) + mie_coeff*(optical_depthM+optical_depth_lightM) + ozone_coeff*(optical_depthO+optical_depth_lightO)
                #transmittance = rayleigh_coeff*(optical_depthR+optical_depth_lightR) + mie_coeff*(optical_depthM+optical_depth_lightM)
                #transmittance = rayleigh_coeff*(optical_depthR+optical_depth_lightR)
                transmittance = mie_coeff*(optical_depthM+optical_depth_lightM)
                attenuation = np.exp(-transmittance)
                #sumR += attenuation*pr
                sumM += attenuation*pm
                #sumO += attenuation*po

        # total intensity at pixel
        #return sun*(sumR*rayleigh_coeff*phaseR + sumM*mie_coeff*phaseM + sumO*ozone_coeff)
        #return sun*(sumR*rayleigh_coeff*phaseR + sumM*mie_coeff*phaseM)
        #return sun*(sumR*rayleigh_coeff*phaseR)
        return sun*(sumM*mie_coeff*phaseM)
    else:
        return mie_coeff*0


############################################################################

def gathered_light_rayleigh(point, dir):
    # code from the pdf i just shared with you
    N = 16
    lat = np.zeros([N], dtype=np.float)
    rays = np.zeros([N,3], dtype=np.float)
    a = 4*pi/N
    d = sqrt(a)
    Mo = round(pi/d)
    do = pi/Mo
    dy = a/do
    flag = 0
    for m in range(Mo):
        O = pi*(m+0.5)/Mo
        My = round(2*pi*sin(O)/dy)
        for n in range(My):
            y = 2*pi*n/My
            lat[flag] = O
            rays[flag] = np.array([sin(O)*cos(y), sin(O)*sin(y), cos(O)])
            flag += 1

    # integral
    I = 0
    for i in range(N):
        phase = fun.phase_rayleigh(np.dot(dir, rays[i]))
        Is = single_scattering_rayleigh(point, rays[i])
        I += Is*phase

    total = I*(4*pi/N)
    return total

'''
def multiple_scattering(cam_pos, cam_rot):
    # intersection between camera and top of atmosphere
    B = fun.atmosphere_intersection(cam_pos, cam_rot)
    # distance from camera to top of atmosphere
    AB = np.linalg.norm(cam_pos-B)
    # length of each inscattering step
    segment = AB/samples
    # Rayleigh and Mie contribution
    sumR = 0
    Tr = 0
    Gr = 0
    pr = 0
    opticalR = 0
    # for each point along AB
    for i in range(samples):
        # distance between each sample point and A
        step = i*segment+segment/2
        # point for inscattering
        P = cam_pos+step*cam_rot
        # height of P from sea level
        height = np.linalg.norm(con.earth_center-P)-con.Re
        # optical depth
        Gr = gathered_light_rayleigh(P, cam_rot)
        pr = fun.density_rayleigh(height)
        opticalR += pr*segment
        Tr = np.exp(-rayleigh_coeff*opticalR)
        sumR += Gr*pr*Tr
    sumR *= segment
    totalR = (rayleigh_coeff/(4*pi))*sumR
    return totalR
'''

def multiple_scattering(cam_pos, cam_rot):
    # intersection between camera and top of atmosphere
    B = fun.atmosphere_intersection(cam_pos, cam_rot)
    # distance from camera to top of atmosphere
    AB = np.linalg.norm(cam_pos-B)
    # length of each inscattering step
    stepSize = AB/samples
    pr = 0
    totalInscatteringRayleigh = 0
    previousInscatteringRayleigh = 0
    for i in range(samples):
        # distance between each sample point and A
        step = i*stepSize+stepSize/2
        # point for inscattering
        P = cam_pos+step*cam_rot
        # height of P from sea level
        height = np.linalg.norm(con.earth_center-P)-con.Re
        # optical depth
        pr += fun.density_rayleigh(height)
        transmittance = np.exp(-rayleigh_coeff*pr)
        currentInscatteringRayleigh = gathered_light_rayleigh(P, cam_rot)*fun.density_rayleigh(height)*transmittance
        totalInscatteringRayleigh += (currentInscatteringRayleigh+previousInscatteringRayleigh)/2*stepSize
        previousInscatteringRayleigh = currentInscatteringRayleigh
    totalInscatteringRayleigh *= rayleigh_coeff/(4*pi)
    return totalInscatteringRayleigh