# Libraries
from math import ceil, cos, exp, pi, radians, sin, sqrt
import numpy as np
import constants as con
import properties as prop
import functions as fun


# Definitions
density_scale_R = prop.air_density
density_scale_M = prop.dust_density
density_scale_O = prop.ozone_density
# position of camera
cam_pos = np.array([0, 0, con.Re+prop.altitude])
sun_dir = fun.geographical_to_direction(radians(prop.sun_lat), 0)
steps = prop.steps
steps_light = prop.steps_light
rayleigh_coeff = con.rayleigh_coeff
mie_coeff = con.mie_coeff
ozone_coeff = con.ozone_coeff


def ray_optical_depth(ray_origin, ray_dir):
    # optical depth along a ray through the atmosphere
    ray_end = fun.atmosphere_intersection(ray_origin, ray_dir)
    ray_length = fun.distance(ray_origin, ray_end)
    # step along the ray in segments and accumulate the optical depth along each segment
    segment_length = ray_length / steps_light
    segment = segment_length * ray_dir
    optical_depth = np.zeros(3)
    # the density of each segment is evaluated at its middle
    P = ray_origin + 0.5 * segment
    
    for _ in range(steps_light):
        # height above sea level
        height = sqrt(np.dot(P, P)) - con.Re
        # accumulate optical depth of this segment
        density = np.array([fun.density_rayleigh(height), fun.density_mie(height), fun.density_ozone(height)])
        optical_depth += segment_length * density
        # advance along ray
        P += segment

    return optical_depth


def single_scattering(ray_dir):
    # intersection between camera and top of atmosphere
    ray_end = fun.atmosphere_intersection(cam_pos, ray_dir)
    # distance from camera to top of atmosphere
    ray_length = fun.distance(cam_pos, ray_end)
    # to compute the inscattering, we step along the ray in segments and
    # accumulate the inscattering as well as the optical depth along each segment
    segment_length = ray_length / steps
    segment = segment_length * ray_dir
    optical_depth_R = 0
    optical_depth_M = 0
    optical_depth_O = 0
    spectrum = np.zeros(con.num_wavelengths)
    # cosine of angle between camera and sun
    mu = np.dot(ray_dir, sun_dir)
    # phase function for scattering
    phase_function_R = fun.phase_rayleigh(mu)
    phase_function_M = fun.phase_mie(mu)
    # the density and in-scattering of each segment is evaluated at its middle
    P = cam_pos + 0.5 * segment

    for _ in range(steps):
        # height above sea level
        height = sqrt(np.dot(P, P)) - con.Re
        # evaluate and accumulate optical depth along the ray
        density_R = density_scale_R * fun.density_rayleigh(height)
        density_M = density_scale_M * fun.density_mie(height)
        density_O = density_scale_O * fun.density_ozone(height)
        optical_depth_R += segment_length * density_R
        optical_depth_M += segment_length * density_M
        optical_depth_O += segment_length * density_O
        # if the Earth isn't in the way, evaluate inscattering from the sun
        if not fun.surface_intersection(P, sun_dir):
            optical_depth_light = ray_optical_depth(P, sun_dir)
            # attenuation of light
            extinction_density_R = (optical_depth_R + density_scale_R * optical_depth_light[0]) * rayleigh_coeff
            extinction_density_M = (optical_depth_M + density_scale_M * optical_depth_light[1]) * 1.11 * mie_coeff
            extinction_density_O = (optical_depth_O + density_scale_O * optical_depth_light[2]) * ozone_coeff
            attenuation = np.exp(-(extinction_density_R + extinction_density_M + extinction_density_O))
            scattering_density_R = density_R * rayleigh_coeff
            scattering_density_M = density_M * mie_coeff
            # compute spectrum
            spectrum += con.irradiance * attenuation * (phase_function_R * scattering_density_R + phase_function_M * scattering_density_M) * segment_length

        # advance along ray
        P += segment

    # spectrum at pixel
    return spectrum