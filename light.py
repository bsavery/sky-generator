# Libraries
import functions as fun
import numpy as np
from constants import earth_radius, irradiance, mie_coeff, num_wavelengths, ozone_coeff, rayleigh_coeff
from math import ceil, cos, exp, pi, radians, sin, sqrt, dist
from properties import air_density, altitude, dust_density, ozone_density, steps, steps_light, sun_lat


# Definitions
cam_altitude = 1000.0 * fun.clamp(altitude, 0.001, 59.999)
cam_pos = np.array([0, 0, earth_radius + cam_altitude])
sun_dir = fun.geographical_to_direction(radians(sun_lat), 0)


def ray_optical_depth(ray_origin, ray_dir):
    # optical depth along a ray through the atmosphere
    ray_end = fun.atmosphere_intersection(ray_origin, ray_dir)
    ray_length = dist(ray_origin, ray_end)
    # step along the ray in segments and accumulate the optical depth along each segment
    segment_length = ray_length / steps_light
    segment = segment_length * ray_dir
    optical_depth = np.zeros(3)
    # the density of each segment is evaluated at its middle
    P = ray_origin + 0.5 * segment
    
    for _ in range(steps_light):
        # height above sea level
        height = sqrt(np.dot(P, P)) - earth_radius
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
    ray_length = dist(cam_pos, ray_end)
    # to compute the inscattering, we step along the ray in segments and
    # accumulate the inscattering as well as the optical depth along each segment
    segment_length = ray_length / steps
    segment = segment_length * ray_dir
    optical_depth_R = 0
    optical_depth_M = 0
    optical_depth_O = 0
    spectrum = np.zeros(num_wavelengths)
    # cosine of angle between camera and sun
    mu = np.dot(ray_dir, sun_dir)
    # phase function for scattering
    phase_function_R = fun.phase_rayleigh(mu)
    phase_function_M = fun.phase_mie(mu)
    # the density and in-scattering of each segment is evaluated at its middle
    P = cam_pos + 0.5 * segment

    for _ in range(steps):
        # height above sea level
        height = sqrt(np.dot(P, P)) - earth_radius
        # evaluate and accumulate optical depth along the ray
        density_R = air_density * fun.density_rayleigh(height)
        density_M = dust_density * fun.density_mie(height)
        density_O = ozone_density * fun.density_ozone(height)
        optical_depth_R += segment_length * density_R
        optical_depth_M += segment_length * density_M
        optical_depth_O += segment_length * density_O
        # if the Earth isn't in the way, evaluate inscattering from the sun
        if not fun.surface_intersection(P, sun_dir):
            optical_depth_light = ray_optical_depth(P, sun_dir)
            # attenuation of light
            extinction_density_R = (optical_depth_R + air_density * optical_depth_light[0]) * rayleigh_coeff
            extinction_density_M = (optical_depth_M + dust_density * optical_depth_light[1]) * 1.11 * mie_coeff
            extinction_density_O = (optical_depth_O + ozone_density * optical_depth_light[2]) * ozone_coeff
            attenuation = np.exp(-(extinction_density_R + extinction_density_M + extinction_density_O))
            scattering_density_R = density_R * rayleigh_coeff
            scattering_density_M = density_M * mie_coeff
            # compute spectrum
            spectrum += irradiance * attenuation * (phase_function_R * scattering_density_R + phase_function_M * scattering_density_M) * segment_length

        # advance along ray
        P += segment

    # spectrum at pixel
    return spectrum