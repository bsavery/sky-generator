# Libraries
import functions as fun
import numpy as np
from constants import earth_radius, irradiance, mie_coeff, num_wavelengths, ozone_coeff, rayleigh_coeff, rayleigh_scale, mie_scale
from math import ceil, cos, exp, pi, radians, sin, sqrt, dist
from properties import air_density, altitude, dust_density, ozone_density, steps, steps_light, sun_lat


# Definitions
cam_altitude = 1000 * fun.clamp(altitude, 0.001, 59.999)
cam_pos = np.array([0, 0, earth_radius + cam_altitude])
sun_dir = fun.geographical_to_direction(radians(sun_lat), 0)
coefficients = np.array([rayleigh_coeff, 1.11 * mie_coeff, ozone_coeff], dtype=np.object)
density_multipliers = np.array([air_density, dust_density, ozone_density])


def ray_optical_depth(ray_origin, ray_dir):
    # optical depth along a ray through the atmosphere
    ray_origin = np.array(ray_origin)
    ray_dir = np.array(ray_dir)
    ray_end = fun.atmosphere_intersection(ray_origin, ray_dir)
    ray = ray_origin - ray_end
    ray_length = np.sqrt(ray.dot(ray))
    # step along the ray in segments and accumulate the optical depth along each segment
    segment_length = ray_length / steps_light
    segment = ray_dir * segment_length
    # the density of each segment is evaluated at its middle
    P = ray_origin + 0.5 * segment
    
    # list of all P values in steps
    Ps = np.outer(np.arange(steps_light), segment) + P
    # height is lengh of Ps - earth_radius
    heights = np.linalg.norm(Ps, axis=1) - earth_radius
    # calculate depths for each height
    optical_depths = np.column_stack((np.exp(-heights / rayleigh_scale), # fun.density_rayleigh for each height
                                      np.exp(-heights / mie_scale), # fun.density_mie for each height
                                      np.where(np.logical_or(heights < 10000, heights >= 40000), 0, # fun.density_ozone for each height
                                                np.where(heights < 25000, heights * 1 / 15000 - 2/3, -(heights * 1 / 15000 - 8/3)))))
    
    return np.sum(optical_depths) * segment_length


def single_scattering(ray_dir):
    # intersection between camera and top of atmosphere
    ray_end = fun.atmosphere_intersection(cam_pos, ray_dir)
    # distance from camera to top of atmosphere
    ray_length = dist(cam_pos, ray_end)
    # to compute the inscattering, we step along the ray in segments and
    # accumulate the inscattering as well as the optical depth along each segment
    segment_length = ray_length / steps
    segment = segment_length * ray_dir
    optical_depth = np.zeros(3)
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
        densities = np.array([fun.density_rayleigh(height), fun.density_mie(height), fun.density_ozone(height)])
        density = density_multipliers * densities
        optical_depth += density * segment_length

        # if the Earth isn't in the way, evaluate inscattering from the sun
        if not fun.surface_intersection(P, sun_dir):
            optical_depth_light = ray_optical_depth(P, sun_dir)
            # attenuation of light
            extinction_density = (optical_depth + density_multipliers * optical_depth_light) * coefficients
            attenuation = np.exp(-np.sum(extinction_density))
            scattering_density_R = density[0] * rayleigh_coeff
            scattering_density_M = density[1] * mie_coeff
            # compute spectrum
            spectrum += attenuation * (phase_function_R * scattering_density_R + phase_function_M * scattering_density_M)

        # advance along ray
        P += segment

    # spectrum at pixel
    return spectrum * irradiance * segment_length