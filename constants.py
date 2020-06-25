# Libraries
import numpy as np
from math import pi


# Globals
h = 6.62607004e-34                                  # Planck's constant
c = 299792458                                       # speed of light (m/s)
T = 5778                                            # sun's Temperature (k)
k = 1.38064852e-23                                  # Boltzmann constant
n = 1.0002926                                       # IOR of air
Hr = 8000                                           # Rayleigh scale height (m)
Hm = 1200                                           # Mie scale height (m)
g = 0.76                                            # aerosols anisotropy
N = 2.504e25                                        # number density of air (molecules/m^3)
Do = 2.687e20                                       # Dobson unit (molecules/m^2)
Ozone_max = 300 * Do / 15e3                         # Maximum number density of ozone molecules (m^-3)
distance = 149.6e9                                  # average distance Earth-Sun (m)
Rs = 695500e3                                       # radius of Sun (m)
Re = 6360e3                                         # radius of Earth (m)
Ra = 6420e3                                         # radius of atmosphere (m)
num_wavelengths = 21                                # number of wavelengths
wavelengths_step = (num_wavelengths - 1) * 10**-9   # step between wavelengths (m)
max_luminous_efficacy = 683                         # maximum luminous efficacy

# Illuminants
Illuminant_D65 = np.array([[3.2404542, -1.5371385, -0.4985314],
                           [-0.9692660, 1.8760108, 0.0415560],
                           [0.0556434, -0.2040259, 1.0572252]])
Illuminant_E = np.array([[2.3706743, -0.9000405, -0.4706338],
                         [-0.5138850, 1.4253036, 0.0885814],
                         [0.0052982, -0.0146949, 1.0093968]])

# wavelengths every 20nm
lam = np.arange(380., 781., 20) * 10**-9
# CIE colour matching function from 380 to 780nm in 20nm intervals
cmf = np.array([
[0.001368000000, 0.000039000000, 0.006450001000],
[0.014310000000, 0.000396000000, 0.067850010000],
[0.134380000000, 0.004000000000, 0.645600000000],
[0.348280000000, 0.023000000000, 1.747060000000],
[0.290800000000, 0.060000000000, 1.669200000000],
[0.095640000000, 0.139020000000, 0.812950100000],
[0.004900000000, 0.323000000000, 0.272000000000],
[0.063270000000, 0.710000000000, 0.078249990000],
[0.290400000000, 0.954000000000, 0.020300000000],
[0.594500000000, 0.995000000000, 0.003900000000],
[0.916300000000, 0.870000000000, 0.001650001000],
[1.062200000000, 0.631000000000, 0.000800000000],
[0.854449900000, 0.381000000000, 0.000190000000],
[0.447900000000, 0.175000000000, 0.000020000000],
[0.164900000000, 0.061000000000, 0.000000000000],
[0.046770000000, 0.017000000000, 0.000000000000],
[0.011359160000, 0.004102000000, 0.000000000000],
[0.002899327000, 0.001047000000, 0.000000000000],
[0.000690078600, 0.000249200000, 0.000000000000],
[0.000166150500, 0.000060000000, 0.000000000000],
[0.000041509940, 0.000014990000, 0.000000000000]])
# blackbody radiation
sun = (2 * pi * h * c * c) / (lam**5 * (np.exp((h * c) / (k * T * lam)) - 1)) * 10**-9
# irradiance on top of atmosphere
irradiance = sun * ((Rs * Rs) / (distance * distance))
# Rayleigh scattering coefficient (m^-1)
rayleigh_coeff = ((8 * pi**3) * (n * n - 1)**2) / (3 * N * lam**4)
# Mie scattering coefficient (m^-1)
mie_coeff = 2e-5
# Ozone cross section (in cm^2/molecule) to absorption coefficient (m^-1)
ozone_cross = np.array([
6.04999720618854E-24,
1.0893103182854E-23,
3.67917967079418E-23,
1.36017282525383E-22,
3.73735792971477E-22,
7.51469261186479E-22,
1.18257044868557E-21,
1.79953556347172E-21,
2.88048754046166E-21,
3.88981479760571E-21,
4.57997871538082E-21,
5.09027352924287E-21,
4.0030863998631E-21,
2.95965464815758E-21,
2.09073684368918E-21,
1.36820899679147E-21,
8.64349280941688E-22,
6.15855599572905E-22,
4.18917236558952E-22,
2.76872520775773E-22,
3.13148927506362E-22,
])
ozone_coeff = ozone_cross * 10**-4 * Ozone_max

def read_filmic_look(path):
	nums = []
	with open(path) as filmic_file:
		for line in filmic_file:
			nums.append(float(line))
	return nums

contrast_high = read_filmic_look("looks/high_contrast.txt")