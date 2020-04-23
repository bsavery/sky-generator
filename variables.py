# Libraries
import numpy as np


# Globals
h = 6.62607004*10**-34  # Planck's constant
c = 299792458           # speed of light (m/s)
T = 5778                # sun's Temperature (k)
k = 1.38064852*10**-23  # Boltzmann constant
n = 1.0002926           # IOR of air
Hr = 8000               # Rayleigh scale height (m)
Hm = 1200               # Mie scale height (m)
g = 0.76                # aerosols anisotropy
N = 2.504*10**25        # number density of air (molecules/m^3)
pn = 0.035              # depolarization factor
distance = 149.6*10**9  # average distance Earth-Sun (m)
Rs = 695500*10**3       # radius of Sun (m)
Re = 6360*10**3         # radius of Earth (m)
Ra = 6420*10**3         # radius of atmosphere (m)
# center of Earth
earth_center = np.array([0, 0, 0])

# Illuminants
D65 = np.array([[3.2404542, -1.5371385, -0.4985314],
                [-0.9692660, 1.8760108, 0.0415560],
                [0.0556434, -0.2040259, 1.0572252]])
E = np.array([[2.3706743, -0.9000405, -0.4706338],
                [-0.5138850, 1.4253036, 0.0885814],
                [0.0052982, -0.0146949, 1.0093968]])