# camera height from sea level (in meters, max: 60km)
altitude = 0
# sun latitude and longitude (in degrees)
sun_lat = 60
sun_lon = 0
# divisions of the rays (more divisions make more accurate results)
samples = 32
samples_light = 16
# density values
air_density = 1
dust_density = 1
ozone_density = 1
# number of processes (squared number must be near the number of logic processors of the CPU)
nprocess = 3
# image size (in pixels)
pixelsx = 128
pixelsy = int(pixelsx/2)
# exposure
exposure = 300
# save PNG image
save_img = False
# image name
img_name = "sky"