# camera altitude from sea level (in meters, max: 60km)
height = 10
# sun latitude and longitude (in degrees)
sun_lat = 60
sun_lon = 0
# divisions of the rays (more divisions make more accurate results)
samples = 32
samples_light = 16
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