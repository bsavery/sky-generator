# Libraries
import math
import numpy as np
from PIL import Image
import spectrum as spec
import functions as fun
import variables as var
import threading
from multiprocessing import Pool


# Geometry
height = 0                                                  # position of camera from sea level (m) (max: 60km)
cam_pos = np.array([0, 0, var.Re+height], dtype=np.int64)   # position of camera
earth_center = np.array([0, 0, 0])                          # center of Earth
samples = 30                                                # number of divisions for the rays
# Sun rotation
sun_lat = math.radians(40)
sun_lon = math.radians(0)
# normalize sun rotation
sun_rot = fun.normalize(sun_lat, sun_lon)
# number of threads (the real number will be nthreads**2)
nthreads = 4


# image size
pixelsx = 128
pixelsy = 64

img = Image.new('RGB', (pixelsx, pixelsy), "black")
pixels = img.load()


def calc_pixel(xmin, xmax, ymin, ymax):
    for i in range(xmin, xmax):
        for j in range(ymin, ymax):
            # camera rotation
            cam_lat = math.radians((1-j/pixelsy)*90)
            cam_lon = math.radians(i/pixelsx*360-180)
            # normalize camera rotation
            cam_rot = fun.normalize(cam_lat, cam_lon)
            rgb = fun.get_rgb(sun_rot, cam_pos, cam_rot, earth_center, samples)
            # print to the final pixels
            pixels[i, j] = rgb


def multithread():
    threads = []
    for i in range(nthreads):
        for j in range(nthreads):
            t = threading.Thread(target=calc_pixel, args=(
                int((pixelsx/nthreads)*j), int((pixelsx/nthreads)*(j+1)), int((pixelsy/nthreads)*i), int((pixelsy/nthreads)*(i+1))
            ))
            threads.append(t)
            t.start()
            t.join()
    # open image
    img.show()


# multiprocessing
if __name__ == '__main__':
    multithread()