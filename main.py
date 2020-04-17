# Libraries
import math
import numpy as np
from PIL import Image
import functions as fun
import variables as var
from multiprocessing import Process, Array


# Geometry
height = 1                                                  # position of camera from sea level (m) (max: 60km)
cam_pos = np.array([0, 0, var.Re+height], dtype=np.int64)   # position of camera
earth_center = np.array([0, 0, 0])                          # center of Earth
# Sun rotation
sun_lat = math.radians(60)
sun_lon = math.radians(0)
# normalize sun rotation
sun_rot = fun.normalize(sun_lat, sun_lon)
# divisions of the rays
samples = 20
# number of processes (final number is squared)
nproc = 3
# image size
pixelsx = 128
pixelsy = 64

# image definition
img = Image.new('RGB', (pixelsx, pixelsy), "black")
pixels = img.load()


def calc_pixel(xmin, xmax, ymin, ymax, pix):
    for i in range(xmin, xmax):
        for j in range(ymin, ymax):
            # camera rotation
            cam_lat = math.radians((1-j/pixelsy)*90)
            cam_lon = math.radians(i/pixelsx*360-180)
            # normalize camera rotation
            cam_rot = fun.normalize(cam_lat, cam_lon)
            rgb = fun.get_rgb(sun_rot, cam_pos, cam_rot, earth_center, samples)
            # print to pixels array
            for l in range(3):
                pix[i*3*pixelsy+j*3+l] = int(rgb[l])


def multithread():
    processes = []
    pix = Array('i', pixelsx*pixelsy*3)
    for i in range(nproc):
        for j in range(nproc):
            p = Process(target=calc_pixel, args=(
                int((pixelsx/nproc)*j), int((pixelsx/nproc)*(j+1)), int((pixelsy/nproc)*i), int((pixelsy/nproc)*(i+1)), pix,
            ))
            processes.append(p)
            p.start()
    
    for p in processes:
        p.join()

    # print to final pixels
    for i in range(pixelsx):
        for j in range(pixelsy):
            pixels[i, j] = (pix[i*3*pixelsy+j*3], pix[i*3*pixelsy+j*3+1], pix[i*3*pixelsy+j*3+2])
    
    # open image
    img.show()


# multiprocessing
if __name__ == '__main__':
    multithread()