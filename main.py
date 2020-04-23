# Libraries
import math
from multiprocessing import Process, Array
import numpy as np
from PIL import Image
import functions as fun
import variables as var


# Properties
# camera altitude from sea level (in m, max: 60km)
height = 1
# sun rotation (latitude and longitude in degrees)
sun_lat = 3
sun_lon = 0
# divisions of the rays (more divisions make more accurate results)
samples = 32
samples_light = 16
# number of processes (squared number must be near the number of logic processors of the CPU)
nprocess = 3
# image size (in pixels)
pixelsx = 128
pixelsy = 64
# save image to file in local directory
save_img = False
img_name = "sky"


# Definitions
# position of camera
cam_pos = np.array([0, 0, var.Re+height], dtype=np.int64)
# normalize sun rotation
sun_rot = fun.normalize(math.radians(sun_lat), 0)
# image definition
img = Image.new('RGB', (pixelsx, pixelsy), "black")
pixels_shifted = img.load()


def multiprocess():
    processes = []
    halfx = int(pixelsx/2)
    # create shared memory array (that can be accessed by multiple processes at the same time)
    pix = Array('i', halfx*pixelsy*3)
    # split the image in nprocess**2 processes
    for i in range(nprocess):
        for j in range(nprocess):
            p = Process(target=calc_pixel, args=(
                int((halfx/nprocess)*j), int((halfx/nprocess)*(j+1)), int((pixelsy/nprocess)*i), int((pixelsy/nprocess)*(i+1)), pix,
            ))
            processes.append(p)
            p.start()
    # wait until all processes end
    for p in processes:
        p.join()
    # print to final pixels
    pixels = np.zeros([pixelsx, pixelsy, 3], dtype=np.int)
    for i in range(halfx):
        for j in range(pixelsy):
            pixels[i][j] = [pix[i*3*pixelsy+j*3], pix[i*3*pixelsy+j*3+1], pix[i*3*pixelsy+j*3+2]]
            pixels[pixelsx-i-1][j] = [pix[i*3*pixelsy+j*3], pix[i*3*pixelsy+j*3+1], pix[i*3*pixelsy+j*3+2]]
    # shift pixels with sun lon change
    shift = sun_lon/360*pixelsx
    s = 0
    for x in range(pixelsx):
        for y in range(pixelsy):
            if x+shift<pixelsx:
                pixels_shifted[x+shift, y] = tuple(pixels[x][y])
            else:
                if s<shift:
                    pixels_shifted[s, y] = tuple(pixels[x, y])
                    if y==pixelsy-1:
                        s += 1
    # open image
    img.show()
    # save image
    if save_img:
        img.save(img_name+".png","PNG")


def calc_pixel(xmin, xmax, ymin, ymax, pix):
    for i in range(xmin, xmax):
        for j in range(ymin, ymax):
            # camera rotation
            cam_lat = math.radians((1-j/pixelsy)*90)
            cam_lon = math.radians(i/pixelsx*360-180)
            # normalize camera rotation
            cam_rot = fun.normalize(cam_lat, cam_lon)
            # get pixel rgb
            I = fun.intensity(sun_rot, cam_pos, cam_rot, samples, samples_light)
            # convert to srgb
            rgb = fun.spec_to_srgb(I)
            # print to pixels array in shared memory
            for l in range(3):
                pix[i*3*pixelsy+j*3+l] = int(rgb[l]*255)


# multiprocessing
if __name__ == '__main__':
    multiprocess()