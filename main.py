# Libraries
import functions as fun
import light
import os
import sys
from math import radians
from multiprocessing import Process, Array
import numpy as np
from PIL import Image, ImageOps, ImageChops
from properties import exposure, pixels_x, pixels_y, save_image, sun_lon


# image definition
halfx = int(pixels_x / 2)
img = Image.new('RGB', (pixels_x, pixels_y), "black")
pixels = img.load()


# get CPU threads number
def n_threads():
    if sys.platform == 'win32':
        return (int)(os.environ['NUMBER_OF_PROCESSORS'])
    else:
        return (int)(os.popen('grep -c cores /proc/cpuinfo').read())

# number of threads
nprocess = n_threads()


# calculate chunk of image
def calc_pixel(xmin, xmax, pix_array):
    for i in range(xmin, xmax):
        for j in range(pixels_y):
            # camera rotation
            cam_lat = radians((j / pixels_y) * (90 + 0.5))
            cam_lon = radians(i / pixels_x * 360 - 180)
            # normalize camera rotation
            cam_dir = fun.geographical_to_direction(cam_lat, cam_lon)
            # get spectrum from camera direction
            spectrum = light.single_scattering(cam_dir)
            # convert spectrum to xyz
            xyz = fun.spec_to_xyz(spectrum)
            # convert xyz to rgb
            rgb = fun.xyz_to_rgb(xyz, exposure) * 255
            # print to pixels array in shared memory
            rgb_int = rgb.astype(int)
            pos = i * 3 * pixels_y + j * 3
            pix_array[pos] = rgb_int[0]
            pix_array[pos + 1] = rgb_int[1]
            pix_array[pos + 2] = rgb_int[2]


def multiprocess():
    processes = []
    # create shared memory array that can be accessed by multiple processes at the same time
    pix_array = Array('i', halfx * pixels_y * 3)
    # split the image in nprocess vertical chunks
    for i in range(nprocess):
        start_x = int((halfx / nprocess) * i)
        end_x = int((halfx / nprocess) * (i + 1))
        process = Process(target = calc_pixel, args = (start_x, end_x, pix_array))
        processes.append(process)
        process.start()

    # wait until all processes end
    for p in processes:
        p.join()
    
    # print to final pixels
    for i in range(halfx):
        for j in range(pixels_y):
            pos = i * 3 * pixels_y + j * 3
            # convert RGB to int
            rgb = tuple([pix_array[pos], pix_array[pos + 1], pix_array[pos + 2]])
            # store pixels
            pixels[i, j] = rgb
            # mirror pixels
            pixels[pixels_x - i - 1, j] = rgb

    # offset image by sun longitude
    offset = int(sun_lon/360*pixels_x)
    img_offset = ImageChops.offset(img, offset, 0)
    # flip image along x axis
    img_flip = ImageOps.flip(img_offset)
    # show image
    img_flip.show()
    # save image
    if save_image:
        img_flip.save("sky.png", "PNG")


# multiprocessing
if __name__ == '__main__':
    multiprocess()