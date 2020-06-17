# Libraries
from math import radians
from multiprocessing import Process, Array
import numpy as np
from PIL import Image
import functions as fun
import properties as prop
import light


# Definitions
pixels_x = prop.pixels_x
pixels_y = prop.pixels_y
# image definition
halfx = int(pixels_x / 2)
img = Image.new('RGB', (pixels_x, pixels_y), "black")
pixels_shifted = img.load()
# number of processes (number of logic processors of the CPU)
nprocess = fun.n_threads()


# calculate chunk of image
def calc_pixel(xmin, xmax, pix_array):
    for i in range(xmin, xmax):
        for j in range(pixels_y):
            # camera rotation
            cam_lat = radians((1 - j / pixels_y) * 90)
            cam_lon = radians(i / pixels_x * 360 - 180)
            # normalize camera rotation
            cam_dir = fun.geographical_to_direction(cam_lat, cam_lon)
            # get spectrum from camera direction
            spectrum = light.single_scattering(cam_dir)
            # convert spectrum to xyz
            xyz = fun.spec_to_xyz(spectrum)
            # convert xyz to rgb
            rgb = fun.xyz_to_rgb(xyz, prop.exposure)
            # print to pixels array in shared memory
            pos = i * 3 * pixels_y + j * 3
            pix_array[pos] = int(rgb[0] * 255)
            pix_array[pos + 1] = int(rgb[1] * 255)
            pix_array[pos + 2] = int(rgb[2] * 255)


def multiprocess():
    processes = []
    # create shared memory array that can be accessed by multiple processes at the same time
    pix_array = Array('f', halfx * pixels_y * 3)
    # split the image in nprocess vertical chunks
    for i in range(nprocess):
        process = Process(target = calc_pixel, args = (int((halfx / nprocess) * i), int((halfx / nprocess) * (i + 1)), pix_array))
        processes.append(process)
        process.start()

    # wait until all processes end
    for p in processes:
        p.join()
    
    pixels = np.zeros([pixels_x, pixels_y, 3], dtype = np.int)
    # print to final pixels
    for i in range(halfx):
        for j in range(pixels_y):
            pos = i * 3 * pixels_y + j * 3
            r = pix_array[pos]
            g = pix_array[pos + 1]
            b = pix_array[pos + 2]
            # store pixels
            pixels[i][j] = [r, g, b]
            # mirror pixels
            pixels[pixels_x - i - 1][j] = [r, g, b]

    # shift pixels with sun lon change
    shift = prop.sun_lon / 360 * pixels_x
    s = 0
    for x in range(pixels_x):
        for y in range(pixels_y):
            if (x + shift) < pixels_x:
                pixels_shifted[x + shift, y] = tuple(pixels[x][y])
            else:
                if s < shift:
                    pixels_shifted[s, y] = tuple(pixels[x][y])
                    if y == pixels_y - 1:
                        s += 1

    # show image
    img.show()
    # save image
    if prop.save_image:
        img.save(prop.image_name + ".png", "PNG")


# multiprocessing
if __name__ == '__main__':
    multiprocess()