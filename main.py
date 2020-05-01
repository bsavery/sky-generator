# Libraries
from math import radians
from multiprocessing import Process, Array
import numpy as np
from PIL import Image
import imageio
import functions as fun
import properties as prop
import light


# Definitions
pixelsx = prop.pixelsx
pixelsy = prop.pixelsy
nprocess = prop.nprocess
linear = prop.linear
# image definition
halfx = int(pixelsx/2)
halfy = int(pixelsy/2)
if prop.half:
    img = Image.new('RGB', (pixelsx, halfy), "black")
else:
    img = Image.new('RGB', (pixelsx, pixelsy), "black")
pixels_shifted = img.load()


# calculate chunk of image
def calc_pixel(xmin, xmax, ymin, ymax, pix):
    for i in range(xmin, xmax):
        for j in range(ymin, ymax):
            # camera rotation
            cam_lat = radians((1-j/halfy)*90)
            cam_lon = radians(i/pixelsx*360-180)
            # normalize camera rotation
            cam_dir = fun.normalize(cam_lat, cam_lon)
            # get spectrum from camera direction
            spectrum = light.single_scattering(cam_dir)
            # convert spectrum to xyz
            xyz = fun.spec_to_xyz(spectrum)
            # convert xyz to rgb
            rgb = fun.xyz_to_rgb(xyz, linear, prop.exposure)
            # print to pixels array in shared memory
            if linear:
                for k in range(3):
                    pix[i*3*halfy+j*3+k] = rgb[k]
            else:
                for k in range(3):
                    pix[i*3*halfy+j*3+k] = int(rgb[k]*255)


def multiprocess():
    processes = []
    # create shared memory array (that can be accessed by multiple processes at the same time)
    pix = Array('f', halfx*halfy*3)
    # split the image in nprocess^2 processes
    for i in range(nprocess):
        for j in range(nprocess):
            p = Process(target=calc_pixel, args=(
                int((halfx/nprocess)*j), int((halfx/nprocess)*(j+1)), int((halfy/nprocess)*i), int((halfy/nprocess)*(i+1)), pix,
            ))
            processes.append(p)
            p.start()
    # wait until all processes end
    for p in processes:
        p.join()
    # print to final pixels
    if linear:
        pixels = np.zeros([pixelsy, pixelsx, 3], dtype=np.float32)
    else:
        pixels = np.zeros([pixelsx, pixelsy, 3], dtype=np.int)
    for i in range(halfx):
        for j in range(halfy):
            if linear:
                pixels[j][i] = [pix[i*3*halfy+j*3], pix[i*3*halfy+j*3+1], pix[i*3*halfy+j*3+2]]
                pixels[j][pixelsx-i-1] = [pix[i*3*halfy+j*3], pix[i*3*halfy+j*3+1], pix[i*3*halfy+j*3+2]]
            else:
                pixels[i][j] = [pix[i*3*halfy+j*3], pix[i*3*halfy+j*3+1], pix[i*3*halfy+j*3+2]]
                pixels[pixelsx-i-1][j] = [pix[i*3*halfy+j*3], pix[i*3*halfy+j*3+1], pix[i*3*halfy+j*3+2]]
    # shift pixels with sun lon change
    if not linear:
        shift = prop.sun_lon/360*pixelsx
        s = 0
        for x in range(pixelsx):
            for y in range(halfy):
                if x+shift<pixelsx:
                    pixels_shifted[x+shift, y] = tuple(pixels[x][y])
                else:
                    if s<shift:
                        pixels_shifted[s, y] = tuple(pixels[x, y])
                        if y==pixelsy-1:
                            s += 1
    # show image
    if not linear:
        img.show()
    # save image
    if linear:
        imageio.imwrite(prop.img_name+'.exr', pixels)
    else:
        if prop.save_img:
            img.save(prop.img_name+".png","PNG")


# multiprocessing
if __name__ == '__main__':
    multiprocess()