# Sky Generator
This is a sky generator that calculates a physically based sky based on an improved version of [Nishita](https://www.scratchapixel.com/lessons/procedural-generation-virtual-worlds/simulating-sky/simulating-colors-of-the-sky) 1993 single scattering model.

### Install
* Download this repo
* Install **Python 3.8** or a later version
* Install pip
* Install numpy library
```
pip install numpy
```
* Install pillow library
```
pip install pillow
```

### Run
To run the code just paste `python main.py` in the main folder.

### How it works
When finishing, the program will show the rendered sky image with the OS image visualizer.
To change the sky parameters, just change them in the `properties.py` file, there you will find the sun rotation, camera position altitude along others.
