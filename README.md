# CSC417 Final Project: Position Based Fluids
![](results/2k_waves.gif)

## Introduction

This is the Final Project CSC417: Physics based animation. This repository contains my code for fluid simulation using [Position Based Fluids] (https://mmacklin.com/pbf_sig_preprint.pdf)

### Project Layout

The directory and file layout is the following:
	README.md
	CMakeLists.txt
	main.cpp
	setup.h
	include/
		Fluid.h
		SpatialHashGrid.h
		kernel.h
		viscocity.h
		vorticity.h
	src/
		Fluid.cpp
		SpatialHashGrid.cpp
		kernel.cpp
		viscocity.cpp
		vorticity.cpp
	shared/
		include/
			visualization.h
		src/
			visualization.cpp
	data/
		Experiment
			frames/
				...
			gif/
				...
	results/
	scripts/
		create_gif.py

The `main.cpp` file is the entry point of the program. It contains all the constants for the algorithm, launches a simulation thread, and handles callbacks for interactive simulation.

The `setup.h` file involves building the initial fluid state from a selection of simulation scenes.

The `Fluid` class stores all the physical attributes of the fluid such as particle position, particle velocities, density estimates and constraints. This is where the main algorithm resides.

The `SpatialHashGrid` class is used for neighborhood search. It stores a hashtable that finds all particle ids within a specified radius of a query point in 3d.

The `kernel` file defines the smoothin kernels for density estimates and gradient estimation. 

The `Viscocity` file implements XSPH viscocity.

The `Vorticity` file implements vorticity confinement.

The `visualization` file handles all visualizatoin components including displaying the fluid, plotting density in real time, and writing frames to disk for replay.

The `data` directory stores written frames and gifs created from these frames.

The `results` directory stores some samples of the method, as well as a writeup of the algorithm.

The `create_gif.py` is a simple python script which takes a directory of frames, builds a gif from these frames, and writes the gif back to the same directory.



### Compilation

### Executation

### TODO

* ~~Setup Simulation Framework~~
* ~~Simulate gravity and collision detection with bounding box~~
* ~~Implement Spatial Hashing Grid for efficient Neighbourhood Finding~~
* ~~Barebones PBF Implementation (dam break)~~
* ~~Implement multiple scenes~~
* ~~UI Toggles~~
* ~~Apply External Forces with mouse~~
* ~~Viscocity / Vorticity working~~
* ~~Optimize~~
* ~~Better command line args~~
* ~~Density plot~~
* ~~Save frames~~
* ~~Clean up and Document~~
* Aesthetics
* GPU

### References

