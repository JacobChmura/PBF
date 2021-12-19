# CSC417 Final Project: Position Based Fluids
![](results/2k_waves.gif)

## Introduction

[Final project](https://github.com/JacobChmura/PBF)  for *CSC417: Physics based animation*. This repository contains my code for fluid simulation using [Position Based Fluids](https://mmacklin.com/pbf_sig_preprint.pdf)

## Project Layout

The directory and file layout is the following:

```bash
├── CMakeLists.txt
├── data
│   └── debug
│       ├── frames
│       └── gifs
├── include
│   ├── Fluid.h
│   ├── kernel.h
│   ├── SpatialHashGrid.h
│   ├── viscocity.h
│   └── vorticity.h
├── main.cpp
├── README.md
├── results
├── scripts
│   └── create_gif.py
├── setup.h
├── shared
│   ├── include
│   │   └── visualization.h
│   └── src
│       └── visualization.cpp
└── src
    ├── Fluid.cpp
    ├── kernel.cpp
    ├── SpatialHashGrid.cpp
    ├── viscocity.cpp
    └── vorticity.cpp

```


The `main.cpp` file is the entry point of the program. It contains all the constants for the algorithm, launches a simulation thread, and handles callbacks for interactive simulation.

The `setup.h` file involves building the initial fluid state from a selection of simulation scenes.

The `Fluid.cpp` class stores all the physical attributes of the fluid such as particle position, particle velocities, density estimates and constraints. This is where the main algorithm resides.

The `SpatialHashGrid.cpp` class is used for neighborhood search. It stores a hashtable that finds all particle ids within a specified radius of a query point in 3d.

The `kernel.cpp` file defines the smoothing kernels for density estimates and gradient estimation. 

The `Viscocity.cpp` file implements XSPH viscocity.

The `Vorticity.cpp` file implements vorticity confinement.

The `visualization.cpp` file handles all visualizatoin components including displaying the fluid, plotting density in real time, and writing frames to disk for replay.

The `data/` directory stores written frames and gifs created from these frames.

The `results/` directory stores some samples of the method, as well as a writeup of the algorithm.

The `create_gif.py` is a simple python script which takes a directory of frames, builds a gif from these frames, and writes the gif back to the same directory.



## Compilation

Starting in this directory, issue:

```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
```

Followed by:
```bash
make
```



## Executation

The compilation should create an object file *Position-Based-Fluids*.


## Program Interaction

## Saving an Interaction


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

