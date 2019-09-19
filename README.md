# Synopsis
**rebl-AMR** 
![alt text](https://github.com/GEM3D/rebl-AMR-V2/blob/master/f16_bunny.png)

This is an Open-Source Software for Binarized-Octree Mesh Generation around Immersed Geometries. <br/>
The main goal is to generate octrees for adaptive mesh refinement (AMR) of Cartesian domains with immersed complex geometries.

## Features
  * MPI parallelization
  * User-friendly interface   
  * Geometry encoding

## Configuration 
rebl-AMR provides config.sh which can be used to automatically configure the code for Stampede2 (https://www.tacc.utexas.edu/systems/stampede2), Comet (https://www.sdsc.edu/services/hpc/hpc_systems.html) and Bridges (https://www.psc.edu/resources/computing/bridges) clusters 

## Linux 
```
source config.sh 
```
## Installation
rebl-AMR requires the following libraries
  * CMAKE
  * MPI (MPI 3.0 standard compliant version)
  * Zoltan
  * ParMetis
  * HDF5

##  Build  
rebl-AMR uses *CMakeLists.txt* and *CMakeModules* folder to detect the library paths. <br/>
These two components are crucial for complilation of rebl-AMR.
Perform the following steps
```
  cd build
  cmake ..
  make 
```
The executable will be placed in the /bin folder


## Run
```
mpirun -np Nprocessor ./bin/amrGem ./input/Geometry.stl procLevel meshLevel 
```
  * NProcessor: Number of processes (squared number)
  * proclevel : topology adaptation level
  * meshlevel : gradient based adaptation level
 
## Visualization
  * The output is written to the /soln folder 
  * Paraview can be used to visualize the solution
  * Simply open the file ending with xdmf in soln/ 
## Paraview Visualization:
  * Poisson solution with Neumann and Dirichlet Boundary conditions
![alt text](https://github.com/GEM3D/rebl-AMR-V2/blob/master/poisson.png)



## Directory structure
```
rebl-AMR
│   README.md
│   CMakeLists.txt    
│   LICENSE
│   config.sh
│
└─── CMakeModules
│   │   FindMYMPI.cmake: cmake script to find MPI
│   │   FindMYHDF5.cmake: cmake script to find HDF5
|   |   FindParmetis.cmake
|   |   FindZoltan.cmake
|   |   ProcessorCount.cmake
|   |   ResolveCompilerPaths.cmake
|   |   CorrectWindowsPaths.cmake
└─── src
│   │   ReblAmr.cpp     
│   │   templateForest.cpp
│   │   communicate.cpp  
│   │   interpolate.cpp
│   │   solver.cpp
│   │   hdf5xmf.cpp
│   │   geomSTL.cpp
│   │   datatype.cpp
|   |   tree.cpp
│   │   templatePhdf5.cpp
│   │   params.cpp 
│   │   definitions.cpp 
│   │ 
│   └─── include (headers)
│       │   ReblAmr.hpp     
│       │   templateForest.hpp
│       │   communicate.hpp  
│       │   interpolate.hpp
│       │   solver.hpp
│       │   hdf5xmf.hpp
│       │   geomSTL.hpp
│       │   datatype.hpp
|       |   tree.hpp
│       │   templatePhdf5.hpp
│       │   params.hpp 
│       │   definitions.hpp 
│   
└─── bin
│       amrGem (executable)  
│  
│
└─── build   
│       Will be populated by cmake   
│  
│
└─── soln 
│   │   Pxdmf3d1.h5: outputs the file in hdf5 format 
│   │   Pxdmf3d1.xmf: meta data to be used by ParaView
│   
└─── archives
 
```

## Details
For complete documentation visit www.reblamr.com

## Notes 
We welcome any feedbacks by the users and developers <br/>
Please read the LICENSE file for how to use this software

## Acknowledgements
**rebl-AMR** is developed as part of the NSF-GEM3D Award No. 1440638. <br/>
It is developed at the Department of Mechanical Engineering at University Pittsburgh, PA, USA. 


## Contributors
  * Jaber J. Hasbestan (jaber@pitt.edu)
  * Shamsulhaq Basir   (shb105@pitt.edu)
  * Inanc Senocak      (senocak@pitt.edu)


