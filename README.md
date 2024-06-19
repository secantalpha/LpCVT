# Lp CVT  
Lp Centroidal Voronoi Tessellation and its Applications

## Brief
Here is the source code repository of lpcvt, downloaded from https://xueyuhanlang.github.io and some modifications to be compatible with CGAL 5.0 or above

## How to compile

### Pre-requisites
- CMake 3.20 or higher
- CGAL 5.0 or higher
  - Boost 1.85
  - GMP
  - MPFR

### Environment
set environment variable
   ```
    BOOST_LIBRARYDIR=path/to/your/lib/boost_1_85_0/lib64-msvc-14.3

    BOOST_INCLUDEDIR=path/to/your/lib/boost_1_85_0

    CGAL_DIR=path/to/your/lib/cgal

    GMP_DIR=path/to/your/lib/gmp

    MPFR_DIR=path/to/your/lib/mpfr
   ```
### Build
```
cd LpCVT/sources/LpCVT
mkdir build && cd build
cmake ..
cmake --build . --target install --config Release
```

### Run
```
cd path/LpCVT/source/LpCVT/Release
LpCVT mesh_filename pts_filename
```