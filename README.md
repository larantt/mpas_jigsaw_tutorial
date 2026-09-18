## A Simple JIGSAW-based Workflow for MPAS-Atmosphere

This repository accompanies the "Generating meshes for MPAS-Atmosphere" mini-tutorial, presented at the 2026 MPAS/WRF Users Workshop.

### Prerequisites

Before generating meshes with the workflow provided by this repository, you'll first need to install [JIGSAW]([url](https://github.com/dengwirda/jigsaw)). The JIGSAW `README.md` file provides installation guidance, though the following should generally be sufficient:
```
git clone https://github.com/dengwirda/jigsaw.git
cd jigsaw && mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<WHERE_TO_INSTALL_JIGSAW>
make -j 4 install
```
(Be sure to set the installation directory by replacing `WHERE_TO_INSTALL_JIGSAW` in the `cmake` command above.)

Additionally, you'll also need a Python environment with at least the `numpy` and `netCDF4` packages. If you don't already have a suitable environemnt, setting up a Python virtual environment may be an easiest:
```
python -m venv mpas_mesh
source mpas_mesh/bin/activate
pip install --upgrade pip
pip install numpy netCDF4
```

Finally, in order to compile the `mkgrid` program, which derives MPAS's required mesh geometry and connectivity information from generating points and their triangulation, you will need an MPI implementation (OpenMPI and MPICH are common options) plus the [PnetCDF library]([url](https://parallel-netcdf.github.io/)). After ensuring that the `mpicc` command is in your `PATH`, and the `PNETCDF` environment variable points to the installation path of the PnetCDF library, you can simply run
```
make
```
to build the `mkgrid` program from `mkgrid.c`.

### Generating your first variable-resolution mesh with JIGSAW

Generating a variable-resolution mesh with a 12-km circular refinement region centered at 38 N, 95 W, relaxing to 60-km grid spacing over a distance of 1600 km is accomplished with the following steps.

1. Run the `create_hfun.py` script to generate an `HFUN.msh` file
2. Run `jigsaw`, specifying `MESH.jig` as its command-line argument, to produce a `MESH.msh` file
3. Run `convert_jigsaw.py` with `MESH.msh` as the arguments to produce `SaveVertices` and `SaveTriangles` files from the `MESH.msh` file
4. Run `create_density.py` to produce a `SaveDensity` file
5. Copy the `hfun.py` file to `SaveCode`
6. Run `srun ./mkgrid`, specifying the minimum mesh resolution in meters (e.g. `3000.0` for 3 km inner region) as its command-line argument, to produce `grid.nc` and `graph.info` files
7. Run `mesh_quality.py` and follow the checks below

#### Notes from MMM tutorial video: https://www.youtube.com/watch?v=Qh3bx5sw_rY
* Note that all meshes have one based indexing to match FORTRAN
* JIGSAW tends to create octagons in the mesh which is not the case with the NCAR produced meshes.
    * *NOTE:* If I am likely to be working on the custom meshes, would be interesting to compare the results/runtimes on these custom meshes to successful simulations on the NCAR generated meshes to see if there is any way I can help out?
* When running `mesh_quality.py`, ensure the following:
    * MPAS cannot have obtuse triangles
    * the Maximum nominal cell size gradient is only a few percent
    * The sum of the cell areas etc is close to 4 $\pi$

#### Compilation notes
* When compiling with the intel OneAPI compilers, the optimization routine for a large grid like used for MPAS-A runs into floating point errors that lead to the JIGSAW optimizer "getting stuck" due to floating point errors. 
* Compiling with the `-fp-model=precise` flag fixed this. [This flag](https://www.intel.com/content/www/us/en/docs/dpcpp-cpp-compiler/developer-guide-reference/2023-1/fp-model-fp.html) "disables optimizations that are not value-safe on floating point data". The [-fPIC compiler flag](https://notes.rdu.im/programming/cpp/compiler/fpic_flag/) is just to tell the compiler to generate code which does not depend on being located at a specific memory address and can be used as a shared module. 

Full compiler command for Tiger3:
```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_C_FLAGS="fp-model=precise -fPIC" \
         -DCMAKE_CXX_FLAGS="-fp-model=precise -fPIC" \
         -DCMAKE_INSTALL_PREFIX="/home/GVILLARI/software/mpas/deps/jigsaw/1.0.0"
make -j 4 install
```