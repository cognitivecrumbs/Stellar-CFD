WIP Incomplete documentation

Project to recreate stellar like formuations using base compressible navier stokes equations with self gravitional effects.
Navier stokes formulated with fvm method and roe averaging to determine boundary fluxes.
Gravitational potential solved using poisson equation using fft method. Accelerations found using the gradient of the gravitational potential.

Visualization using opengl in c++. Size of points represent density and color represent temperature.


![demo](demo.gif)

After build with cmake can run main_most_updated.cpp for the following 

Compile and Runing the program 
create build directory
```bash
mkdir build && cd build
``` 
Build a `Makefile` using `cmake`.
```bash
cmake ..
```
Compile the program using the `Makefile`.
```bash
make
```
Run the compiled program.
```bash
./OpenGLApp
```


Built off prior work in Nbody simulations
Implemented fvm code for compressible flow

Change initializers in main_most_updated.cpp
    computes gravitational attraction using fast fft poisson equation
        Implementation may be wrong?