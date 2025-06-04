WIP No documentation yet

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