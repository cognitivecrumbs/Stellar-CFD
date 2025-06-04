
#ifndef CFD_MESH_H
#define CFD_MESH_H

#include <cmath>
#include <cstdio>
#include <algorithm>


#include <algorithm>
#include <iostream>
#include <thread>
#include <future>
#include <thread_pool.h> // Hypothetical or third-party thread pool header
#include <Eigen/Sparse>
#include <Eigen/Dense>
// #include <fftw3.h> // Ensure you have FFTW installed and linked


class CFDInitializerFull;


class CFDTools {
public:
    // Constructor
    CFDTools() {}

    // Find maximum of three numbers
    float maxs(float a, float b, float c);

    // Find minimum of three numbers
    float mins(float a, float b, float c);

    // Sign(x) function
    int sgn(float v);

    // Minmod Limiter
    float minmod(float a, float b, float c);

    // HLLE Solver
    template <size_t size_t>
    void HLLE(float (&flux)[size_t], float qL[4], float qR[4], float gamma, int normal_x, int normal_y);
    // tii (0 is xx), (1 is yy), (2 is xy); dtdi (0 is x), (1 is y)
    template <size_t size_t>
    void VISC(float (&flux)[size_t], float tii[3], float dtdi[2], float qL[4], float qR[4], int normal_x, int normal_y);

    // Flux solver
    // void FLUXSOLVER(float*** U, float*** dUdt, float*** dUdx, float*** dUdy, int nx, int ny, int ghost_cell, float dx, float dy, float gamma);
    void FLUXSOLVER(float*** U, float** e, float*** dUdt, float*** dUdx, float*** dUdy, float*** tii, float*** dtdi, int nx, int ny, int ghost_cell, float dx, float dy, float gamma);
};


class CFDSimulation {
private:
    CFDTools tools;
    CFDInitializerFull* initializer;
    int nx,ny,ghost_cell;
    float dx,dy;
    float gamma,CFL;
    float** e;
    float*** U;
    float*** U_buffer;
    float*** tii;
    float*** dtdi;
    float*** dUdt;
    float*** dUdx;
    float*** dUdy;
    float dt, dt_buffer;
    void solvePoissonEquation(float***& U, float**& phi, int nx, int ny, int ghost_cell, float dx, float dy, float tol, int max_iters);
    void solvePoissonEquationFFT(float***& U, float**& phi, int nx, int ny, int ghost_cell, float dx, float dy);
    void computeGravitationalAcceleration(float***& U, float**& phi, int nx, int ny, int ghost_cell, float dx, float dy);

public:
    float** phi;
    float* xc;
    float* yc;
    CFDSimulation(CFDInitializerFull* initializer, int nx, int ny, int ghost_cell, float dx, float dy, float gamma);

    // void update(float***& q, int nx, int ny, float dx, float dy, float gamma, float dt);
    void update();
    float*** UPtr();
};


#endif