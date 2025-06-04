#ifndef INITIALIZATION_H
#define INITIALIZATION_H

// #include <particles.h>
#include <vector>
#include <particles.h>
#include <mesh.h>
// #include <functional>

// Define a function pointer type
// typedef void (*FuncPtr)(PhysEngine*);

// takes mesh (which contains particles) and performs physics on it (leapfrog integration)
// includes handlers for boundary conditions

// // forward declare to avoid circular include
// class CFDMesh;

class Initialization
{
public:
    Initialization() {};
    virtual void init(std::vector <Particle*> *particles);
    virtual void reInitParticle(Particle* particle) {};
};

class DiskInitialization: public Initialization
{
public:
    DiskInitialization(float base_mass, float radius_max, float radius_min, float disk_ang_mom, float base_energy, float base_ang_momentum);

    void init(std::vector <Particle*> *particles) override;
    void reInitParticle(Particle* particle) override;
private:
    float base_mass;          // Default mass for particles
    // float disk_ang_mom = sqrt(center_mass+0.5*nparticles*base_mass);       // Angular momentum scaling factor
    float radius_max;          // Maximum radius of the disk
    float radius_min;
    float disk_ang_mom;       // Angular momentum scaling factor
    float base_energy;         // Default energy
    float base_ang_momentum;   // Default angular momentum 
};

class CFDInitializer
{
public:
CFDInitializer() {};
virtual void init(std::vector <CFDMesh*> &mesh_list, int nx, int ny, int ghost_cell) {};
virtual void updateBC(std::vector <CFDMesh*> mesh_list, int nx, int ny, int ghost_cell){};
};

class KelvinHelholtz:public CFDInitializer
{
public:
KelvinHelholtz() {};
void init(std::vector <CFDMesh*> &mesh_list, int nx, int ny, int ghost_cell) override;
void updateBC(std::vector <CFDMesh*> mesh_list, int nx, int ny, int ghost_cell) override;
};


class ObliqueShock:public CFDInitializer
{
public:
ObliqueShock() {};
void init(std::vector <CFDMesh*> &mesh_list, int nx, int ny, int ghost_cell) override;
void updateBC(std::vector <CFDMesh*> mesh_list, int nx, int ny, int ghost_cell) override;

};

class GravPressure:public CFDInitializer
{
public:
GravPressure() {};
void init(std::vector <CFDMesh*> &mesh_list, int nx, int ny, int ghost_cell) override;
void updateBC(std::vector <CFDMesh*> mesh_list, int nx, int ny, int ghost_cell) override;

};



class CFDInitializerFull
{
public:
CFDInitializerFull() {};
virtual void init(float *** &U, float xc[], float yc[], int nx, int ny, int ghost_cell, float gamma) {};
virtual void updateBC(float ***&U, int nx, int ny, int ghost_cell){};
};

class KelvinHelholtzFull:public CFDInitializerFull
{
public:
KelvinHelholtzFull() {};
void init(float *** &U, float xc[], float yc[], int nx, int ny, int ghost_cell, float gamma) override;
void updateBC(float *** &U, int nx, int ny, int ghost_cell) override;
};

class ShockTube:public CFDInitializerFull
{
public:
ShockTube() {};
void init(float *** &U, float xc[], float yc[], int nx, int ny, int ghost_cell, float gamma) override;
void updateBC(float *** &U, int nx, int ny, int ghost_cell) override;
};


class GravPressureFull:public CFDInitializerFull
{
public:
GravPressureFull() {};
void init(float *** &U, float xc[], float yc[], int nx, int ny, int ghost_cell, float gamma) override;
void updateBC(float *** &U, int nx, int ny, int ghost_cell) override;

};

class ConvectionCurrent:public CFDInitializerFull
{
public:
ConvectionCurrent() {};
void init(float *** &U, float xc[], float yc[], int nx, int ny, int ghost_cell, float gamma) override;
void updateBC(float *** &U, int nx, int ny, int ghost_cell) override;

};

class ShearStress:public CFDInitializerFull
{
public:
ShearStress() {};
void init(float *** &U, float xc[], float yc[], int nx, int ny, int ghost_cell, float gamma) override;
void updateBC(float *** &U, int nx, int ny, int ghost_cell) override;

};


#endif