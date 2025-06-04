#ifndef PHYS_ENGINE_UNI_H
#define PHYS_ENGINE_UNI_H

#include <glm/glm.hpp>
#include <GLFW/glfw3.h>
#include <mesh.h>
// #include <particles.h>
#include <vector>
#include <iostream>
#include <future>
#include <thread_pool.h>
// #include <functional>

// Define a function pointer type
// typedef void (*FuncPtr)(PhysEngine*);

// takes mesh (which contains particles) and performs physics on it (leapfrog integration)
// includes handlers for boundary conditions
class PhysEngine
{
public:
    PhysEngine(float width, float height, int nmesh_hrz, int nmesh_vrt, int nparticles, glm::vec2 center,int threads);
    void update();
    ~PhysEngine();
    std::vector<Particle*> particles;
    Mesh* mesh;
    // NbodyMesh* mesh;
    // BarnesHutMesh* mesh;
    std::vector<int> collision_ind;
    int threads;
private:
    float width;
    float height;
    float dt;
    int nmesh_hrz, nmesh_vrt, nparticles;

    // FuncPtr eom;
    // std::function<glm::vec2(PhysEngine&)> eom;
    

    // void assign_particles();
    // void update_state();
    void kick_drift_kick(Particle* particle);
    virtual void bc_handler(Particle* particle){};
    virtual void collision_handler(Particle* particle1,Particle* particle2){};
    virtual void checkCollision(float dist_squared,int particle1_ind,int particle2_ind){};
    virtual void mesh_tree_calc(glm::vec2* acc, Particle* particle, NbodyMesh* mesh){};
    // virtual void mesh_tree_calc(glm::vec2* acc, Particle* particle, BarnesHutMesh* mesh){};
    virtual glm::vec2 calc_eom(Particle* particle){return glm::vec2(0,0);};
    void postUpdateProcess(int i);

};

class GravEngine:public PhysEngine
{
public:
    GravEngine(float width, float height, int nmesh_hrz, int nmesh_vrt, int nparticles,glm::vec2 center,int threads);    
private:
    void bc_handler(Particle* particle) override;
    void collision_handler(Particle* particle1,Particle* particle2) override;
    void mesh_tree_calc(glm::vec2* acc, Particle* particle, NbodyMesh* mesh) override;
    // void mesh_tree_calc(glm::vec2* acc, Particle* particle, BarnesHutMesh* mesh) override;
    void checkCollision(float dist_squared,int particle1_ind,int particle2_ind) override;
    // void calcForce(float dist_squared);

    glm::vec2 calc_eom(Particle* particle) override;

    float G=1;
};




#endif