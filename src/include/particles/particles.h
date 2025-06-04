#ifndef PARTICLES_H
#define PARTICLES_H

#include <glm/glm.hpp>
#include <GLFW/glfw3.h>
// #include <mesh.h>

// forward declaration of mesh so compiler knows how to work with mesh pointer in particle
class NbodyMesh;

class Particle
{
public:
    // glm::vec3 pos;
    // glm::vec3 vel;
    // glm::vec3 acc;
    glm::vec2 pos;
    glm::vec2 vel;
    glm::vec2 acc;

    glm::vec2 buffer_pos;
    glm::vec2 buffer_vel;
    glm::vec2 buffer_acc;


    float mass;
    float internal_energy;
    float ang_mom;
    NbodyMesh *mesh;
    int ind;
    float render_size;
    glm::vec3 render_color;
    Particle(int ind, glm::vec2 pos, glm::vec2 vel, glm::vec2 acc, float mass, float internal_energy, float ang_mom);
    
    void calcRenderParams();
    void updateFromBuffer();
};



#endif