#ifndef RENDER_H
#define RENDER_H

#include <shader_class.h>
#include <GLFW/glfw3.h>
#include <vector>
// #include <particles.h>
#include <physics_engine.h>


// takes mesh (which contains particles) and performs physics on it (leapfrog integration)
// includes handlers for boundary conditions
class RenderHandler
{
public:
    RenderHandler(const char* vertexPath, const char* fragmentPath, glm::vec3 background_color);
    // void update();
    // ~RenderHandler();
    void initBuffers(std::vector<Particle*>* particles_ptr);
    void updateBuffers();

    void render(GLFWwindow* window,  PhysEngine* engine);
    void updateFPS(float deltaTime, GLFWwindow* window, PhysEngine* engine); // New method to update FPS
    
    void terminateClean();

    std::chrono::steady_clock::time_point lastTime;
    int frameCount = 0;
    float fps = 0.0f;

    glm::vec3 background_color;
private:
    // PhysEngine* engine;
    unsigned int VAO,VBO_pos,VBO_color,VBO_size;
    Shader shader;
    std::vector<Particle*>* particles_ptr;

    std::vector<glm::vec2> positions;
    std::vector<glm::vec3> colors;
    std::vector<float> sizes;

    void updateRenderParams();
    void changeEngineParticles(std::vector<Particle*>* particles_ptr);
    void draw();


};

class MeshRenderHandler
{
public:
    MeshRenderHandler(const char* vertexPath, const char* fragmentPath, PhysEngine* engine);

    // updates buffers and draws
    void render();
    void terminateClean();
    
private:
    // PhysEngine* engine;
    std::vector<std::vector<glm::vec2>> vertices;
    std::vector<glm::vec2> com;
    Shader shader;
    unsigned int VAO,VBO,VAO_com,VBO_com;
    PhysEngine* engine;

};


#endif