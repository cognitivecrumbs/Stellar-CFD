#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <filesystem>

#include <chrono>
// #include <shader_class.h>
#include <glfw_handler.h>

#include <render.h>
#include <physics_engine.h>
// #include <physics_engine_uniform.h>


// settings
const unsigned int SCR_WIDTH = 800;
// const unsigned int SCR_HEIGHT = 600;
const unsigned int SCR_HEIGHT = 800;

// bool to render mesh
bool mesh_render=false;

void processInput(GLFWwindow *window);

int main()
{
    std::filesystem::path cwd = std::filesystem::current_path();

    std::cout << "Current working directory: " << cwd << std::endl;

    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();

    GLFW_handler glfw_handler = GLFW_handler(SCR_WIDTH,SCR_HEIGHT);
    int glfw_handler_return = glfw_handler.init_screen();
    if (glfw_handler_return == -1){
        std::cout << "Screen initialization failed" << std::endl;
        return -1;
    }
    GLFWwindow* window = glfw_handler.window;

    RenderHandler render = RenderHandler("src/include/shaders/vertexShader.glsl","src/include/shaders/fragShader.glsl",glm::vec3(0.1,0.1,0.1));

    // GravEngine engine = GravEngine(SCR_WIDTH*4, SCR_HEIGHT*4, 1, 1, 2, glm::vec2(0,0));
    // GravEngine engine = GravEngine(SCR_WIDTH, SCR_HEIGHT, 1, 1, 1000, glm::vec2(0,0));
    // GravEngine engine = GravEngine(9, 9, 51, 51, 4000, glm::vec2(0,0),8); //15 fps with barnes hut, 6 with uniform no collisions
    // GravEngine engine = GravEngine(9, 9, 51, 51, 10000, glm::vec2(0,0),8); 
    // GravEngine engine = GravEngine(9, 9, 51, 51, 2000, glm::vec2(0,0),8); 
    GravEngine engine = GravEngine(15, 15, 51, 51, 2000, glm::vec2(0,0),8); 
    // GravEngine engine = GravEngine(2.5, 2.5, 2, 2, 3, glm::vec2(0,0),8);
    render.initBuffers(&engine.particles);

    // mesh render 
    // MeshRenderHandler mesh_renderer;
    MeshRenderHandler mesh_renderer = MeshRenderHandler("src/include/shaders/vertex_mesh.glsl","src/include/shaders/frag_mesh.glsl",&engine);
    // mesh render stup end

    constexpr float MAX_FPS = 60.f;
    constexpr float FRAME_TIME = 1.f / MAX_FPS;
    // float elapsed_ms;
    int i =0;
    float last_time = glfwGetTime();
    float last_update_time = glfwGetTime();
    float accumulated_time = 0.f;

    while (!glfwWindowShouldClose(window))
    {   
        float current_time = glfwGetTime();
        float dt = current_time - last_time;
        last_time = current_time;
        accumulated_time += dt;


        if (accumulated_time >= FRAME_TIME)
            {
                float update_dt = current_time - last_update_time;
                last_update_time = current_time;
                accumulated_time = 0.f;

            i ++ ;
            // Process user input

            std::cout << "iteration: " << i << std::endl;
            processInput(window);
            // glEnable(GL_PROGRAM_POINT_SIZE);

            engine.update();

            // std::cout << engine.particles[0]->pos.x << ", " << engine.particles[0]->pos.y << "|" << engine.particles[1]->pos.x << ", " << engine.particles[1]->pos.y << std::endl;
            // glm::vec2 tot_mom = engine.particles[0]->vel*engine.particles[0]->mass + engine.particles[1]->vel*engine.particles[1]->mass;
            // std::cout << engine.particles[1]->pos.x << ", " << engine.particles[1]->pos.y << std::endl;
            // std::cout << engine.particles[0]->render_size << ", " << engine.particles[0]->render_size << std::endl;
            // std::cout << engine.particles[0] << ", " << engine.particles[1] << std::endl;

            // glm::vec2 rel_vel = engine.particles[0]->vel - engine.particles[1]->vel;
            // glm::vec2 rel_pos = engine.particles[0]->pos - engine.particles[1]->pos;

            // // float energy = dot(rel_vel,rel_vel)/2.0f - (engine.particles[0]->mass+engine.particles[1]->mass)/sqrt(dot(rel_pos,rel_pos));
            // // std::cout << energy << std::endl;

            // float kinetic_energy = 0.5f * engine.particles[0]->mass * dot(engine.particles[0]->vel, engine.particles[0]->vel)
            //          + 0.5f * engine.particles[1]->mass * dot(engine.particles[1]->vel, engine.particles[1]->vel);

            // float distance = sqrt(dot(rel_pos, rel_pos));
            // float potential_energy = -(1.0f * engine.particles[0]->mass * engine.particles[1]->mass) / distance;

            // float total_energy = kinetic_energy + potential_energy;
            // std::cout << total_energy << std::endl;
    
            // if (i%2==0){
            render.updateBuffers();
            render.render(window,&engine);

            if (mesh_render){
                mesh_renderer.render();
            }

            // Swap buffers and poll events
            glfwSwapBuffers(window);
            glfwPollEvents();
            }
        
    }


    // Cleanup
    render.terminateClean();
    mesh_renderer.terminateClean();


    glfwTerminate();
    return 0;
}


void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}
