#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <filesystem>

#include <chrono>
#include <glfw_handler.h>

#include <render.h>
#include <physics_engine.h>
// #include <mesh.h>

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 800;

// bool to render mesh
bool mesh_render=false;

void processInput(GLFWwindow *window);




int main()
{
    float l = 1;
    int nx = 128; 
    int ny = 128;
    // nx = 240; ny = 240;

    int ghost = 2;
    CFDMeshWrapper mesh_wrapper = CFDMeshWrapper(nx,ny,ghost,l/nx);

    
    // std::filesystem::path cwd = std::filesystem::current_path();

    // std::cout << "Current working directory: " << cwd << std::endl;

    // std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();

    GLFW_handler glfw_handler = GLFW_handler(SCR_WIDTH,SCR_HEIGHT);
    int glfw_handler_return = glfw_handler.init_screen();
    if (glfw_handler_return == -1){
        std::cout << "Screen initialization failed" << std::endl;
        return -1;
    }
    GLFWwindow* window = glfw_handler.window;
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    Shader shader = Shader("src/include/shaders/vertex_cfd_mesh.glsl","src/include/shaders/frag_cfd_mesh.glsl");

    unsigned int vao,vbo;

    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);

    float last_time = glfwGetTime();
    float last_update_time = glfwGetTime();
    float accumulated_time = 0.f;

    std::chrono::steady_clock::time_point lastTime;
    int frameCount = 0;
    float fps = 0.0f;

    int i = 0;
    while (!glfwWindowShouldClose(window))
    {   
        float current_time = glfwGetTime();
        float dt = current_time - last_time;
        last_time = current_time;
        accumulated_time += dt;

        // glfwWaitEvents();  // Wait for events (this pauses the loop)

        // // Check if the spacebar is pressed
        // if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
        // {
        //     std::cout << "Spacebar pressed!" << std::endl;
        //     // break; // Break the loop after spacebar is pressed
        // }
        // if (i==7){
        //     std::cout << "before nan" << endl;
        // }
        // std::cout << "iteration: " << i << endl;
        i++;

        vector<glm::vec2> vertices;
        vector<float> densities;
        vector<float> temperature;

        // for (int i = 0; i < N; i++) {
        //     for (int j = 0; j < N; j++) {
        //         float x = j * 2.0f / N - 1.0f;
        //         float y = i * 2.0f / N - 1.0f;
        //         vertices.push_back(x);
        //         vertices.push_back(y);
        //         densities.push_back(rho(i, j)); // Add density for rendering
        //     }
        // }
        mesh_wrapper.update_mesh();

        float update_dt = current_time - last_update_time;
        last_update_time = current_time;
        accumulated_time = 0.f;


        frameCount++;

        // Update FPS once per second
        auto currentTime = std::chrono::steady_clock::now();
        float timeInterval = std::chrono::duration<float>(currentTime - lastTime).count();

        if (timeInterval >= 1.0f) {
            fps = frameCount / timeInterval;
            lastTime = currentTime;
            frameCount = 0;

            // Update the window title with the current FPS
            std::ostringstream title;
            title << "FPS: " << std::fixed << std::setprecision(2) << fps;
            glfwSetWindowTitle(window, title.str().c_str());
        }

        for (CFDMesh* mesh: mesh_wrapper.mesh){
            vertices.push_back(mesh->center);
            densities.push_back(mesh->U[0]);
            // temperature.push_back((mesh->temp-50)/70);
            // temperature.push_back(mesh->vel_squared*10);
            // temperature.push_back(mesh->prim[1]*50);
            // temperature.push_back(mesh->U[3]);

            float p = (mesh->gamma - 1) * (mesh->U[3] - 0.5 * (mesh->U[1] * mesh->U[1] + mesh->U[2] * mesh->U[2]) / mesh->U[0]);
            // temperature.push_back((mesh->U[0]-2));

            float u = mesh->U[1]/mesh->U[0];
            float v = mesh->U[2]/mesh->U[0];

            temperature.push_back((mesh->U[3]-0.5*(u*u-v*v)-5)/8);

            // temperature.push_back((mesh->U[2]/mesh->U[0])*100+1);



            // temperature.push_back((mesh->U[1]/mesh->U[0]+0.5));
            // temperature.push_back(float(mesh->ind)/((nx+2*ghost)*(ny+2*ghost)));

            // temperature.push_back(mesh->U[1]/mesh->U[0]/10);

            // cout << mesh->U[0] << ", ";
        }
        // cout << endl;

        glBindVertexArray(vao);

        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec2), vertices.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (void*)0);
        glEnableVertexAttribArray(0);

        GLuint densityBuffer;
        glGenBuffers(1, &densityBuffer);
        glBindBuffer(GL_ARRAY_BUFFER, densityBuffer);
        glBufferData(GL_ARRAY_BUFFER, densities.size() * sizeof(float), densities.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
        glEnableVertexAttribArray(1);

        GLuint temperatureBuffer;
        glGenBuffers(1, &temperatureBuffer);
        glBindBuffer(GL_ARRAY_BUFFER, temperatureBuffer);
        glBufferData(GL_ARRAY_BUFFER, temperature.size() * sizeof(float), temperature.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
        glEnableVertexAttribArray(2);

        shader.use();
        glClear(GL_COLOR_BUFFER_BIT);
        glDrawArrays(GL_POINTS, 0, (nx+2*ghost) * (ny+2*ghost));

    

        // Swap buffers and poll events
        glfwSwapBuffers(window);
        glfwPollEvents();
            
        
    }


    // // Cleanup
    // render.terminateClean();
    // mesh_renderer.terminateClean();


    // glfwTerminate();
    // return 0;
}


void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}
