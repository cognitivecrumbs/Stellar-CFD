#include <render.h>

RenderHandler::RenderHandler(const char* vertexPath, const char* fragmentPath, glm::vec3 background_color):
shader(Shader(vertexPath,fragmentPath)), background_color(background_color)
{
    // this->shader = Shader(vertexPath,fragmentPath);
}

void RenderHandler::initBuffers(std::vector<Particle*>* particles_ptr)
{
    this->particles_ptr = particles_ptr;
    // initialize vector for storage
    // this->positions = std::vector<glm::vec2> positions(glm::vec2(0,0));
    // this->color = std::vector<glm::vec2> color(glm::vec3(0,0,0));
    // this->size = std::vector<float> size(0);
    this->positions.resize(this->particles_ptr->size(),glm::vec2(0,0));
    this->colors.resize(this->particles_ptr->size(),glm::vec3(0,0,0));
    this->sizes.resize(this->particles_ptr->size(),0.f);

    this->updateRenderParams();

    std::cout << this->positions[0].x << ", " << this->positions[0].y << "|" <<this->positions[1].x << ", " << this->positions[1].y << std::endl;
    std::cout << this->sizes[0] << ", " << this->sizes[1] << std::endl;

    // unsigned int VAO,VBO_pos,VBO_color,VBO_size;
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    glGenVertexArrays(1, &VAO);
    // glGenVertexArrays(1, &VAO2);
    glGenBuffers(1, &VBO_pos);
    glGenBuffers(1, &VBO_color);
    glGenBuffers(1, &VBO_size);

    // glGenBuffers(1, &VBO2);

    // Points setup
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_pos);
    // glBufferData(GL_ARRAY_BUFFER, sizeof(this->positions), this->positions.data(), GL_DYNAMIC_DRAW);
    glBufferData(GL_ARRAY_BUFFER, this->positions.size()*sizeof(glm::vec2), this->positions.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (void*)0);
    glEnableVertexAttribArray(0);

    
    glBindBuffer(GL_ARRAY_BUFFER, VBO_color);
    // glBufferData(GL_ARRAY_BUFFER, sizeof(this->colors), this->colors.data(), GL_DYNAMIC_DRAW);
    glBufferData(GL_ARRAY_BUFFER, this->colors.size()*sizeof(glm::vec3), this->colors.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(1);


    // Create separate buffer for point sizes
    // unsigned int VBO_point_size;
    // glGenBuffers(1, &VBO_point_size);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_size);
    // glBufferData(GL_ARRAY_BUFFER, sizeof(point_size), point_size, GL_STATIC_DRAW);
    glBufferData(GL_ARRAY_BUFFER, sizeof(this->sizes), this->sizes.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
    glEnableVertexAttribArray(2);

    // clear buffer binds
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void RenderHandler::updateBuffers(){
    this -> updateRenderParams();

    // Update positions
    glBindBuffer(GL_ARRAY_BUFFER, this->VBO_pos);
    glBufferData(GL_ARRAY_BUFFER, this->positions.size() * sizeof(glm::vec2), this->positions.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (void*)0);
    glEnableVertexAttribArray(0);

    // Update colors
    glBindBuffer(GL_ARRAY_BUFFER, this->VBO_color);
    glBufferData(GL_ARRAY_BUFFER, this->colors.size() * sizeof(glm::vec3), this->colors.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(1);

    // Update sizes
    glBindBuffer(GL_ARRAY_BUFFER, this->VBO_size);
    glBufferData(GL_ARRAY_BUFFER, this->sizes.size() * sizeof(float), this->sizes.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
    glEnableVertexAttribArray(2);


}

void RenderHandler::changeEngineParticles(std::vector<Particle*>* particles_ptr){
    this->particles_ptr = particles_ptr;
}

void RenderHandler::render(GLFWwindow* window,PhysEngine* engine){
    glClearColor(this->background_color.r,this->background_color.b,this->background_color.g, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    this->draw();

    // Calculate deltaTime for FPS
    static auto lastTime = std::chrono::steady_clock::now();
    auto currentTime = std::chrono::steady_clock::now();
    float deltaTime = std::chrono::duration<float>(currentTime - lastTime).count();
    lastTime = currentTime;

    // Update FPS in window title
    updateFPS(deltaTime, window, engine);

}
void RenderHandler::draw(){
    this->shader.use();

    // bind buffer
    glBindVertexArray(this->VAO);
    // draw points
    glDrawArrays(GL_POINTS, 0, this->particles_ptr->size());
}

void RenderHandler::updateRenderParams(){
    for (size_t i=0; i<this->particles_ptr->size(); i++){
        if ((*this->particles_ptr)[i] != nullptr){
            positions[i] = (*this->particles_ptr)[i]->pos;
            colors[i]    = (*this->particles_ptr)[i]->render_color;
            sizes[i]     = (*this->particles_ptr)[i]->render_size;
        }
        else {
            sizes[i]     = 0;
            colors[i]    = this->background_color;
        }
    }
    
}

void RenderHandler::updateFPS(float deltaTime, GLFWwindow* window, PhysEngine* engine) {
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
        title << "FPS: " << std::fixed << std::setprecision(2) << fps << ", Nbodies: " << std::fixed << std::setprecision(2) << engine->mesh->npart_active;
        glfwSetWindowTitle(window, title.str().c_str());
    }
}

void RenderHandler::terminateClean(){
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO_pos);
    glDeleteBuffers(1, &VBO_color);
    glDeleteBuffers(1, &VBO_size);
    // glDeleteVertexArrays(1, &VAO2);
    // glDeleteBuffers(1, &VBO2);
    // glDeleteBuffers(1, &VBO_point_size);
    glDeleteProgram(shader.ID);
}

MeshRenderHandler::MeshRenderHandler(const char* vertexPath, const char* fragmentPath,PhysEngine* engine):
shader(Shader(vertexPath,fragmentPath))
{
    // this -> shader = Shader(vertexPath,fragmentPath);

    std::vector<std::vector<glm::vec2>> vertices;
    std::vector<glm::vec2> com;

    this -> engine = engine;
    this -> engine -> mesh->return_corners(&vertices,&com);

    // OpenGL buffers
    // unsigned int VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glGenVertexArrays(1, &VAO_com);
    glGenBuffers(1, &VBO_com); 
}

void MeshRenderHandler::render(){
    this->com.clear();
    this->vertices.clear();
    this->engine->mesh->return_corners(&this->vertices,&this->com);

    this->shader.use();
    glBindVertexArray(VAO);
    for (const auto& shape : this->vertices) {
        // Update vertex buffer
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, shape.size() * sizeof(glm::vec2), shape.data(), GL_DYNAMIC_DRAW);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);


        // Draw shape as a line loop
        glDrawArrays(GL_LINE_LOOP, 0, shape.size());
    }


    glBindVertexArray(VAO_com);
    // Update vertex buffer
    glBindBuffer(GL_ARRAY_BUFFER, VBO_com);
    glBufferData(GL_ARRAY_BUFFER, this->com.size() * sizeof(glm::vec2), this->com.data(), GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glPointSize(100.0f);
    glDrawArrays(GL_POINTS, 0, this->com.size());

    glBindVertexArray(0);
}


void MeshRenderHandler::terminateClean(){
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shader.ID);
}