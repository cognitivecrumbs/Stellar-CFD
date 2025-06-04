#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <glad/glad.h>  // Include GLAD header
#include <glm/glm.hpp>
#include <GLFW/glfw3.h>
#include <glfw_handler.h>


using namespace std;
using namespace Eigen;

// Constants
const double v = 1.0;  // Velocity field (constant)
const double L = 1.0;  // Length of the domain
const int N = 100;     // Number of grid points
const double T_max = 1000.0;  // Max simulation time
const double CFL = 0.5; // Courant number for stability

// Function to initialize the solution (Gaussian pulse in the center)
void initialize(VectorXd& u) {
    for (int i = 0; i < N; i++) {
        double x = double(i) * L / N;
        u[i] = exp(-100.0 * (x - 0.5) * (x - 0.5));
    }
}

// Function to update the solution using finite volume method
void update(VectorXd& u, double dt, double dx) {
    VectorXd u_new = u;
    
    // Apply upwind scheme for the advection term
    for (int i = 1; i < N-1; i++) {
        u_new[i] = u[i] - CFL * (u[i] - u[i-1]);  // Upwind scheme
    }
    
    // Boundary conditions: Periodic boundary conditions
    u_new[0] = u_new[N-2];  // Left boundary
    u_new[N-1] = u_new[1];  // Right boundary

    // Update the solution
    u = u_new;
}

// Function to write the solution to a file
void writeSolution(const VectorXd& u, double time, const string& filename) {
    ofstream file(filename, ios::app);
    file << time;
    for (int i = 0; i < N; i++) {
        file << " " << u[i];
    }
    file << endl;
}

// OpenGL Shader program
GLuint vbo, vao, shaderProgram;
GLFWwindow* window;

const char* vertexShaderSource = R"(
#version 330 core
layout(location = 0) in vec3 aPos;  // Vertex position attribute
void main()
{
    gl_PointSize = 5;
    gl_Position = vec4(aPos, 1.0);  // Pass the position directly to the fragment shader
}
)";

const char* fragmentShaderSource = R"(
#version 330 core
out vec4 FragColor;  // Output color
void main()
{
    FragColor = vec4(1.0, 1.0, 1.0, 1.0);  // Set color to white
}
)";

// Function to compile a shader and check for errors
GLuint compileShader(GLenum type, const char* source) {
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);
    
    // Check for compilation errors
    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetShaderInfoLog(shader, 512, nullptr, infoLog);
        std::cerr << "ERROR: Shader compilation failed\n" << infoLog << std::endl;
    }
    return shader;
}

// Function to create a shader program
GLuint createShaderProgram() {
    GLuint vertexShader = compileShader(GL_VERTEX_SHADER, vertexShaderSource);
    GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);
    
    // Link the shaders into a program
    GLuint program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);
    
    // Check for linking errors
    GLint success;
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetProgramInfoLog(program, 512, nullptr, infoLog);
        std::cerr << "ERROR: Shader program linking failed\n" << infoLog << std::endl;
    }
    
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    
    return program;
}

// Function to initialize OpenGL
void initOpenGL() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW!" << std::endl;
        exit(EXIT_FAILURE);
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    #ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif
    // Create an OpenGL window with GLFW
    window = glfwCreateWindow(800, 600, "FVM Visualization", NULL, NULL);
    if (!window) {
        glfwTerminate();
        std::cerr << "Failed to create GLFW window!" << std::endl;
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);

    // Initialize GLAD to load OpenGL extensions
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Failed to initialize GLAD!" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "OpenGL version: " << glGetString(GL_VERSION) << std::endl;
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    // Create the shader program
    shaderProgram = createShaderProgram();

    // Set up the buffer objects
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
}

// Function to draw the solution using OpenGL
void drawOpenGL(const VectorXd& u) {
    glClear(GL_COLOR_BUFFER_BIT);

    // Map the vector data to a GLfloat array
    GLfloat* data = new GLfloat[N * 3];  // 3 for x, y, z positions
    for (int i = 0; i < N; i++) {
        // The x position corresponds to the index, and the y position corresponds to u[i]
        data[i * 3] = (i * L / N) * 2.0f - 1.0f;  // Normalize to [-1, 1] for x
        data[i * 3 + 1] = u[i] * 2.0f - 1.0f;    // Normalize u values to [-1, 1] for y
        data[i * 3 + 2] = 0.0f;                  // z = 0 (2D)
    }

    // Bind the VAO and VBO
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * N * 3, data, GL_DYNAMIC_DRAW);

    // Set the vertex position attribute pointer
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);
    glEnableVertexAttribArray(0);

    // Use the shader program
    glUseProgram(shaderProgram);

    // Draw the points
    glDrawArrays(GL_POINTS, 0, N);

    delete[] data;

    glfwSwapBuffers(window);
    glfwPollEvents();
}


int main() {
    // Initialize OpenGL
    initOpenGL();

    // Discretization
    double dx = L / N;             // Grid spacing
    double dt = CFL * dx / v;     // Time step based on CFL condition
    int steps = int(T_max / dt);  // Number of time steps

    // Create the vector to store solution
    VectorXd u(N);
    initialize(u);  // Initialize solution

    // Write initial condition to file
    writeSolution(u, 0.0, "solution.dat");

    // OpenGL render loop
    int timeStep = 0;
    while (!glfwWindowShouldClose(window) && timeStep <= steps) {
        // Update solution
        update(u, dt, dx);
        
        // Write solution to file every 10 steps
        if (timeStep % 10 == 0) {
            writeSolution(u, timeStep * dt, "solution.dat");
        }

        // Visualize the solution
        drawOpenGL(u);

        timeStep++;
    }

    glfwTerminate();
    std::cout << "Simulation finished!" << std::endl;
    return 0;
}