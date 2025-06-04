#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

using namespace std;
using namespace Eigen;

// Constants
const int N = 128;          // Grid size (N x N)
const double L = 1.0;       // Domain length
const double dt = 0.0005;     // Time step
const double T_max = 2.0;   // Maximum simulation time
const double nu = 0.01;     // Kinematic viscosity
const double R = 287.0;     // Specific gas constant for air
const double gamma = 1.4;   // Ratio of specific heats

// Fields
MatrixXd u(N, N), v(N, N), p(N, N), rho(N, N), T(N, N);
MatrixXd u_new(N, N), v_new(N, N), p_new(N, N), rho_new(N, N), T_new(N, N);
MatrixXd u_source(N, N), v_source(N, N), rho_source(N, N), T_source(N, N);

// OpenGL objects
GLuint vbo, vao, shaderProgram;
GLFWwindow* window;

// const char* vertexShaderSource = R"(
// #version 330 core
// layout(location = 0) in vec2 aPos;
// out vec2 fragPos;
// void main() {
//     fragPos = aPos;
//     gl_Position = vec4(aPos, 0.0, 1.0);
// }
// )";

// const char* fragmentShaderSource = R"(
// #version 330 core
// in vec2 fragPos;
// out vec4 FragColor;
// void main() {
//     float intensity = (1.0 + fragPos.y) / 2.0;
//     FragColor = vec4(intensity, intensity, 1.0, 1.0);
// }
// )";

const char* vertexShaderSource = R"(
#version 330 core
layout(location = 0) in vec2 aPos;
layout(location = 1) in float density;
out float fragDensity;

void main() {
    fragDensity = density;
    gl_PointSize = 1.0 * density; // Scale point size based on density
    gl_Position = vec4(aPos, 0.0, 1.0);
}
)";

const char* fragmentShaderSource = R"(
#version 330 core
in float fragDensity;
out vec4 FragColor;

void main() {
    float intensity = fragDensity / 2.0; // Normalize density for color mapping
    FragColor = vec4(intensity, 0.0, 1.0 - intensity, 1.0); // Color gradient
}
)";



// Function to compile a shader
GLuint compileShader(GLenum type, const char* source) {
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);

    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetShaderInfoLog(shader, 512, nullptr, infoLog);
        cerr << "ERROR: Shader compilation failed\n" << infoLog << endl;
    }
    return shader;
}

// Function to create the shader program
GLuint createShaderProgram() {
    GLuint vertexShader = compileShader(GL_VERTEX_SHADER, vertexShaderSource);
    GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);

    GLuint program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);

    GLint success;
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetProgramInfoLog(program, 512, nullptr, infoLog);
        cerr << "ERROR: Shader program linking failed\n" << infoLog << endl;
    }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    return program;
}

// Initialize OpenGL
void initOpenGL() {
    if (!glfwInit()) {
        cerr << "Failed to initialize GLFW!" << endl;
        exit(EXIT_FAILURE);
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    window = glfwCreateWindow(800, 800, "Compressible Navier-Stokes Visualization", NULL, NULL);
    if (!window) {
        glfwTerminate();
        cerr << "Failed to create GLFW window!" << endl;
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        cerr << "Failed to initialize GLAD!" << endl;
        exit(EXIT_FAILURE);
    }

    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
    shaderProgram = createShaderProgram();

    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
}

// // Update visualization with OpenGL
// void visualize() {
//     vector<float> vertices;
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < N; j++) {
//             float x = j * 2.0f / N - 1.0f;
//             float y = i * 2.0f / N - 1.0f;
//             float intensity = T(i, j);
//             vertices.push_back(x);
//             vertices.push_back(y);
//         }
//     }

//     glBindVertexArray(vao);
//     glBindBuffer(GL_ARRAY_BUFFER, vbo);
//     glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_DYNAMIC_DRAW);

//     glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
//     glEnableVertexAttribArray(0);

//     glUseProgram(shaderProgram);
//     glClear(GL_COLOR_BUFFER_BIT);
//     glDrawArrays(GL_POINTS, 0, N * N);

//     glfwSwapBuffers(window);
//     glfwPollEvents();
// }

void visualize() {
    vector<float> vertices;
    vector<float> densities;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float x = j * 2.0f / N - 1.0f;
            float y = i * 2.0f / N - 1.0f;
            vertices.push_back(x);
            vertices.push_back(y);
            densities.push_back(rho(i, j)); // Add density for rendering
        }
    }

    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    GLuint densityBuffer;
    glGenBuffers(1, &densityBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, densityBuffer);
    glBufferData(GL_ARRAY_BUFFER, densities.size() * sizeof(float), densities.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);

    glUseProgram(shaderProgram);
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawArrays(GL_POINTS, 0, N * N);

    glfwSwapBuffers(window);
    glfwPollEvents();
}

// Function to initialize fields
// void initializeFields() {
//     u.setZero();
//     v.setZero();
//     p.setZero();
//     rho.setConstant(1.0); // Initial density
//     T.setConstant(300.0); // Initial temperature

//     for (int i = 0; i < N; i++) {
//         double y = i * L / N;
//         for (int j = 0; j < N; j++) {
//             double x = j * L / N;

//             if (y < 0.25 || y > 0.75) {
//                 u(i, j) = -0.5;  // Lower and upper regions
//             } else {
//                 u(i, j) = 0.5;   // Middle region
//             }

//             v(i, j) = 0.01 * sin(2.0 * M_PI * x); // Small perturbation
//         }
//     }

//     p = rho.cwiseProduct(T) * R; // Initial pressure using ideal gas law
// }

void initializeFields() {
    u.setZero();
    v.setZero();
    p.setZero();
    T.setConstant(300.0); // Initial temperature

    for (int i = 0; i < N; i++) {
        double y = i * L / N;
        for (int j = 0; j < N; j++) {
            double x = j * L / N;

            if (y < 0.25) {
                u(i, j) = -0.5; // Lower region
                rho(i, j) = 1.2; // Higher density
            } else if (y > 0.75) {
                u(i, j) = -0.5; // Upper region
                rho(i, j) = 1.2; // Lower density
            } else {
                u(i, j) = 0.5;  // Middle region
                rho(i, j) = 1.; // Intermediate density
            }

            v(i, j) = 0.01 * sin(2.0 * M_PI * x); // Small perturbation
        }
    }

    p = rho.cwiseProduct(T) * R; // Initial pressure using ideal gas law
}


// Apply periodic boundary conditions
void applyBoundaryConditions() {
    for (int i = 0; i < N; i++) {
        u(i, 0) = u(i, N - 2);
        u(i, N - 1) = u(i, 1);

        v(i, 0) = v(i, N - 2);
        v(i, N - 1) = v(i, 1);

        rho(i, 0) = rho(i, N - 2);
        rho(i, N - 1) = rho(i, 1);

        T(i, 0) = T(i, N - 2);
        T(i, N - 1) = T(i, 1);
    }

    for (int j = 0; j < N; j++) {
        u(0, j) = u(N - 2, j);
        u(N - 1, j) = u(1, j);

        v(0, j) = v(N - 2, j);
        v(N - 1, j) = v(1, j);

        rho(0, j) = rho(N - 2, j);
        rho(N - 1, j) = rho(1, j);

        T(0, j) = T(N - 2, j);
        T(N - 1, j) = T(1, j);
    }
}

// Compute advection
void computeAdvection() {
    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            double dx = L / N;

            u_source(i, j) = u(i, j) - dt * (
                (u(i, j) * (u(i, j) - u(i, j - 1)) / dx) +
                (v(i, j) * (u(i, j) - u(i - 1, j)) / dx));

            v_source(i, j) = v(i, j) - dt * (
                (u(i, j) * (v(i, j) - v(i, j - 1)) / dx) +
                (v(i, j) * (v(i, j) - v(i - 1, j)) / dx));

            rho_source(i, j) = rho(i, j) - dt * (
                u(i, j) * (rho(i, j) - rho(i, j - 1)) / dx +
                v(i, j) * (rho(i, j) - rho(i - 1, j)) / dx);

            T_source(i, j) = T(i, j) - dt * (
                u(i, j) * (T(i, j) - T(i, j - 1)) / dx +
                v(i, j) * (T(i, j) - T(i - 1, j)) / dx);
        }
    }
}

// Compute pressure and diffusion
void computeDiffusionAndPressure() {
    double dx = L / N;
    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            u_new(i, j) = u_source(i, j) + dt * nu * (
                (u(i + 1, j) - 2.0 * u(i, j) + u(i - 1, j)) / (dx * dx) +
                (u(i, j + 1) - 2.0 * u(i, j) + u(i, j - 1)) / (dx * dx));

            v_new(i, j) = v_source(i, j) + dt * nu * (
                (v(i + 1, j) - 2.0 * v(i, j) + v(i - 1, j)) / (dx * dx) +
                (v(i, j + 1) - 2.0 * v(i, j) + v(i, j - 1)) / (dx * dx));

            p_new(i, j) = (gamma - 1.0) * rho_source(i, j) * T_source(i, j); // Update pressure
        }
    }
}

// Main simulation loop
void simulate() {
    double time = 0.0;
    initializeFields();

    while (time < T_max && !glfwWindowShouldClose(window)) {
        applyBoundaryConditions();
        computeAdvection();
        computeDiffusionAndPressure();

        // Swap buffers
        u = u_new;
        v = v_new;
        p = p_new;
        rho = rho_source;
        T = T_source;

        time += dt;

        visualize();

        // Output progress
        cout << "Time: " << time << " / " << T_max << endl;

        cout << rho.maxCoeff() << endl;
    }

    glfwTerminate();
}

int main() {
    initOpenGL();
    simulate();
    cout << "Simulation completed!" << endl;
    return 0;
}
