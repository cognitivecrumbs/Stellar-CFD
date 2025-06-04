#ifndef GLFW_H
#define GLFW_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>

class GLFW_handler
{
public:
    GLFWwindow* window;
    int SCR_WIDTH,SCR_HEIGHT;

    GLFW_handler(int SCR_WIDTH, int SCR_HEIGHT);
    
    int init_screen();
    // static void framebuffer_size_callback(GLFWwindow* window, int width, int height);

private:
    static void framebuffer_size_callback(GLFWwindow* window, int width, int height);
    
};

#endif