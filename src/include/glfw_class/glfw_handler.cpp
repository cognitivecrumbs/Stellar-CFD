#include <glfw_handler.h>


GLFW_handler::GLFW_handler(int SCR_WIDTH, int SCR_HEIGHT)
    : SCR_WIDTH(SCR_WIDTH), SCR_HEIGHT(SCR_HEIGHT) {}
    // {this->SCR_WIDTH = SCR_WIDTH; this->SCR_HEIGHT=SCR_HEIGHT;}

int GLFW_handler::init_screen(){
    // glfw: initialize and configure
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    #ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif

    // glfw window creation
    this->window = glfwCreateWindow(this->SCR_WIDTH, this->SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (this->window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(this->window);
    // glfwSetFramebufferSizeCallback(window, this->framebuffer_size_callback);
    
    // Associate the window with this instance
    // glfwSetWindowUserPointer(window, this);

    // Set the static callback
    glfwSetFramebufferSizeCallback(this->window, framebuffer_size_callback);


    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }
    return 0;
    }

void GLFW_handler::framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
    // glViewport(0, 0, 100, 100);
    // return;
}
