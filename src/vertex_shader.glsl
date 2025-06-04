#version 330 core

layout(location = 0) in vec3 aPos;  // Vertex position attribute

void main()
{
    gl_Position = vec4(aPos, 1.0);  // Pass the position directly to the fragment shader
}
