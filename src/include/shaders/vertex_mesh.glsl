#version 330 core
layout (location = 0) in vec2 aPos;

void main() {
    gl_Position = vec4(aPos/5.,0., 1.0);
    gl_PointSize = 4.f;

}
