#version 330 core
layout(location = 0) in vec2 aPos;
layout(location = 1) in float density;
layout(location = 2) in float temp;
out float fragDensity;

void main() {
    fragDensity = temp;
    // gl_PointSize = 1.0 * density * 6; // Scale point size based on density
    gl_PointSize = 1.0 * density * 4 + 6; // Scale point size based on density
    // gl_PointSize = 1.0 * density * 10; // Scale point size based on density

    // gl_PointSize = 10;
    // gl_PointSize = 2; // Scale point size based on density
    gl_Position = vec4((aPos-0.5)*1.5, 0.0, 1.0);
    // gl_Position = vec4((aPos-0.5)/4-0.5, 0.0, 1.0);
}
