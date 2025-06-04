#version 330 core
in float fragDensity;
out vec4 FragColor;

void main() {
    // float intensity = (fragDensity-0.5)*100+0.5; // Normalize density for color mapping
    float intensity = fragDensity; // Normalize density for color mapping
    FragColor = vec4(intensity, 0.0, 1.0 - intensity, 1.0); // Color gradient
    // FragColor = vec4(0., 0.0, 1.0 , 1.0); // Color gradient
}
