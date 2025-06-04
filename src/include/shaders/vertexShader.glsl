#version 330 core
layout (location = 0) in vec2 aPos;
// layout (location = 0) in vec3 aPos;
// layout (location = 1) in float pointSize;

layout (location = 1) in vec3 pointColor;
layout (location = 2) in float pointSize;
out float particleSize;
out vec3 particleColor;
void main()
{
   // pointSize = 100;
   gl_Position = vec4(aPos/5., 0., 1.0);
   // gl_Position = vec4(0.,0.,0.,1.0);
   gl_PointSize = pointSize;
   // gl_PointSize = 10.f;
   particleSize = pointSize;
   // particleSize = 10.f;

   particleColor = pointColor;
   // particleColor = vec3(1.0f,0.0f,0.0f);
}