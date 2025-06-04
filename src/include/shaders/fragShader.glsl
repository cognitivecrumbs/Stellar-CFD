#version 330 core
out vec4 FragColor;
in float particleSize;
in vec3 particleColor;

void main()
{
   // convert point to circle PointCoord.xy gives display points normalized between 0,1 for point
   // 0.5 radius squared -> 0.25
   vec2 pos = gl_PointCoord.xy-0.5;
   if(dot(pos,pos) > 0.25 && particleSize > 3.0){
      discard;
   }
   else if (particleSize == 0){
      discard;
   }
   FragColor = vec4(particleColor, 1.0f);
   // FragColor = vec4(1.0f, 0.f, 0.f, 1.0f);
}