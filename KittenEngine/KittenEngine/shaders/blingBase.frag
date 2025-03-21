#version 430 core

#include "kittenCommonFrag.glsl"
#include "kittenLighting.glsl"

in vec2 uv;
in vec3 norm;
in vec3 wPos;
out vec4 fragColor;

void main() {
	fragColor = vec4(getLightCol(), 1) * matColor * texture(tex_d, uv);
}