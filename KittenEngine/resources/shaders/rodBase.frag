#version 430 core

#include "kittenCommonFrag.glsl"
#include "kittenLighting.glsl"

out vec4 fragColor;
in vec3 col;

void main() {
	fragColor = vec4(matColor.xyz * col, 1);
}