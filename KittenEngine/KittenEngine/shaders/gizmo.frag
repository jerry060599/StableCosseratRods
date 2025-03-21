#version 430 core

#include "kittenCommonFrag.glsl"

out vec4 fragColor;
in vec4 col;

void main() {
	fragColor = col;
}