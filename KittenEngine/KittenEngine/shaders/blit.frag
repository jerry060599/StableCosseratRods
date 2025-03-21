#version 430 core

#include "kittenCommonFrag.glsl"

in vec2 uv;
out vec4 fragColor;

void main() {
	fragColor = matColor * texture(tex_d, uv);
}