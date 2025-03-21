#version 430 core

#include "kittenCommonVert.glsl"

out vec3 dir;

void main() {
    gl_Position = vec4(vUv.x * 2 - 1, vUv.y * 2 - 1, 1, 1);
	dir = (viewMatInv * vec4((projMatInv * vec4(gl_Position.xyz, 1)).xyz, 0)).xyz;
}