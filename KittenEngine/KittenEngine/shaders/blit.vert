#version 430 core

#include "kittenCommonVert.glsl"

out vec2 uv;

void main() {
    uv = vUv;
    gl_Position = vec4(2 * vPos - 1, 1);
}