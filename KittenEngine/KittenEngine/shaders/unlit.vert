#version 430 core

#include "kittenCommonVert.glsl"

void main() {
    gl_Position = vpMat * (modelMat * vec4(vPos, 1));
}