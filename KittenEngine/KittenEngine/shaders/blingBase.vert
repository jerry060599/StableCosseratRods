#version 430 core

#include "kittenCommonVert.glsl"

out vec3 norm;
out vec3 wPos;
out vec2 uv;

void main() {
	uv = vUv;
    gl_Position = vpMat * (modelMat * vec4(vPos, 1));
    norm = mat3(modelMat_n) * vNorm;
    wPos = vec3(modelMat * vec4(vPos, 1));
}