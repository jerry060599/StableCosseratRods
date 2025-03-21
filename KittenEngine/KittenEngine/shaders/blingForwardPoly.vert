#version 430 core

#include "kittenCommonVert.glsl"

out vec3 v_norm;
out vec3 v_wPos;
out vec2 v_uv;

void main() {
	v_uv = vUv;
    gl_Position = vpMat * (modelMat * vec4(vPos, 1));
    v_wPos = vec3(modelMat * vec4(vPos, 1));
}