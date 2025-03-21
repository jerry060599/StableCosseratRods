#version 430 core

#include "kittenCommonVert.glsl"

struct GizmoData {
    mat4 model;
    vec4 color;
};

layout(binding = 3, std430) buffer dataBlock {
    GizmoData data[];
};

out vec4 col;

void main() {
    GizmoData d = data[gl_InstanceID];
    col = d.color;
    gl_Position = d.model * vec4(vPos, 1);
}