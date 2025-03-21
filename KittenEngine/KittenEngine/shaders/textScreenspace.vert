#version 430 core

#include "kittenCommonVert.glsl"

layout (binding = 3, std140) uniform uvUBO {
    vec4 atlasUvs[129];
};

layout(binding = 4, std430) buffer textBlock {
    vec4 text[];
};

out vec3 norm;
out vec3 wPos;
out vec2 uv;

void main() {
    norm = mat3(modelMat_n) * vNorm;

    vec4 data = text[gl_InstanceID];
    vec4 uvo = atlasUvs[int(data.x)];
    
    vec3 pos = vec3(vUv.x * uvo.y, -vUv.y, 0);
    pos.xy *= data.y * uvo.w;
    pos.xy += data.zw;

    gl_Position = modelMat * vec4(pos, 1);
    wPos = vec3(modelMat * vec4(pos, 1));

    uv = vec2(vUv.x * uvo.z, vUv.y * uvo.w);
    uv.x += uvo.x;
}