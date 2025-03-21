#version 430 core

layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

in vec3 v_wPos[];
in vec2 v_uv[];

out vec3 norm;
out vec2 uv;
out vec3 wPos;

uniform mat4 uModelMatrix;

void main() {
    vec3 edge1 = v_wPos[1] - v_wPos[0];
    vec3 edge2 = v_wPos[2] - v_wPos[0];
    vec3 normal = normalize(cross(edge1, edge2));

    for (int i = 0; i < 3; ++i) {
        norm = normal;
        uv = v_uv[i];
        wPos = v_wPos[i];
        gl_Position = gl_in[i].gl_Position;

        EmitVertex();
    }

    EndPrimitive();
}