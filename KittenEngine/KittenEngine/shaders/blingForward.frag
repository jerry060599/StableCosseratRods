#version 430 core

#include "kittenCommonFrag.glsl"
#include "kittenLighting.glsl"

in vec2 uv;
in vec3 norm;
in vec3 wPos;
out vec4 fragColor;

void main() {
	vec3 n = norm / length(norm);
	vec3 lDir = getLightDir(wPos);
	vec3 halfVec = 0.5 * (-lDir + getViewDir(wPos));

	float d = max(dot(-lDir, n), 0);
	float s = pow(max(dot(halfVec, n), 0), matParams0.x);

	fragColor = vec4(getLightCol(wPos) * (
					 matColor1.xyz * texture(tex_s, uv).xyz * s
				   + matColor.xyz * texture(tex_d, uv).xyz * d
				   	), 0);
}