
layout (binding = 0, std140) uniform globalUBO {
	mat4 projMat;
	mat4 projMatInv;
	mat4 viewMat;
	mat4 viewMatInv;
	mat4 vpMat;
	mat4 vpMatInv;
	mat4 viewMat_n;
};

layout (binding = 1, std140) uniform modelUBO {
	mat4 modelMat;
	mat4 modelMatInv;
	mat4 modelMat_n;
};

layout (binding = 2, std140) uniform materialUBO {
	vec4 matColor;
	vec4 matColor1;
	vec4 matParams0;
	vec4 matParams1;
};

#include "kittenUtils.glsl"

vec3 camWPos() {
	return vec3(viewMatInv[3]);
}

vec3 getViewDir(vec3 worldPos) {
	if (projMat[3][3] > 0.5) // Check if its orthographic or not
		return vpMatInv[2].xyz;
	else
		return normalize(camWPos() - worldPos);
}