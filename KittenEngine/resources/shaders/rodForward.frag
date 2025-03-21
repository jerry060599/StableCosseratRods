#version 430 core

#include "kittenCommonFrag.glsl"
#include "kittenLighting.glsl"

in vec3 col;
in vec3 norm;
in vec3 wPos;
in vec3 mPos;

out vec4 fragColor;

void main() {
	vec3 n = norm / length(norm);
	vec3 lDir = getLightDir(wPos);
	vec3 halfVec = 0.5 * (-lDir + getViewDir(wPos));

	float d = max(dot(-lDir, n), 0);

	vec4 mCol = matColor;
	mCol.yz *= -1 < mPos.x && mPos.x < 1 ? 0.5 : 1;
	mCol.xz *= -1 < mPos.y && mPos.y < 1 ? 0.5 : 1;
	mCol.xy *= -1 < mPos.z && mPos.z < 1 ? 0.5 : 1;

	fragColor = vec4(getLightCol(wPos) * max(d, 0.1) * mCol.xyz * col, 0);
}