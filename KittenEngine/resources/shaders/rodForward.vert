#version 430 core

#include "kittenCommonVert.glsl"
#include "rod.glsl"

out vec3 col;
out vec3 norm;
out vec3 wPos;
out vec3 mPos;

uniform int numVerts;
uniform float radius;
uniform float density;
uniform float renderByDensity;

void main() {
	Seg seg = segs[gl_InstanceID];
	Vert n0 = verts[seg.i0];
	Vert n1 = verts[seg.i1];

	vec3 axis = n1.pos - n0.pos;
	float l = length(axis);
	
	mat3 basis = orthoBasisY(axis / l);

	vec3 lPos = vPos;

	float m1 = (n1.invMass <= 0) ? 0.5 : 0.5 / n1.invMass;
	float m0 = (n0.invMass <= 0) ? m1 : 0.5 / n0.invMass;

	// Quick radius based on estimating desity
	// Approximately based on m0 + m1 = pi * len / 3 * (r0^2 + r1^2 + r0 * r1) * density
	// Assume r0 * r1 = 0 because r0 must of a sole function of m0 and the same for r1
	float invd = 2 / (3.1415 * density * seg.l);
	float r0 = sqrt(m0 * invd);
	float r1 = sqrt(m1 * invd);
	float rd = mix(r0, r1, vPos.y);

	// Blend with just normal radius
	float r = mix(radius, rd, renderByDensity);
	lPos *= vec3(r, l, r);

	wPos = (modelMat * vec4(basis * lPos + n0.pos, 1)).xyz;
	gl_Position = vpMat * vec4(wPos, 1);

	float k = l / seg.l - 1;
	col = vec3(1, 1, 1);
	col.xy *= 1 - 0.2 * (gl_InstanceID % 4);

	norm = (modelMat * vec4(basis * (vec3(1, 0, 1) * vPos), 0)).xyz;

	// Rotation frame
	mPos = 4 * rotorMatrix(vec4(-seg.q.xyz, seg.q.w)) * basis * (vec3(1, l / radius, 1) * (vPos + vec3(0, -0.5, 0)));
}