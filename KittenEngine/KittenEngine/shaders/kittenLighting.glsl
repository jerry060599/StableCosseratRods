
layout (binding = 3, std140) uniform lightUBO {
	vec4 light_col;

	vec3 light_dir;
	float light_radius;

	float light_bias;
	int light_hasShadow;
	float light_spread;
	float light_focus;

	vec3 light_pos;
	int light_type;

	mat4 light_shadowProj;
};

layout(binding = 6) uniform sampler2DShadow tex_shadow;

vec3 getLightDir(vec3 worldPos) {
	if (light_type == 3)
		return light_dir;
	return normalize(worldPos - light_pos);
}

vec3 getLightCol() {
	return light_col.xyz * light_col.w;
}

float getLightSpread(vec3 worldPos) {
	float d = dot(getLightDir(worldPos), light_dir);
	d = clamp((d - light_spread) / (1 - light_spread), 0, 1);
	d = pow(d, light_focus);
	return d;
}

float getLightShadow(vec3 worldPos) {
	vec4 p = mat4(
				0.5, 0.0, 0.0, 0.0,
				0.0, 0.5, 0.0, 0.0,
				0.0, 0.0, 0.5, 0.0,
				0.5, 0.5, 0.5, 1.0
			) * (light_shadowProj * vec4(worldPos, 1));
	p.z -= light_bias;
	return textureProj(tex_shadow, p);
}

vec3 getLightCol(vec3 worldPos){
	vec3 c = getLightCol();
	float d2 = length2(worldPos - light_pos);
	float h = d2 + pow2(light_radius);
	float distSqrInv = 2 / (h + sqrt(h * d2));

	switch(light_type) {
	case 0: // Ambient
		return c;
	case 1: // Point
		return c * distSqrInv;
	case 2: // Spot
		return c * distSqrInv * getLightSpread(worldPos) * getLightShadow(worldPos);
	case 3: // Dir
		return c * getLightShadow(worldPos);
	}
}