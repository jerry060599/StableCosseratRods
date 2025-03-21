
const float pi = 3.14159265359;
const float e = 2.71828182846;

float cross2d(vec2 a, vec2 b) {
	return a.x * b.y - a.y * b.x;
}

float length2(vec3 v) {
	return dot(v, v);
}

float length2(vec2 v) {
	return dot(v, v);
}

float pow2(float v) {
	return v * v;
}

float pow3(float v) {
	return v * v * v;
}

float pow4(float v) {
	v *= v;
	return v * v;
}

vec3 reflect(vec3 v, vec3 norm) {
	return v - 2 * dot(norm, v) * norm;
}

vec3 hue2rgb(float h) {
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(h + K.xyz) * 6.0 - K.www);
    return clamp(p - K.xxx, 0.0, 1.0);
}


// http://psgraphics.blogspot.com/2014/11/making-orthonormal-basis-from-unit.html
mat3 orthoBasisX(vec3 n) {
	mat3 basis;
	basis[0] = n;
	if (n.z >= n.y) {
		const float a = 1.0f / (1.0f + n.z);
		const float b = -n.x * n.y * a;
		basis[1] = vec3(b, 1.0f - n.y * n.y * a, -n.y);
		basis[2] = -vec3(1.0f - n.x * n.x * a, b, -n.x);
	}
	else {
		const float a = 1.0f / (1.0f + n.y);
		const float b = -n.x * n.z * a;
		basis[1] = vec3(1.0f - n.x * n.x * a, -n.x, b);
		basis[2] = -vec3(b, -n.z, 1.0f - n.z * n.z * a);
	}
	return basis;
}

mat3 orthoBasisY(vec3 n) {
	mat3 basis;
	basis[1] = n;
	if (n.z >= n.y) {
		const float a = 1.0f / (1.0f + n.z);
		const float b = -n.x * n.y * a;
		basis[0] = vec3(b, 1.0f - n.y * n.y * a, -n.y);
		basis[2] = vec3(1.0f - n.x * n.x * a, b, -n.x);
	}
	else {
		const float a = 1.0f / (1.0f + n.y);
		const float b = -n.x * n.z * a;
		basis[0] = vec3(1.0f - n.x * n.x * a, -n.x, b);
		basis[2] = vec3(b, -n.z, 1.0f - n.z * n.z * a);
	}
	return basis;
}

mat3 orthoBasisZ(vec3 n) {
	mat3 basis;
	basis[2] = n;
	if (n.z >= n.y) {
		const float a = 1.0f / (1.0f + n.z);
		const float b = -n.x * n.y * a;
		basis[0] = vec3(1.0f - n.x * n.x * a, b, -n.x);
		basis[1] = vec3(b, 1.0f - n.y * n.y * a, -n.y);
	}
	else {
		const float a = 1.0f / (1.0f + n.y);
		const float b = -n.x * n.z * a;
		basis[0] = vec3(b, -n.z, 1.0f - n.z * n.z * a);
		basis[1] = vec3(1.0f - n.x * n.x * a, -n.x, b);
	}
	return basis;
}

mat3 abT(vec3 a, vec3 b) {
	return mat3(
		b.x * a.x, b.y * a.x, b.z * a.x,
		b.x * a.y, b.y * a.y, b.z * a.y,
		b.x * a.z, b.y * a.z, b.z * a.z
	);
}

mat3 crossMatrix(vec3 v) {
	return mat3(
		0, v.z, -v.y,
		-v.z, 0, v.x,
		v.y, -v.x, 0
		);
}

vec3 applyRotor(vec4 r, vec3 v) {
	// Calculate v * ab
	vec3 a = r.w * v + cross(r.xyz, v);	// The vector
	float c = dot(v, r.xyz);			// The trivector

	// Calculate (r.w - r.xyz) * (a + c). Ignoring the scaler-trivector parts
	return r.w * a			// The scaler-vector product
		+ cross(r.xyz, a)	// The bivector-vector product
		+ c * r.xyz;		// The bivector-trivector product
}

mat3 rotorMatrix(vec4 r) {
	mat3 cm = crossMatrix(r.xyz);
	return abT(r.xyz, r.xyz) + mat3(r.w * r.w) + 2 * r.w * cm + cm * cm;
}

// Cubic Catmull-Rom Spline
vec3 cmrSpline(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float t) {
	vec3 a0 = mix(p0, p1, t + 1);
	vec3 a1 = mix(p1, p2, t + 0);
	vec3 a2 = mix(p2, p3, t - 1);

	vec3 b0 = mix(a0, a1, (t + 1) * 0.5f);
	vec3 b1 = mix(a1, a2, (t + 0) * 0.5f);

	return mix(b0, b1, t);
}

// Tangent of the Cubic Catmull-Rom Spline
vec3 cmrSplineTangent(vec3 p0, vec3 p1, vec3 p2, vec3 p3, float t) {
	return 0.5f * (
		(-(t - 1) * (3 * t - 1)) * p0 +
		(t * (9 * t - 10)) * p1 +
		-(t - 1) * (9 * t + 1) * p2 +
		t * (3 * t - 2) * p3
		);
}
