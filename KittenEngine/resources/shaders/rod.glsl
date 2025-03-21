
struct Vert {
	vec3 pos;				// Node position
	float invMass;			// Inverse nodal mass
};

struct Seg {
	vec4 q;
	int i0;
	int i1;
	float kss;
	float l;
};

// In in all shader types
layout(binding = 5, std430) buffer vertBlock {
	Vert verts[];
};

// In in all shader types
layout(binding = 6, std430) buffer segBlock {
	Seg segs[];
};

