#include "../CosseratExamples.h"

namespace Cosserat {
	Sim* generateCantileverExample() {
		// Example of a simple rod
		Sim* sim = new Cosserat::Sim();

		// Create a rod
		const int numSegs = 10;
		const float segLen = 0.1f;

		// Create vertices
		for (int i = 0; i <= numSegs; i++)
			sim->addVertex(vec3(i * segLen, 0, 0));

		// Link vertices with segments
		for (int i = 0; i < numSegs; i++)
			sim->addSegment(ivec2(i, i + 1));

		// Link segments with rest angles
		for (int i = 0; i < numSegs - 1; i++)
			sim->addRestAngle(ivec2(i, i + 1), 4.f);

		// Pin the first vertex and segment
		sim->pinVertex(0);
		sim->pinSegment(0);

		sim->meshHeightLen = 0.1f * segLen;
		sim->meshRadialSegs = 16;

		sim->meta.drag = 3.0f;

		return sim;
	}

	Sim* generateVBDFailureExample() {
		// Example of a failure case for VBD
		// When elastic material is stretched, the hessian becomes negative definite
		// This is described by the stability condition (N is number of connected bending constraints)
		// length(x_1 - x_0) - l_0 < 2 * N * kbt / (kss * l_0)
		// This comes from smallest eigen value being N * 4 kbt / l - 2 * kss * l * (length(v) / l - 1)
		Sim* sim = new Cosserat::Sim();

		// Create a rod
		const int numSegs = 10;
		const float segLen = 1.f;

		// Create vertices
		srand(0);
		for (int i = 0; i <= numSegs; i++)
			sim->addVertex(vec3(i * segLen, 0.5f * segLen * (i % 2), 0.1f * segLen * (i % 3)));

		// Link vertices with segments
		for (int i = 0; i < numSegs; i++)
			sim->addSegment(ivec2(i, i + 1), 1e4, 0.2f);

		// Link segments with rest angles
		for (int i = 0; i < numSegs - 1; i++)
			sim->addRestAngle(ivec2(i, i + 1), 1e2);

		// Pin the first vertex and segment
		sim->pinVertex(0);
		sim->pinVertex(10);
		sim->pinSegment(0);
		sim->maxH = 1e-3;
		sim->meta.damping = 1e-4f;
		sim->verts[10].pos.x += 2 * segLen;

		sim->renderRadius = 0.3 * segLen;

		sim->meshHeightLen = 0.1f * segLen;
		sim->meshRadialSegs = 16;

		return sim;
	}
}