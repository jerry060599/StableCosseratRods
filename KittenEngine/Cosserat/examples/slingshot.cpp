#include "../CosseratExamples.h"
#include "fastPRGN.h"

namespace Cosserat {
	inline float randf(fastPRNG::fastXS64& rng) {
		return glm::clamp((uint32_t)rng.xoroshiro128p() / (float)UINT32_MAX, 0.f, 1.f);
	}

	Sim* generateSlingshotExample() {
		Sim* sim = new Cosserat::Sim();

		// Geometry
		const int numHandleSegs = 8;
		const float forkArc = 0.4 * glm::pi<float>();
		const float forkRadius = 1.f;
		const float segLen = 0.3f;	// 1 cm
		const float cableSegLen = 0.15f;

		// Material
		const float density = 3;
		const float densityRubber = 0.2;
		const float kssFork = 1e6;
		const float kbtFork = 1e5;
		const float kssRubber = 4e2;
		const float kbtRubber = 0.1;

		// Create handle
		for (int i = 0; i <= numHandleSegs; i++)
			sim->addVertex(vec3(0, i * segLen, 0));
		for (int i = 0; i < numHandleSegs; i++)
			sim->addSegment(ivec2(i, i + 1), kssFork, density);
		for (int i = 0; i < numHandleSegs - 1; i++)
			sim->addRestAngle(ivec2(i, i + 1), kbtFork);

		// Create fork
		int numForkSegs = (int)ceil(forkArc * forkRadius / segLen);

		// One fork
		int vOffset = sim->verts.size();
		int sOffset = sim->segs.size();
		for (int i = 1; i <= numForkSegs; i++) {
			float angle = forkArc * i / (float)numForkSegs;
			sim->addVertex(vec3(forkRadius * sin(angle), numHandleSegs * segLen + forkRadius * (1 - cos(angle)), 0));
		}
		sim->addSegment(ivec2(numHandleSegs, vOffset), kssFork, density);
		sim->addRestAngle(ivec2(numHandleSegs - 1, sim->segs.size() - 1), kbtFork);
		for (int i = 1; i < numForkSegs; i++)
			sim->addSegment(ivec2(i - 1, i) + ivec2(vOffset), kssFork, density);
		for (int i = 0; i < numForkSegs - 1; i++)
			sim->addRestAngle(ivec2(i, i + 1) + ivec2(sOffset), kbtFork);

		// Two fork
		vOffset = sim->verts.size();
		sOffset = sim->segs.size();
		for (int i = 1; i <= numForkSegs; i++) {
			float angle = forkArc * i / (float)numForkSegs;
			sim->addVertex(vec3(-forkRadius * sin(angle), numHandleSegs * segLen + forkRadius * (1 - cos(angle)), 0));
		}
		sim->addSegment(ivec2(numHandleSegs, vOffset), kssFork, density);
		sim->addRestAngle(ivec2(numHandleSegs - 1, sim->segs.size() - 1), kbtFork);
		for (int i = 1; i < numForkSegs; i++)
			sim->addSegment(ivec2(i - 1, i) + ivec2(vOffset), kssFork, density);
		for (int i = 0; i < numForkSegs - 1; i++)
			sim->addRestAngle(ivec2(i, i + 1) + ivec2(sOffset), kbtFork);

		printf("Handle indices end at %d\n", sim->verts.size());

		// Add cables
		int centerVertex = sim->verts.size();
		sim->addVertex(vec3(0, numHandleSegs * segLen + forkRadius, 0));

		using namespace fastPRNG;
		fastXS64 prng(1);

		auto connectCable = [&](int i, int j) {
			vec3 a = sim->verts[i].pos;
			vec3 b = sim->verts[j].pos;
			float len = length(b - a);
			int numSegs = (int)ceil(len / cableSegLen);
			numSegs = glm::max(numSegs, 3);

			const int cableVOffset = sim->verts.size();
			const int cableSOffset = sim->segs.size();
			for (int i = 1; i < numSegs; i++)
				sim->addVertex(a + (i / (float)numSegs) * (b - a) + 0.2f * cableSegLen * (vec3(randf(prng), randf(prng), randf(prng)) * 2.f - 1.f));
			sim->addSegment(ivec2(i, cableVOffset), kssRubber, densityRubber);
			for (int i = 0; i < numSegs - 2; i++)
				sim->addSegment(ivec2(cableVOffset + i, cableVOffset + i + 1), kssRubber, densityRubber);
			sim->addSegment(ivec2(cableVOffset + numSegs - 2, j), kssRubber, densityRubber);
			for (int i = 0; i < numSegs - 1; i++)
				sim->addRestAngle(ivec2(i, i + 1) + ivec2(cableSOffset), kbtRubber);
			};

		ivec4 lastSegs;
		ivec4 nearCenterVertices;
		lastSegs[0] = sim->segs.size();
		nearCenterVertices[0] = sim->verts.size();
		connectCable(centerVertex, numHandleSegs + numForkSegs);
		lastSegs[1] = sim->segs.size();
		nearCenterVertices[1] = sim->verts.size();
		connectCable(centerVertex, numHandleSegs + numForkSegs - 1);
		lastSegs[2] = sim->segs.size();
		nearCenterVertices[2] = sim->verts.size();
		connectCable(centerVertex, numHandleSegs + 2 * numForkSegs);
		lastSegs[3] = sim->segs.size();
		nearCenterVertices[3] = sim->verts.size();
		connectCable(centerVertex, numHandleSegs + 2 * numForkSegs - 1);

		// Add rest angle between pairs of cables
		sim->addRestAngle(ivec2(lastSegs[0], lastSegs[3]), kbtRubber);
		sim->restAngles.back().qRest = Rotor(0, 1, 0, 0);
		sim->addRestAngle(ivec2(lastSegs[1], lastSegs[2]), kbtRubber);
		sim->restAngles.back().qRest = Rotor(0, 1, 0, 0);

		// Add rest angle between cables and handle
		sim->addRestAngle(ivec2(lastSegs[1] - 1, numHandleSegs + numForkSegs - 1), kbtRubber);
		sim->addRestAngle(ivec2(lastSegs[2] - 1, numHandleSegs + numForkSegs - 1), kbtRubber);
		sim->addRestAngle(ivec2(lastSegs[3] - 1, numHandleSegs + 2 * numForkSegs - 1), kbtRubber);
		sim->addRestAngle(ivec2(sim->segs.size() - 1, numHandleSegs + 2 * numForkSegs - 1), kbtRubber);

		// Tighten the rubber a little from rest
		for (int i = lastSegs[0]; i < sim->segs.size(); i++)
			sim->segs[i].l *= 0.94f;

		// Pin the first vertex and segment
		sim->pinVertex(0);
		sim->pinSegment(0);

		printf("Slingshot center vertex is vertex %d\n", centerVertex);
		printf("Slingshot near center vertices %d, %d, %d, %d\n",
			nearCenterVertices[0], nearCenterVertices[1], nearCenterVertices[2], nearCenterVertices[3]);
		
		sim->verts[centerVertex].pos -= vec3(0, 1.3f, 5.f);
		for (int k = 0; k < 4; k++) {
			sim->verts[nearCenterVertices[k]].pos -= 0.98f * vec3(0, 1.3f, 5.f);
			sim->pinVertex(nearCenterVertices[k]);
		}
		sim->pinVertex(centerVertex);

		sim->maxH = 1e-3;
		sim->meta.numItr = 4;
		sim->meta.drag = 1.f;
		sim->meta.damping = 5e-3;
		sim->meta.useProportionalDamping = false;

		sim->renderByDensity = 1;
		sim->renderDensity = 90;

		sim->meshJunctionPenalty = 0.2;
		sim->meshRadialSegs = 20;
		sim->meshPhongThreshold = 0.5f;

		return sim;
	}
}