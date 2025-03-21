#include "../CosseratExamples.h"

namespace Cosserat {
	constexpr vec3 teapotPoints[] = {
		vec3(0, 28, 0), vec3(2, 21, 0), vec3(7, 19, 0), vec3(12, 18, 0), vec3(20, 18, 0),
		vec3(29, 18, 0), vec3(31, 13, 0), vec3(39, 13, 0), vec3(48, 11, 0), vec3(55, 8, 0),
		vec3(52, 3, 0), vec3(59, 1, 0), vec3(66, 3, 0), vec3(63, 9, 0), vec3(72, 11, 0),
		vec3(80, 13, 0), vec3(88, 14, 0), vec3(93, 22, 0), vec3(95, 28, 0), vec3(98, 34, 0),
		vec3(104, 31, 0), vec3(108, 26, 0), vec3(109, 20, 0), vec3(115, 14, 0), vec3(126, 14, 0),
		vec3(119, 19, 0), vec3(115, 25, 0), vec3(113, 32, 0), vec3(111, 39, 0), vec3(106, 45, 0),
		vec3(98, 50, 0), vec3(93, 56, 0), vec3(86, 61, 0), vec3(77, 62, 0), vec3(67, 62, 0),
		vec3(55, 62, 0), vec3(43, 62, 0), vec3(31, 61, 0), vec3(25, 56, 0), vec3(21, 50, 0),
		vec3(12, 46, 0), vec3(6, 41, 0), vec3(3, 36, 0)
	};

	Sim* generateSlinkyExample() {
		// Example of a simple rod
		Sim* sim = new Cosserat::Sim();

		// Geometry settings
		const int numSegs = 32 * COSSERAT_THREAD_BATCH;
		const float pitch = 0.1f;	// Spacing between coils
		const float scale = 5e-3f;

		// Material settings
		const float density = 0.15;
		const float kss = 1e3;
		const float kbt = 1e3;

		constexpr int numPoints = sizeof(teapotPoints) / sizeof(vec3);
		float totalLen = 0;
		for (int i = 0; i < numPoints; i++)
			totalLen += glm::length(teapotPoints[i] - teapotPoints[(i + 1) % numPoints]);

		float ts[numPoints] = { 0 };

		float runningLen = 0;
		for (int i = 0; i < numPoints; i++) {
			ts[i] = runningLen / totalLen;
			runningLen += glm::length(teapotPoints[i] - teapotPoints[(i + 1) % numPoints]);
		}

		for (int i = 0; i <= numSegs; i++) {
			float y = (i / numPoints + ts[i % numPoints]) * pitch;
			vec3 point = teapotPoints[i % numPoints] * scale;
			point.y = -point.y;
			sim->addVertex(vec3(point.x, point.y, y));
		}

		for (int i = 0; i < numSegs; i++)
			sim->addSegment(ivec2(i, i + 1), kss, density);

		for (int i = 0; i < numSegs - 1; i++)
			sim->addRestAngle(ivec2(i, i + 1), kbt);

		// Pin the first vertex and segment
		sim->pinVertex(0);
		sim->pinSegment(0);

		sim->pinVertex(numSegs);
		sim->pinSegment(numSegs - 1);

		sim->maxH = 3e-4;
		sim->meta.numItr = 8;
		sim->meta.damping = 5e-8;
		sim->meta.drag = 0;

		sim->renderRadius = scale * 0.4f * totalLen / numPoints;

		return sim;
	}
}