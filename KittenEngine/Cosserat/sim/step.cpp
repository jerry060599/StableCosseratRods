#include "../Cosserat.h"

namespace Cosserat {
	void Sim::step(float h) {
		float temp = maxH;
		maxH = h;
		advance(h);
		maxH = temp;
	}

	float Sim::advance(float h) {
		if (h <= 0) return 0;

		int steps = glm::max(1, (int)ceil(h / maxH));
		meta.lastH = meta.h;
		meta.h = h / steps;

		for (int s = 0; s < steps; s++, stepCounter++) {
			if (method == IntegrationMethod::XPBD) startIterateXPBD();
			else startIterate();

			for (int itr = 0; itr < meta.numItr; itr++)
				iterate();

			if (method == IntegrationMethod::XPBD) endIterateXPBD();
			else endIterate();
			meta.time += meta.h;
		}

		// Perform breaking
		destroyByStrainMagnitude(meta.segStrainThreshold, meta.angleStrainThreshold);
		return h;
	}
}