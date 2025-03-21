#include "../Cosserat.h"
#include "KittenEngine/includes/modules/SymMat.h"

namespace Cosserat {
	void Sim::startIterate() {
		const float h = meta.h;
		const float invLastH = 1 / meta.lastH;
		const vec3 g = meta.gravity;

#pragma omp parallel for if(multithreaded) schedule(static, COSSERAT_THREAD_BATCH_STATIC)
		for (int tid = 0; tid < verts.size(); tid++) {
			auto vel = vels[tid];
			vec3 dx = vel * h;
			vec3 y = dx;

			float invMass = verts[tid].invMass;
			if (invMass > 0) {
				vec3 aExt = g;
				if (externalForce)
					aExt += externalForce(this, verts[tid], meta.time) * invMass;

				y += (h * h) * aExt;

				float ae2 = length2(aExt);
				if (ae2 > 0) {
					vec3 a = (vel - lastVels[tid]) * invLastH;
					float s = glm::clamp(dot(a, aExt) / ae2, 0.f, 1.f);
					dx += (h * h * s) * aExt;
				}
			}

			lastVels[tid] = vel;
			vels[tid] = y;
			dxs[tid] = dx;
		}
	}

	void Sim::endIterate() {
		const float h = meta.h;
		const float invH = 1 / h;

#pragma omp parallel for if(multithreaded) schedule(static, COSSERAT_THREAD_BATCH_STATIC)
		for (int tid = 0; tid < verts.size(); tid++) {
			vec3 dx = dxs[tid];
			if (verts[tid].invMass > 0)
				vels[tid] = dx * (invH * (1 - meta.drag * h));
			else
				vels[tid] = dx * invH;
			verts[tid].pos += dx;
		}
	}

	void Sim::iterate() {
		switch (method) {
		case IntegrationMethod::StableCosserat:
			iterateVertVBD();
			iterateSegLambda();
			break;
		case IntegrationMethod::StableCosseratExact:
			iterateVertVBD();
			iterateSegLambdaExact();
			break;
		case IntegrationMethod::VBD:
			iterateVertVBD();
			iterateSegVBD();
			break;
		case IntegrationMethod::LCVBD:
			iterateVertVBD();
			iterateSegLCVBD();
			break;
		case IntegrationMethod::PSDPLCVBD:
			iterateVertVBD();
			iterateSegPSDPLCVBD();
			break;
		case IntegrationMethod::XPBD:
			iterateVertXPBD();
			iterateSegXPBD();
			break;
		}
	}
}