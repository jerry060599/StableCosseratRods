#include "../Cosserat.h"

namespace Cosserat {
	void Sim::startIterateXPBD() {
		const float h = meta.h;
		const vec3 g = meta.gravity;

		for (size_t i = 0; i < verts.size(); i++) {
			vec3 dx = vels[i] * h;
			if (verts[i].invMass > 0)
				dx += (h * h) * g;
			dxs[i] = dx;
		}

		for (size_t i = 0; i < xpbdSegLambdas.size(); i++) {
			xpbdSegLambdas[i] = vec3(0);
			Rotor q = segs[i].q;
			xpbdLastQs[i] = q;

			q = normalize(q.v + (0.5f * h) * (Rotor(xpbdSegWs[i], 0) * q).v);
			segs[i].q = q;
		}
		for (size_t i = 0; i < xpbdAngleLambdas.size(); i++)
			xpbdAngleLambdas[i] = vec4(0);
	}

	void Sim::endIterateXPBD() {
		const float h = meta.h;

		for (size_t i = 0; i < verts.size(); i++) {
			vec3 dx = dxs[i];
			verts[i].pos += dx;
			if (verts[i].invMass > 0)
				vels[i] = dx * ((1 / h) * (1 - meta.drag * h));
		}

		for (size_t i = 0; i < xpbdSegLambdas.size(); i++) {
			Rotor q = segs[i].q;
			Rotor lastQ = xpbdLastQs[i];
			xpbdSegWs[i] = (q * lastQ.inverse()).q * (2 / h);
		}
	}

	inline float getInvMoment(Seg& seg) {
		if (seg.kss < 0) return 0;		// If orientation pinned
		return 12.f / Kit::pow3(seg.l);
	}

	void Sim::iterateVertXPBD() {
		const float h = meta.h;

		for (int tid = 0; tid < segs.size(); tid++) {
			auto seg = segs[tid];
			float invl = 1 / seg.l;

			// Calculate constraint
			auto v0 = verts[seg.i.x];
			auto v0dx = dxs[seg.i.x];
			auto v1 = verts[seg.i.y];
			auto v1dx = dxs[seg.i.y];
			vec3 c = ((v1.pos - v0.pos) + (v1dx - v0dx)) * invl - seg.q * vec3(1, 0, 0);

			// Update lambda
			vec3 lambda = xpbdSegLambdas[tid];
			const float ALPHA_STRETCH = 1 / (glm::abs(seg.kss) * h * h);
			vec3 dl = -c - ALPHA_STRETCH * lambda;

			float invMoment = getInvMoment(seg);
			float lSqr = seg.l * seg.l;
			dl *= lSqr / (v0.invMass + v1.invMass + (4 * invMoment + ALPHA_STRETCH) * lSqr);

			// Update lambda
			lambda += dl;
			xpbdSegLambdas[tid] = lambda;

			// Update position
			vec3 dx = dl / seg.l;
			dxs[seg.i.x] += -v0.invMass * dx;
			dxs[seg.i.y] += v1.invMass * dx;

			seg.q += invMoment * (Rotor(dl, 0) * seg.q * Rotor(2)).v;
			seg.q = normalize(seg.q);
			segs[tid].q = seg.q;
		}
	}

	void Sim::iterateSegXPBD() {
		const float h = meta.h;

		for (int tid = 0; tid < restAngles.size(); tid++) {
			auto angle = restAngles[tid];

			auto seg0 = segs[angle.i.x];
			auto seg1 = segs[angle.i.y];

			Rotor qq = seg0.q.inverse() * seg1.q;
			float phi = (dot(qq.v, angle.qRest.v) > 0) ? 1 : -1;
			vec4 c = qq.v - phi * angle.qRest.v;

			// Compute delta lambda
			const float ALPHA_BEND = 1 / (angle.kbt * h * h);
			vec4 lambda = xpbdAngleLambdas[tid];

			vec4 dl = -c - ALPHA_BEND * lambda;
			float invI0 = getInvMoment(seg0);
			float invI1 = getInvMoment(seg1);
			dl *= 1 / (invI0 + invI1 + ALPHA_BEND);

			// Lambda update
			lambda += dl;
			xpbdAngleLambdas[tid] = lambda;

			// Rotation update
			segs[angle.i.x].q = normalize(seg0.q.v + invI0 * (seg1.q * ((Rotor)dl).inverse()).v);
			segs[angle.i.y].q = normalize(seg1.q.v + invI1 * (seg0.q * (Rotor)dl).v);
		}
	}
}