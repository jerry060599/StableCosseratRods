#include "../Cosserat.h"

namespace Cosserat {
	void Sim::iterateSegLambda() {
#pragma omp parallel for if(multithreaded) schedule(dynamic, COSSERAT_THREAD_BATCH)
		for (int tid = 0; tid < segs.size(); tid++) {
			Seg seg = segs[tid];

			if (seg.kss < 0) continue;

			const vec3 v0 = verts[seg.i.x].pos;
			const vec3 dx = dxs[seg.i.x];
			const vec3 p1 = verts[seg.i.y].pos;
			const vec3 p1dx = dxs[seg.i.y];
			const vec3 v = (-2 * seg.kss / seg.l) * ((p1 - v0) + (p1dx - dx));

			vec4 b(0);

			// Angle matching energy.
			const int start = segAngleIndices[tid];
			const int end = segAngleIndices[tid + 1];
			for (int s = start; s < end; s++) {
				int angId = segAngleIds[s];
				RestAngle ang = restAngles[angId];

				const int oid = (ang.i.x == tid) ? ang.i.y : ang.i.x;
				const auto oq = segs[oid].q;
				auto qq = oq.inverse() * seg.q;
				if (ang.i.x == tid) qq = qq.inverse();

				float phi = (dot(qq.v, ang.qRest.v) > 0) ? 1 : -1;

				if (ang.i.x == tid) ang.qRest = ang.qRest.inverse();
				b += (ang.kbt * phi) * (oq * ang.qRest).v;
			}

			if (start == end)
				seg.q = Rotor::fromTo(vec3(1, 0, 0), normalize(-v));
			else {
				float lambda = length(v) + length(b);
				seg.q = Rotor(normalize((Rotor(v) * Rotor(b) * Rotor(1)).v + lambda * b));
			}
			segs[tid].q = seg.q;
		}
	}

	void Sim::iterateSegLambdaExact() {
#pragma omp parallel for if(multithreaded) schedule(dynamic, COSSERAT_THREAD_BATCH)
		for (int tid = 0; tid < segs.size(); tid++) {
			Seg seg = segs[tid];

			if (seg.kss < 0) continue;

			const vec3 v0 = verts[seg.i.x].pos;
			const vec3 dx = dxs[seg.i.x];
			const vec3 p1 = verts[seg.i.y].pos;
			const vec3 p1dx = dxs[seg.i.y];
			const vec3 v = (-2 * seg.kss / seg.l) * ((p1 - v0) + (p1dx - dx));

			vec4 b(0);

			// Angle matching energy.
			const int start = segAngleIndices[tid];
			const int end = segAngleIndices[tid + 1];
			for (int s = start; s < end; s++) {
				int angId = segAngleIds[s];
				RestAngle ang = restAngles[angId];

				const int oid = (ang.i.x == tid) ? ang.i.y : ang.i.x;
				const auto oq = segs[oid].q;
				auto qq = oq.inverse() * seg.q;
				if (ang.i.x == tid) qq = qq.inverse();

				float phi = (dot(qq.v, ang.qRest.v) > 0) ? 1 : -1;

				if (ang.i.x == tid) ang.qRest = ang.qRest.inverse();
				b += (ang.kbt * phi) * (oq * ang.qRest).v;
			}

			if (start == end)
				seg.q = Rotor::fromTo(vec3(1, 0, 0), normalize(-v));
			else {
				const float vl = length(v);
				const float bl = length(b);
				float lambda = vl + bl * glm::clamp(lambdaTradGammas[tid], 1e-3f, 1.f);

				lambda = sqrt(length((Rotor(v) * Rotor(b) * Rotor(1)).v + lambda * b) + vl * vl);
				lambdaTradGammas[tid] = (lambda - vl) / bl;
				seg.q = Rotor(normalize((Rotor(v) * Rotor(b) * Rotor(1)).v + lambda * b));
			}
			segs[tid].q = seg.q;
		}
	}
}