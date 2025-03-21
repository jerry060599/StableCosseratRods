#include "../Cosserat.h"
#include "KittenEngine/includes/modules/StopWatch.h"

namespace Cosserat {
	float Sim::getIncrementalPotential() {
		double momentum = 0;
		for (int i = 0; i < verts.size(); i++) {
			Vert v0 = verts[i];
			if (v0.invMass > 0)
				momentum += length2(vels[i] - dxs[i]) / (2 * v0.invMass * meta.h * meta.h);
		}
		return momentum + getPotentialEnergy();
	}

	float Sim::getPotentialEnergy() {
		double energy = 0;

		for (int i = 0; i < segs.size(); i++) {
			Seg seg = segs[i];
			vec3 p0 = verts[seg.i.x].pos + dxs[seg.i.x];
			vec3 p1 = verts[seg.i.y].pos + dxs[seg.i.y];

			vec3 c = (p1 - p0) / seg.l - seg.q * vec3(1, 0, 0);
			energy += 0.5f * abs(seg.kss) * dot(c, c);
		}

		for (int i = 0; i < restAngles.size(); i++) {
			RestAngle ang = restAngles[i];
			auto q0 = segs[ang.i.x].q;
			auto q1 = segs[ang.i.y].q;

			auto qq = q0.inverse() * q1;
			float phi = (dot(qq.v, ang.qRest.v) > 0) ? 1 : -1;
			vec4 c = qq.v - phi * ang.qRest.v;
			energy += 0.5f * ang.kbt * dot(c, c);
		}

		return (float)energy;
	}

	void Sim::checkSegConvergence() {
		// Run it till rest
		auto oldMethod = method;
		method = IntegrationMethod::StableCosseratExact;
		multithreaded = false;
		meta.angleStrainThreshold = meta.segStrainThreshold = INFINITY;
		advance(3);
		meta.damping = 0;
		meta.drag = 0;
		method = oldMethod;

		meta.lastH = meta.h;

		if (method == IntegrationMethod::XPBD) startIterateXPBD();
		else startIterate();

		printf("\nIncremental potential: \n");

		printf("%f\n", getIncrementalPotential());

		for (int itr = 0; itr < 1024; itr++) {
			iterate();
			printf("%f\n", getIncrementalPotential());
		}

		if (method == IntegrationMethod::XPBD) endIterateXPBD();
		else endIterate();
	}

	void Sim::checkSegConvergenceOverTime() {
		vector<Vert> vertsCopy = verts;
		vector<Seg> segsCopy = segs;
		vector<vec3> velsCopy = vels;
		vector<vec3> lastVelsCopy = lastVels;

		FILE* file = fopen("segConvergence.csv", "w");

		float duration = 2.f;
		for (float t = 0; t <= duration; t += maxH) {
			// Save state
			std::copy(verts.begin(), verts.end(), vertsCopy.begin());
			std::copy(segs.begin(), segs.end(), segsCopy.begin());
			std::copy(vels.begin(), vels.end(), velsCopy.begin());
			std::copy(lastVels.begin(), lastVels.end(), lastVelsCopy.begin());

			// Check energy of each method
			meta.numItr = 4;
			float energies[(int)IntegrationMethod::NUM_METHODS];
			for (int k = 0; k < (int)IntegrationMethod::NUM_METHODS; k++) {
				method = (IntegrationMethod)k;
				step(maxH);
				energies[k] = getPotentialEnergy();
				std::copy(vertsCopy.begin(), vertsCopy.end(), verts.begin());
				std::copy(segsCopy.begin(), segsCopy.end(), segs.begin());
				std::copy(velsCopy.begin(), velsCopy.end(), vels.begin());
				std::copy(lastVelsCopy.begin(), lastVelsCopy.end(), lastVels.begin());
			}

			// Run a ridiculous number of iterations for ground truth
			method = IntegrationMethod::StableCosseratExact;
			meta.numItr = 512;
			step(maxH);

			printf("Frame %f\n", t);
			fprintf(file, "%f, %.7f", t, getPotentialEnergy());
			for (int k = 0; k < (int)IntegrationMethod::NUM_METHODS; k++)
				fprintf(file, ", %.7f", energies[k]);
			fprintf(file, "\n");
		}

		fclose(file);
	}

	void Sim::checkExactLambdaConvergence() {
		double mse = 0;
		double mseExact = 0;
		double averageGamma = 0;
		int numFreeSegs = 0;

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

			if (start != end) {
				const float vl = length(v);
				const float bl = length(b);

				float lambda = vl + bl;
				float len = length((Rotor(v) * Rotor(b) * Rotor(1)).v + lambda * b) / (lambda * lambda - vl * vl);

				float gamma = lambdaTradGammas[tid];
				averageGamma += gamma;
				lambda = vl + bl * gamma;
				float lenExact = length((Rotor(v) * Rotor(b) * Rotor(1)).v + lambda * b) / (lambda * lambda - vl * vl);

				mse += pow2(len - 1);
				mseExact += pow2(lenExact - 1);
				numFreeSegs++;
			}
		}

		mse /= numFreeSegs;
		mseExact /= numFreeSegs;
		averageGamma /= numFreeSegs;

		printf("MSE: %e\n", mse);
		printf("MSE Exact: %e\n", mseExact);
		printf("Average Gamma: %f\n", averageGamma);
	}

	void Sim::performSelfCheck() {
		// Check vertex-segment indices
		{
			if (vertSegIndices.size() != verts.size() + 1)
				printf("Vertex segment indices size is %d, expected %d\n", vertSegIndices.size(), verts.size() + 1);
			for (size_t i = 0; i < verts.size(); i++) {
				auto vert = verts[i];
				if (!glm::all(glm::isfinite(vert.pos))) printf("Vertex %d has a non-finite position\n", i);

				int start = vertSegIndices[i];
				int end = vertSegIndices[i + 1];

				if (start > end) printf("Vertex %d has a negative number of segments\n", i);
				if (start < 0 || start > vertSegIds.size()) printf("Vertex %d has an invalid segment start index %d\n", i, start);
				if (end < 0 || end > vertSegIds.size()) printf("Vertex %d has an invalid segment end index %d\n", i, end);
				if (start == end) printf("Vertex %d:(%f %f %f) has no segments\n", i, vert.pos.x, vert.pos.y, vert.pos.z);

				for (int j = start; j < end; j++) {
					int segId = vertSegIds[j];
					ivec2 ind = segs[segId].i;
					if (ind.x != i && ind.y != i)
						printf("Vertex %d references segment %d:(%d %d) which does not connect to this vertex!\n", i, segId, ind.x, ind.y);
				}
			}
			if (vertSegIndices.back() != vertSegIds.size())
				printf("Last vertex segment index is %d, expected %d\n", vertSegIndices.back(), vertSegIds.size());
		}

		// Check segment-angle indices
		{
			if (segAngleIndices.size() != segs.size() + 1)
				printf("Segment angle indices size is %d, expected %d\n", segAngleIndices.size(), segs.size() + 1);
			for (size_t i = 0; i < segs.size(); i++) {
				auto seg = segs[i];
				if (!glm::all(glm::isfinite(seg.q.v))) printf("Segment %d:(%d %d) has a non-finite orientation\n", i, seg.i.x, seg.i.y);
				if (abs(length(seg.q.v) - 1) > 0.001f) printf("Segment %d:(%d %d) has a non-unit quaternion\n", i, seg.i.x, seg.i.y);
				if (!glm::isfinite(seg.l)) printf("Segment %d:(%d %d) has a non-finite rest length\n", i, seg.i.x, seg.i.y);

				int start = segAngleIndices[i];
				int end = segAngleIndices[i + 1];

				if (start > end) printf("Segment %d has a negative number of angles\n", i);
				if (start < 0 || start > segAngleIds.size()) printf("Segment %d has an invalid angle start index %d\n", i, start);
				if (end < 0 || end > segAngleIds.size()) printf("Segment %d has an invalid end angle index %d\n", i, end);

				for (int j = start; j < end; j++) {
					int angId = segAngleIds[j];
					ivec2 ind = restAngles[angId].i;
					if (ind.x != i && ind.y != i)
						printf("Segment %d references angle %d:(%d %d) which does not connect to this segment!\n", i, angId, ind.x, ind.y);
				}
			}
			if (segAngleIndices.back() != segAngleIds.size())
				printf("Last segment angle index is %d, expected %d\n", segAngleIndices.back(), segAngleIds.size());
		}

		// Check segment-vertex indices
		{
			for (size_t i = 0; i < segs.size(); i++) {
				auto seg = segs[i];
				for (int k = 0; k < 2; k++) {
					int v = seg.i[k];
					if (v < 0 || v >= verts.size()) printf("Segment %d:(%d %d) has an invalid vertex index %d\n", i, seg.i.x, seg.i.y, v);

					int start = vertSegIndices[v];
					int end = vertSegIndices[v + 1];
					bool found = false;
					for (int j = start; j < end; j++) {
						if (vertSegIds[j] == i) {
							found = true;
							break;
						}
					}
					if (!found) printf("Segment %d:(%d %d) is not in the list of connected segments for vertex %d\n", i, seg.i.x, seg.i.y, v);
				}
			}
		}

		// Check angle-segment indices
		{
			for (size_t i = 0; i < restAngles.size(); i++) {
				auto ang = restAngles[i];
				if (!glm::all(glm::isfinite(ang.qRest.v))) printf("Angle %d:(%d %d) has a non-finite rest orientation\n", i, ang.i.x, ang.i.y);

				for (int k = 0; k < 2; k++) {
					int s = ang.i[k];
					if (s < 0 || s >= segs.size()) printf("Angle %d:(%d %d) has an invalid segment index %d\n", i, ang.i.x, ang.i.y, s);

					int start = segAngleIndices[s];
					int end = segAngleIndices[s + 1];
					bool found = false;
					for (int j = start; j < end; j++) {
						if (segAngleIds[j] == i) {
							found = true;
							break;
						}
					}
					if (!found) printf("Angle %d:(%d %d) is not in the list of connected angles for segment %d\n", i, ang.i.x, ang.i.y, s);
				}
			}
		}
	}
}