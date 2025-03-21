#include "../Cosserat.h"
#include <iostream>

namespace Cosserat {
	void Sim::exportMesh(const char* path) {
		FILE* file = fopen(path, "w");
		if (!file)
			throw std::runtime_error("Failed to open file for writing");

		const float densityScale = 3 / (3.1415 * renderDensity);
		const float uvScale = 1.0f;

		int numVerts = 0;
		for (size_t i = 0; i < segs.size(); i++) {
			auto seg = segs[i];
			auto v0 = verts[seg.i.x];
			auto v1 = verts[seg.i.y];

			int angStart = segAngleIndices[i];
			int angEnd = segAngleIndices[i + 1];

			// Radius estimation
			vec2 masses((v0.invMass > 0) ? 1 / v0.invMass : 1, (v1.invMass > 0) ? 1 / v1.invMass : 1);
			masses[0] *= 0.5f;
			masses[1] *= 0.5f;

			// Find any neighboring segments
			ivec2 neighbors(-1);
			Vert nVerts[2]{ v0, v1 };
			Seg nsegs[2]{ seg, seg };
			for (int k = 0; k < 2; k++) {
				int start = vertSegIndices[seg.i[k]];
				int end = vertSegIndices[seg.i[k] + 1];

				// If there are multiple connections, prioritize the most "similar" one
				float minDiff = std::numeric_limits<float>::infinity();

				for (int j = start; j < end; j++) {
					int segId = vertSegIds[j];
					if (segId == i) continue;

					auto os = segs[segId];

					// If the seg length is different
					float score = abs(seg.l - os.l);
					// If the angle is different
					score += 1 - ((os.i[1 - k] == seg.i[k]) ? 1 : -1) * dot(seg.q * vec3(1, 0, 0), os.q * vec3(1, 0, 0));

					if (abs(abs(seg.kss / seg.l) - abs(os.kss / os.l)) > 10) score += 0.3;

					// Grab the min score
					if (score < minDiff) {
						minDiff = score;
						neighbors[k] = segId;
					}
				}

				if (end - start > 2)
					minDiff += meshJunctionPenalty;

				// If the most similar segment is too different, ignore it
				if (end - start > 1 + meshExcludeInlinePhong && minDiff > meshPhongThreshold) {
					neighbors[k] = -1;
					masses[k] = masses[1 - k];
				}

				if (neighbors[k] != -1) {
					auto os = segs[neighbors[k]];
					nVerts[k] = (os.i.x == seg.i[k]) ? verts[os.i.y] : verts[os.i.x];

					// If the segments go in opposite directions, flip the tangent frame
					// This happens if k=0 (lower vertex of this segment) and os.i.x == seg.i.x (butt against butt) OR
					// if k=1 (upper vertex of this segment) and os.i.y == seg.i.y (tip against tip)
					nsegs[k] = os;
					if ((k == 0 && os.i.x == seg.i.x) || (k == 1 && os.i.y == seg.i.y))
						nsegs[k].q = Rotor(0, 1, 0, 0) * nsegs[k].q; // Rotate 180 degrees around y

					if (dot(nsegs[k].q.v, seg.q.v) < 0) nsegs[k].q.v = -nsegs[k].q.v;
				}
			}

			if (v0.invMass == 0) masses[0] = masses[1];
			if (v1.invMass == 0) masses[1] = masses[0];

			float r0 = sqrt(masses[0] * densityScale / (0.5f * (seg.l + nsegs[0].l)));
			float r1 = sqrt(masses[1] * densityScale / (0.5f * (seg.l + nsegs[1].l)));
			r0 = mix(renderRadius, r0, renderByDensity);
			r1 = mix(renderRadius, r1, renderByDensity);

			Rotor lastFrame;
			int numHeightSegments = (int)ceil(seg.l / meshHeightLen);
			if (numHeightSegments < 1) numHeightSegments = 1;
			for (int hi = 0; hi <= numHeightSegments; hi++) {
				// CMR3 interpolation to get position
				float t = hi / (float)numHeightSegments;
				vec3 pos = Kit::cmrSpline(nVerts[0].pos, v0.pos, v1.pos, nVerts[1].pos, t);
				Rotor frame = seg.q;
				if (t < 0.5f) frame = normalize(mix(nsegs[0].q.v, seg.q.v, t + 0.5f));
				else frame = normalize(mix(seg.q.v, nsegs[1].q.v, t - 0.5f));

				float radius = mix(r0, r1, t);

				// Write out a circle in obj format
				for (int ri = 0; ri < meshRadialSegs; ri++) {
					float theta = 2 * glm::pi<float>() * ri / (float)meshRadialSegs;
					vec3 lPos = vec3(0, cos(theta), sin(theta));
					vec3 wPos = pos + frame * (radius * lPos);
					fprintf(file, "v %.8f %.8f %.8f %.3f %.3f %.3f\n", wPos.x, wPos.y, wPos.z, radius * 0.1f, radius * 0.5f, radius);
					fprintf(file, "vt %.3f %.3f %.3f\n", uvScale * (1 - 2 * abs(ri / (float)meshRadialSegs - 0.5f)), t + i, lPos.z);
					numVerts++;
				}

				// Write out faces
				if (hi > 0) {
					// Because the tubes can twist, we want to shift the triangle indices to match the previous frame
					float sinAngle = 2 * asinf((lastFrame.inverse() * frame).x);
					int twist = meshRadialSegs - (int)floor((float)meshRadialSegs * sinAngle / (2 * glm::pi<float>()));

					for (int ri = 0; ri < meshRadialSegs; ri++) {
						int i0 = numVerts - 2 * meshRadialSegs + ri + 1;
						int i1 = numVerts - 2 * meshRadialSegs + (ri + 1) % meshRadialSegs + 1;
						int i2 = numVerts - meshRadialSegs + (ri + twist) % meshRadialSegs + 1;
						int i3 = numVerts - meshRadialSegs + (ri + twist + 1) % meshRadialSegs + 1;

						fprintf(file, "f %d/%d %d/%d %d/%d\n", i0, i0, i1, i1, i2, i2);
						fprintf(file, "f %d/%d %d/%d %d/%d\n", i3, i3, i2, i2, i1, i1);
					}
				}
				lastFrame = frame;
			}

			// Write out endcaps
			for (int k = 0; k < 2; k++)
				if (neighbors[k] < 0) {
					float radius = (k == 0) ? r0 : r1;
					// Write out vertices
					for (int ri = 0; ri < meshRadialSegs; ri++) {
						float theta = 2 * glm::pi<float>() * ri / (float)meshRadialSegs;
						vec3 lPos = vec3(0, cos(theta), sin(theta));
						vec3 wPos = seg.q * (radius * lPos);
						wPos += (k == 0) ? v0.pos : v1.pos;
						fprintf(file, "v %.8f %.8f %.8f %.3f %.3f %.3f\n", wPos.x, wPos.y, wPos.z, radius * 0.1f, radius * 0.5f, radius);
						fprintf(file, "vt %.3f %.3f %.3f\n", lPos.y, lPos.z, i + (float)k);
						numVerts++;
					}

					// Write out triangles
					for (int ri = 0; ri < meshRadialSegs; ri++) {
						int i0 = numVerts - meshRadialSegs + 1;
						int i1 = numVerts - meshRadialSegs + ri + 1;
						int i2 = numVerts - meshRadialSegs + (ri + 1) % meshRadialSegs + 1;
						if (k) std::swap(i1, i2);

						fprintf(file, "f %d/%d %d/%d %d/%d\n", i0, i0, i1, i1, i2, i2);
					}
				}
		}

		fclose(file);
	}
}
