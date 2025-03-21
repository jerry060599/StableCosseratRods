#include "../Cosserat.h"

namespace Cosserat {
	void Sim::setWindFieldExternalForce(std::function<vec3(vec3, float)> windField, float dragCoeff) {
		externalForce = [=](Cosserat::Sim*, Cosserat::Vert& vert, float time) -> vec3 {
			int vid = &vert - &verts[0];
			vec3 pos = vert.pos;
			vec3 wind = windField(pos, time);
			vec3 dv = wind - vels[vid];

			float dvl = length(dv);
			if (dvl == 0) return vec3(0);
			float idvl2 = 1 / (dvl * dvl);

			const float estRadius = sqrt(1 / vert.invMass);

			vec3 force(0);
			int start = vertSegIndices[vid];
			int end = vertSegIndices[vid + 1];
			for (int i = start; i < end; i++) {
				int segId = vertSegIds[i];
				auto seg = segs[segId];
				auto other = verts[seg.i.x == vid ? seg.i.y : seg.i.x];

				vec3 dir = normalize(other.pos - pos);								// Wind direction
				float sinSqr = glm::max(1 - idvl2 * Kit::pow2(dot(dir, dv)), 0.f);	// Angle of segment against wind
				float referenceArea = estRadius * seg.l * glm::sqrt(sinSqr);		// Area of segment against wind

				force += (dragCoeff * dvl * referenceArea) * dv;					// Wind force is squared velocity
			}

			return force;
			};
	}

	void Sim::destroyByStrainMagnitude(float segStrainThreshold, float angleStrainThreshold) {
		std::function<bool(Sim*, Seg&, vec3)> segFunc = nullptr;
		if (glm::isfinite(segStrainThreshold)) {
			segStrainThreshold *= segStrainThreshold;
			segFunc = [segStrainThreshold](Sim* sim, Seg& seg, vec3 c) {
				return length2(c) < segStrainThreshold;
				};
		}

		std::function<bool(Sim*, RestAngle&, vec4)> angleFunc = nullptr;
		if (glm::isfinite(angleStrainThreshold)) {
			angleStrainThreshold *= angleStrainThreshold;
			angleFunc = [angleStrainThreshold](Sim* sim, RestAngle& angle, vec4 c) {
				return length2(c) < angleStrainThreshold;
				};
		}

		if (segFunc != nullptr || angleFunc != nullptr)
			deformDestroy(segFunc, angleFunc);
	}

	void updateId(vector<int>& indices, vector<int>& ids, int index, int from, int to) {
		int start = indices[index];
		int end = indices[index + 1];
		for (int i = start; i < end; i++)
			if (ids[i] == from) {
				ids[i] = to;
				break;
			}
	}

	void Sim::deformDestroy(
		std::function<bool(Sim*, Seg&, vec3)> segPredicate,
		std::function<bool(Sim*, RestAngle&, vec4)> anglePredicate) {
		// Removes any segments or angles that doesnt pass the test
		// Reorders the segments and angles to remove gaps
		// Updates references from vertices to segments and segments to angles
		// Super over complicated because all this is done in place

		// Process segments
		bool segRemoved = false;
		if (segPredicate != nullptr) {
#pragma omp parallel for if(multithreaded) schedule(dynamic, COSSERAT_THREAD_BATCH)
			for (int tid = 0; tid < segs.size(); tid++) {
				auto seg = segs[tid];

				// Compute strain
				auto v0 = verts[seg.i.x].pos;
				auto v1 = verts[seg.i.y].pos;
				float invl = 1 / seg.l;
				vec3 c = (v1 - v0) * invl - seg.q * vec3(1, 0, 0);

				bool keep = segPredicate(this, segs[tid], c);
				if (!keep) {
					segRemoved = true;

					// We dont want to shuffle the data around because removing angles need to reference these indices.
					// Set references from vertices to -1 for removal
					updateId(vertSegIndices, vertSegIds, seg.i.x, tid, -1);
					updateId(vertSegIndices, vertSegIds, seg.i.y, tid, -1);

					// Set referenced angles for removal
					int start = segAngleIndices[tid];
					int end = segAngleIndices[tid + 1];
					for (int i = start; i < end; i++) {
						auto& angle = restAngles[segAngleIds[i]].i;
						if (angle.x == tid) angle.x = -1;
						if (angle.y == tid) angle.y = -1;
					}

					// Mark for deletion
					segs[tid].i.x = -1;
				}
			}
		}

#pragma omp parallel for if(multithreaded) schedule(dynamic, COSSERAT_THREAD_BATCH)
		for (int tid = 0; tid < restAngles.size(); tid++) {
			auto angle = restAngles[tid];

			bool keep = angle.i.x >= 0 && angle.i.y >= 0;
			if (keep && anglePredicate != nullptr) {
				// Compute strain
				auto seg0 = segs[angle.i.x];
				auto seg1 = segs[angle.i.y];

				Rotor qq = seg0.q.inverse() * seg1.q;
				float phi = (dot(qq.v, angle.qRest.v) > 0) ? 1 : -1;
				vec4 c = (qq.v - phi * angle.qRest.v) * (4 / (seg0.l + seg1.l));

				if (!anglePredicate(this, restAngles[tid], c))
					keep = false;
			}

			if (!keep) {
				// Set references to -1 for removal if the segment still exists
				if (angle.i.x >= 0) updateId(segAngleIndices, segAngleIds, angle.i.x, tid, -1);
				if (angle.i.y >= 0) updateId(segAngleIndices, segAngleIds, angle.i.y, tid, -1);
				restAngles[tid].i.x = -1;	// Mark for deletion
			}
		}

		// If segments were removed, all we have to do is the following 8 things:
		// 1. Remove any vertices that dont have any connected segments
		// 2. Compact/move vertex data to remove gaps
		// 3. Compact/move index tables from vertices to segments with the vertex
		// 4. Update vertex indices stored in segments because those changed
		// 5. Compact/move segs data to remove gaps
		// 6. Compact/move index tables from segments to angles with the segment
		// 7. Update index tables from vertices to segments because the segments moved
		// 8. Update segment indices stored in angles because the segments moved

		// If angles were removed:
		// 1. Compact/move angle data to remove gaps
		// 2. Update index tables from segments to angles because the angles changed

		// The following does the 10 easy steps above, all at once, in 3 loops.

		// Process rest angles
		bool angleRemoved = false;
		{
			int curIndex = 0;
			for (int tid = 0; tid < restAngles.size(); tid++) {
				auto angle = restAngles[tid];
				bool keep = angle.i.x >= 0 && angle.i.y >= 0;

				// Deleted objects are just overwritten.
				if (keep) {
					if (curIndex != tid) {
						// Update segment indices to this
						updateId(segAngleIndices, segAngleIds, angle.i.x, tid, curIndex);
						updateId(segAngleIndices, segAngleIds, angle.i.y, tid, curIndex);

						// Overwrite this segment
						restAngles[curIndex] = restAngles[tid];
					}
					curIndex++;
				}
			}
			angleRemoved = restAngles.size() != curIndex;
			restAngles.resize(curIndex);
		}

		// Compact references from vertices to segments
		if (segRemoved) {
			int curIndex = 0;
			int curVertIndex = 0;
			for (int i = 0; i < verts.size(); i++) {
				if (i != curVertIndex) {
					// Move vertex entry to the new index
					verts[curVertIndex] = verts[i];
					vels[curVertIndex] = vels[i];
					lastVels[curVertIndex] = lastVels[i];
					vertSegIndices[curVertIndex] = curIndex;
				}

				int start = vertSegIndices[i];
				int end = vertSegIndices[i + 1];
				vertSegIndices[i] = curIndex;
				int oldCurIndex = curIndex;

				for (int j = start; j < end; j++) {
					int idx = vertSegIds[j];
					if (idx >= 0) {
						if (i != curVertIndex) {
							// Update segment-vertex reference to point to the new index
							auto& seg = segs[idx];
							if (seg.i.x == i) seg.i.x = curVertIndex;
							if (seg.i.y == i) seg.i.y = curVertIndex;
						}

						vertSegIds[curIndex++] = idx;
					}
				}

				// Only keep this vertex if it has segments
				if (oldCurIndex != curIndex)
					curVertIndex++;
			}
			verts.resize(curVertIndex);
			vertSegIndices.resize(verts.size() + 1);
			vertSegIds.resize(curIndex);
			vertSegIndices[verts.size()] = curIndex;
		}

		// Compact segment references everywhere
		if (segRemoved || angleRemoved) {
			int curSegIndex = 0;
			int curAngleIndex = 0;
			for (int i = 0; i < segs.size(); i++) {
				auto seg = segs[i];

				// If we are keeping this segment
				if (seg.i.x >= 0) {
					if (i != curSegIndex) {
						// Move segment entry to the new index
						segs[curSegIndex] = seg;
						xpbdSegWs[curSegIndex] = xpbdSegWs[i];
						lambdaTradGammas[curSegIndex] = lambdaTradGammas[i];

						// Update vertex-segment reference to point to the new index
						updateId(vertSegIndices, vertSegIds, seg.i.x, i, curSegIndex);
						updateId(vertSegIndices, vertSegIds, seg.i.y, i, curSegIndex);
					}

					int start = segAngleIndices[i];
					int end = segAngleIndices[i + 1];
					if (angleRemoved || i != curSegIndex) {
						// Update and move segment to angle reference
						segAngleIndices[curSegIndex] = curAngleIndex;

						// We need to do this sweep if any angles were removed to compact it
						// We also need to do this sweep if the segment was moved to copy it
						for (int j = start; j < end; j++) {
							int idx = segAngleIds[j];
							if (idx >= 0) {
								if (i != curSegIndex) {
									// Update angle-segment reference to point to the new index
									auto& angle = restAngles[idx].i;
									if (angle.x == i) angle.x = curSegIndex;
									if (angle.y == i) angle.y = curSegIndex;
								}

								// Copy ids to new location
								segAngleIds[curAngleIndex++] = idx;
							}
						}
					}
					else curAngleIndex += end - start;	// If the segment was not moved, we can just increment the index

					curSegIndex++;
				}
			}
			segs.resize(curSegIndex);
			segAngleIndices.resize(segs.size() + 1);
			segAngleIds.resize(curAngleIndex);
			segAngleIndices[segs.size()] = curAngleIndex;

			// Resize all the memory
			vels.resize(verts.size(), vec3(0));
			lastVels.resize(verts.size(), vec3(0));
			dxs.resize(verts.size(), vec3(0));

			xpbdSegLambdas.resize(segs.size(), vec3(0));
			xpbdAngleLambdas.resize(restAngles.size(), vec4(0));
			xpbdLastQs.resize(segs.size(), vec3(0));
			xpbdSegWs.resize(segs.size(), vec3(0));
			lambdaTradGammas.resize(segs.size(), 1.f);
		}
	}
}