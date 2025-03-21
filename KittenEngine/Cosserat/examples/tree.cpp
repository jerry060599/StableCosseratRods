#include "../CosseratExamples.h"
#include "fastPRGN.h"

#include <bit>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <KittenEngine/includes/modules/StopWatch.h>
#include <KittenEngine/includes/modules/Algo.h>

namespace Cosserat {
	inline float randf(fastPRNG::fastXS64& rng) {
		return glm::clamp((uint32_t)rng.xoroshiro128p() / (float)UINT32_MAX, 0.f, 1.f);
	}

	float dfsTreeStr(int cur, vector<vec3>& verts, unordered_multimap<int, int>& conMap, vector<float>& masses) {
		masses[cur] = 0;
		float mass = 0;
		vec3 p = verts[cur];

		for (auto it = conMap.equal_range(cur); it.first != it.second; it.first++) {
			int next = it.first->second;
			mass += length(verts[next] - p);

			if (masses[next] < 0) // Unvisited
				mass += dfsTreeStr(next, verts, conMap, masses);
		}

		masses[cur] = mass;
		return mass;
	}

	ivec2 searchAttractionPoints(float r2, float killR2, int lastVertIdx, vec3 lastVert,
		vector<vec3>& attractionPoints, vector<tuple<float, int>>& closestPairs, vector<float>& vertLength) {
		// Get the closest attraction point to the closest vertex
		float minVertLen = std::numeric_limits<float>::infinity();
		ivec2 closestPair(-1);

		for (int i = 0; i < attractionPoints.size(); i++) {
			float d = length2(lastVert - attractionPoints[i]);

			// Kill attraction points within kill radius
			if (d <= killR2) {
				std::swap(attractionPoints[i], attractionPoints.back());
				attractionPoints.pop_back();
				std::swap(closestPairs[i], closestPairs.back());
				closestPairs.pop_back();
				i--;
				continue;
			}

			// Update closest vertex
			auto& [dist, idx] = closestPairs[i];
			if (d < dist) {
				dist = d;
				idx = lastVertIdx;
			}

			// If within attraction distance
			if (dist <= r2) {
				// Pick the one with the highest fitness
				float vertLen = vertLength[idx] + sqrt(dist);
				if (vertLen < minVertLen) {
					minVertLen = vertLen;
					closestPair = ivec2(idx, i);
				}
			}
		}

		// Exit if no more attraction points are close enough
		if (!glm::isfinite(minVertLen))
			return ivec2(-1);

		return closestPair;
	}

	Sim* generateTreeExample(int numTrees) {
		using namespace fastPRNG;

		// Example of trees generated using "Modeling Trees with a Space Colonization Algorithm", 2007.
		Sim* sim = new Cosserat::Sim();

		// Generation settings
		const int numAttractionPoints = 1500;
		const float attractionRadius = 0.6f;
		const float killRadius = 0.3f * attractionRadius;

		const float attractionRadius2 = Kit::pow2(attractionRadius);
		const float killRadius2 = Kit::pow2(killRadius);

		// Tree geometry settings
		const float maxSegLen = 0.2f;
		const float trunkHeight = 3.f;
		const float crownWidth = 4.f;
		const float crownLowerHeight = 4.0f;
		const float crownUpperHeight = 1.5f;
		const float trunkPen = 1.5f;
		vec3 branchBias = vec3(0, 0.7f, 0);

		// Tree location settings
		const float rootRadius = 1.1f * sqrt((float)numTrees) * crownWidth;

		// Material settings

		// Minimum density of branches kg per meter.
		// Assumes branches are 2cm in diameter with a density of 700kg/m^3.
		const float density = 3.1415f * Kit::pow2(0.02f) * 700.f;
		const float kss = 1e4;
		const float kbt = 50;

		// Realistically should be 2, but 1.8 looks like a tree.
		// This is we assume radius is sqrt(m) and moment is r^4 -> moment is m^2.
		const float momentPower = 1.8f;

		const int64_t seed = 0x12345678;

		Kit::StopWatch timer;
		vector<vec3> roots;
		{
			int numRootSample = numTrees * 10;

			fastXS64 prng(seed);
			vector<vec3> rootPoints;
			rootPoints.reserve(numRootSample);
			for (int i = 0; i < numRootSample; i++) {
				float r = sqrt(randf(prng)) * rootRadius;
				float theta = randf(prng) * 2 * glm::pi<float>();
				rootPoints.push_back(vec3(r * cos(theta), 0, r * sin(theta)));
			}

			roots = bluenoiseSample(rootPoints, numTrees);
		}
		timer.time("Root location selection");

#pragma omp parallel for schedule(dynamic, 1)
		for (int tree = 0; tree < numTrees; tree++) {
			fastXS64 prng(seed ^ tree);
			vec3 root = roots[tree];

			// Blue noise sample the tree crown
			vector<vec3> attractionPoints;
			const int numSamples = numAttractionPoints * 5;
			attractionPoints.reserve(numSamples);
			for (int i = 0; i < numSamples; i++) {
				// Random sphere
				float r = pow(randf(prng), 1 / 3.f);
				float theta = randf(prng) * 2 * glm::pi<float>();
				float phi = (2 * randf(prng) - 1) * glm::pi<float>();
				vec3 p = r * vec3(cos(theta) * cos(phi), sin(phi), sin(theta) * cos(phi));

				if (p.y > 0) p.y *= crownUpperHeight;
				else p.y *= crownLowerHeight;
				p.x *= crownWidth;
				p.z *= crownWidth;
				p.y += trunkHeight + crownLowerHeight - trunkPen;
				p += root;

				attractionPoints.push_back(p);
			}

			attractionPoints = bluenoiseSample(attractionPoints, numAttractionPoints);
			vector<tuple<float, int>> closestPairs(attractionPoints.size(),
				{ std::numeric_limits<float>::infinity(), -1 });

			// Build trunk
			vector<vec3> verts;
			vector<float> vertLength;
			vector<ivec2> segs;

			int trunkSegs = (int)ceil(trunkHeight / maxSegLen);
			for (int i = 0; i <= trunkSegs; i++) {
				if (i) searchAttractionPoints(attractionRadius2, killRadius2, verts.size() - 1, verts.back(), attractionPoints, closestPairs, vertLength);
				verts.push_back(vec3(root.x, trunkHeight * i / (float)trunkSegs, root.z));
				vertLength.push_back(0);
			}
			for (int i = 0; i < trunkSegs; i++)
				segs.push_back(ivec2(i, i + 1));

			// Space colonization
			while (true) {
				// Get the closest attraction point to the closest vertex
				ivec2 closestPair = searchAttractionPoints(attractionRadius2, killRadius2, verts.size() - 1, verts.back(),
					attractionPoints, closestPairs, vertLength);

				// Exit if no more attraction points are close enough
				if (closestPair.x < 0) break;

				// Add a segment
				vec3 p0 = verts[closestPair.x];
				vec3 p1 = attractionPoints[closestPair.y];

				vec3 n = p1 - p0;
				float l = length(n);
				n = normalize(n / l + branchBias);

				verts.push_back(p0 + maxSegLen * n);
				vertLength.push_back(vertLength[closestPair.x] + l);
				segs.push_back(ivec2(closestPair.x, verts.size() - 1));
			}

			// Process geometry to get branch mass
			unordered_multimap<int, int> conMap;
			for (auto s : segs) {
				conMap.insert({ s.x, s.y });
				conMap.insert({ s.y, s.x });
			}
			vector<float> masses(verts.size(), -1);
			float totalMass = dfsTreeStr(0, verts, conMap, masses);

			// Add rest angles between all connected segments
			unordered_multimap<int, int> segMap;
			for (int i = 0; i < segs.size(); i++) {
				auto s = segs[i];
				segMap.insert({ s.x, i });
				segMap.insert({ s.y, i });
			}
			vector<ivec2> angles;
			for (int i = 0; i < segs.size(); i++) {
				auto s = segs[i];
				for (auto it = segMap.equal_range(s.x); it.first != it.second; it.first++) {
					int other = it.first->second;
					if (other >= i) continue;
					angles.push_back(ivec2(other, i));
				}
				for (auto it = segMap.equal_range(s.y); it.first != it.second; it.first++) {
					int other = it.first->second;
					if (other >= i) continue;
					angles.push_back(ivec2(other, i));
				}
			}

#pragma omp critical
			{
				// Add to simulation
				int vertOffset = sim->verts.size();
				int segOffset = sim->segs.size();

				for (auto v : verts) sim->addVertex(v);
				for (auto s : segs) {
					float mass = (masses[s.x] + masses[s.y]) * 0.5f;
					float segLen = length(verts[s.y] - verts[s.x]);
					float localDensity = mass / segLen;
					sim->addSegment(ivec2(s.x + vertOffset, s.y + vertOffset),
						kss * localDensity, localDensity * density);
				}
				for (auto a : angles) {
					auto s0 = segs[a.x];
					auto s1 = segs[a.y];

					float avgMass = 0.25f * (masses[s0.x] + masses[s0.y] + masses[s1.x] + masses[s1.y]);
					sim->addRestAngle(ivec2(a.x + segOffset, a.y + segOffset), kbt * pow(avgMass, momentPower));
				}

				// Pin the root
				sim->pinVertex(vertOffset);
				sim->pinSegment(segOffset);

				printf("Generated tree %d with %f len, %d verts, %d segs, and %d angles.\n",
					tree, totalMass, verts.size(), segs.size(), angles.size());
			}
		}

		printf("Generated %d trees with %d verts, %d segs, and %d angles.\n", numTrees, sim->verts.size(), sim->segs.size(), sim->restAngles.size());
		timer.time("Tree generation");
		timer.printTimes();

		sim->maxH = 1 / 1000.f;
		sim->meta.numItr = 4;
		sim->meta.damping = 1e-6f;
		sim->meta.drag = 0.01f;
		sim->meshHeightLen = 0.25f * maxSegLen;
		sim->meshExcludeInlinePhong = 1;

		sim->renderByDensity = 1;
		sim->renderDensity = 8e4;

		return sim;
	}
}