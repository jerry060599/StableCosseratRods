
#include "MeshCCD.h"
#include "../includes/modules/Bound.h"
#include "../opt/polynomial.h"
#include <unordered_set>

namespace Kitten {
	// Determines when the tetrahedra formed by the four points in p becomes degenerate
	int planarMovingPoints(const mat4x3& p, const mat4x3& d, vec3& t) {
		dvec3 e0 = dvec3(p[1]) - dvec3(p[0]);
		dvec3 e1 = dvec3(p[2]) - dvec3(p[0]);
		dvec3 e2 = dvec3(p[3]) - dvec3(p[0]);

		dvec3 d0 = dvec3(d[1]) - dvec3(d[0]);
		dvec3 d1 = dvec3(d[2]) - dvec3(d[0]);
		dvec3 d2 = dvec3(d[3]) - dvec3(d[0]);

		dvec3 e0xe1 = cross(e0, e1);
		dvec3 d0xd1 = cross(d0, d1);
		dvec3 cs = cross(e0, d1) + cross(d0, e1);

		double poly[4];
		poly[0] = dot(e0xe1, e2);
		poly[1] = dot(d2, e0xe1) + dot(e2, cs);
		poly[2] = dot(e2, d0xd1) + dot(d2, cs);
		poly[3] = dot(d0xd1, d2);

		dvec3 dt;
		int n = cy::PolynomialRoots<3>((double*)&dt, poly, -1e-6, 1 + 1e-6, 1e-9);

		t = vec3(dt);
		return n;
	}

	// ccd between triangle formed by points 0,1,2 and point 3.
	bool intMovingTriPoint(const mat4x3& points, const mat4x3& deltas, mat4x3& x, float& t, vec3& bary) {
		vec3 ts;
		int nt = planarMovingPoints(points, deltas, ts);
		for (int i = 0; i < nt; i++) {
			x = points + deltas * ts[i];
			bary = baryCoord(x);
			if (all(greaterThanEqual(bary, vec3(-5 * numeric_limits<float>::epsilon())))) {
				t = ts[i];
				return true;
			}
		}
		return false;
	}

	inline vec3 largestMag(vec3 a, vec3 b) {
		return length2(a) > length2(b) ? a : b;
	}

	// ccd between edges formed by points 0,1 and 2,3.
	bool intMovingEdgeEdge(const mat4x3& points, const mat4x3& deltas, mat4x3& x, float& t, vec2& uv, vec3& norm) {
		vec3 ts;
		int nt = planarMovingPoints(points, deltas, ts);
		for (int i = 0; i < nt; i++) {
			x = points + deltas * ts[i];
			vec3 ad = x[1] - x[0];
			vec3 bd = x[3] - x[2];
			vec3 diff = x[2] - x[0];

			norm = largestMag(cross(ad, bd), cross(ad, diff));
			vec3 an = cross(norm, ad);
			vec3 bn = cross(norm, bd);

			float ax = dot(diff, an);
			ax *= (ax + dot(bd, an));
			float bx = -dot(diff, bn);
			bx *= (bx + dot(ad, bn));

			if (ax <= 0 && bx <= 0 && (ax != 0 || bx != 0)) {
				uv = lineClosestPoints(x[0], x[1], x[2], x[3]);
				return true;
			}
		}
		return false;
	}

	void rtcErrorFunc(void* userPtr, enum RTCError error, const char* str) {
		printf("error %d: %s\n", error, str);
	}

	MeshCCD::MeshCCD() {
		rtcDevice = rtcNewDevice(nullptr);
		if (!rtcDevice)
			printf("error %d: cannot create device\n", rtcGetDeviceError(NULL));
		rtcSetDeviceErrorFunction(rtcDevice, rtcErrorFunc, NULL);

		rtcTriScene = rtcNewScene(rtcDevice);
		rtcEdgeScene = rtcNewScene(rtcDevice);
		rtcVertScene = rtcNewScene(rtcDevice);
	}

	MeshCCD::~MeshCCD() {
		for (auto p : meshes) {
			rtcDetachGeometry(rtcTriScene, p.second->triGeomID);
			rtcDetachGeometry(rtcEdgeScene, p.second->edgeGeomID);
			rtcDetachGeometry(rtcVertScene, p.second->vertGeomID);
			delete[] p.second->edges;
			delete p.second;
		}
		rtcReleaseScene(rtcTriScene);
		rtcReleaseScene(rtcEdgeScene);
		rtcReleaseScene(rtcVertScene);
		rtcReleaseDevice(rtcDevice);
	}

	ivec2* calcEdgeList(Kitten::Mesh* mesh, int& nEdges) {
		std::unordered_set<pair<int, int>> edgeSet;
		vector<ivec2> edges;
		for (size_t i = 0; i < mesh->indices.size() / 3; i++)
			for (int k = 0; k < 3; k++) {
				auto e = std::make_pair(mesh->indices[k + 3 * i], mesh->indices[((k + 1) % 3) + 3 * i]);
				if (e.first > e.second) std::swap(e.first, e.second);
				if (!edgeSet.count(e)) {
					edges.push_back({ e.first, e.second });
					edgeSet.insert(e);
				}
			}

		ivec2* dat = new ivec2[edges.size()];
		memcpy(dat, &edges[0], sizeof(ivec2) * edges.size());
		nEdges = (int)edges.size();
		return dat;
	}

	glm::vec3*& MeshCCD::operator[](Kitten::Mesh* mesh) {
		auto itr = meshes.find(mesh);
		if (itr != meshes.end())
			return itr->second->delta;

		RTCMesh* ptr = new RTCMesh;
		*ptr = {};
		ptr->mesh = mesh;
		ptr->edges = calcEdgeList(mesh, ptr->nEdges);
		meshes[mesh] = ptr;
		return ptr->delta;
	}

	void MeshCCD::attachStatic(Kitten::Mesh* mesh) {
		attach(mesh, nullptr);
	}

	void MeshCCD::attach(Kitten::Mesh* mesh, glm::vec3* delta) {
		(*this)[mesh] = delta;
	}

	void MeshCCD::detach(Kitten::Mesh* mesh) {
		auto itr = meshes.find(mesh);
		rtcDetachGeometry(rtcTriScene, itr->second->triGeomID);
		rtcDetachGeometry(rtcEdgeScene, itr->second->edgeGeomID);
		rtcDetachGeometry(rtcVertScene, itr->second->vertGeomID);
		delete[] itr->second->edges;
		delete itr->second;
		meshes.erase(itr);
	}

	void MeshCCD::dirty(Kitten::Mesh* mesh) {
		meshes[mesh]->dirty = true;
	}

	void triBoundFunc(const struct RTCBoundsFunctionArguments* args) {
		MeshCCD::RTCMesh& data = *(MeshCCD::RTCMesh*)args->geometryUserPtr;
		Kitten::Mesh& mesh = *data.mesh;

		ivec3 i = ivec3(0, 1, 2) + (int)(3 * args->primID);
		i = ivec3(mesh.indices[i[0]], mesh.indices[i[1]], mesh.indices[i[2]]);

		Kitten::Bound<> bound(mesh.vertices[i.x].pos);
		bound.absorb(mesh.vertices[i.y].pos);
		bound.absorb(mesh.vertices[i.z].pos);

		if (data.delta) {
			bound.absorb(mesh.vertices[i.x].pos + data.delta[i.x]);
			bound.absorb(mesh.vertices[i.y].pos + data.delta[i.y]);
			bound.absorb(mesh.vertices[i.z].pos + data.delta[i.z]);
		}

		*(vec3*)&args->bounds_o->lower_x = bound.min;
		*(vec3*)&args->bounds_o->upper_x = bound.max;
	}

	void edgeBoundFunc(const struct RTCBoundsFunctionArguments* args) {
		MeshCCD::RTCMesh& data = *(MeshCCD::RTCMesh*)args->geometryUserPtr;
		Kitten::Mesh& mesh = *data.mesh;

		ivec2 i = data.edges[args->primID];

		Kitten::Bound<> bound(mesh.vertices[i.x].pos);
		bound.absorb(mesh.vertices[i.y].pos);

		if (data.delta) {
			bound.absorb(mesh.vertices[i.x].pos + data.delta[i.x]);
			bound.absorb(mesh.vertices[i.y].pos + data.delta[i.y]);
		}

		*(vec3*)&args->bounds_o->lower_x = bound.min;
		*(vec3*)&args->bounds_o->upper_x = bound.max;
	}

	void vertBoundFunc(const struct RTCBoundsFunctionArguments* args) {
		MeshCCD::RTCMesh& data = *(MeshCCD::RTCMesh*)args->geometryUserPtr;
		Kitten::Mesh& mesh = *data.mesh;

		Kitten::Bound<> bound(mesh.vertices[args->primID].pos);
		if (data.delta) bound.absorb(mesh.vertices[args->primID].pos + data.delta[args->primID]);

		*(vec3*)&args->bounds_o->lower_x = bound.min;
		*(vec3*)&args->bounds_o->upper_x = bound.max;
	}

	void MeshCCD::rebuildBVH() {
		// Loop through every mesh
		for (auto& pair : meshes) {
			// Allocate rtcGeom and attach to scene if needed.
			if (pair.second->rtcTriGeom == nullptr) {
				// Allocate tri
				pair.second->rtcTriGeom = rtcNewGeometry(rtcDevice, RTC_GEOMETRY_TYPE_USER);
				rtcSetGeometryUserPrimitiveCount(pair.second->rtcTriGeom, (unsigned int)(pair.first->indices.size() / 3));
				rtcSetGeometryUserData(pair.second->rtcTriGeom, pair.second);
				rtcSetGeometryBoundsFunction(pair.second->rtcTriGeom, triBoundFunc, pair.first);
				rtcSetGeometryIntersectFunction(pair.second->rtcTriGeom, nullptr);
				rtcSetGeometryOccludedFunction(pair.second->rtcTriGeom, nullptr);

				// Allocate edge
				pair.second->rtcEdgeGeom = rtcNewGeometry(rtcDevice, RTC_GEOMETRY_TYPE_USER);
				rtcSetGeometryUserPrimitiveCount(pair.second->rtcEdgeGeom, (unsigned int)pair.second->nEdges);
				rtcSetGeometryUserData(pair.second->rtcEdgeGeom, pair.second);
				rtcSetGeometryBoundsFunction(pair.second->rtcEdgeGeom, edgeBoundFunc, pair.first);
				rtcSetGeometryIntersectFunction(pair.second->rtcEdgeGeom, nullptr);
				rtcSetGeometryOccludedFunction(pair.second->rtcEdgeGeom, nullptr);

				// Allocate vert
				pair.second->rtcVertGeom = rtcNewGeometry(rtcDevice, RTC_GEOMETRY_TYPE_USER);
				rtcSetGeometryUserPrimitiveCount(pair.second->rtcVertGeom, (unsigned int)pair.first->vertices.size());
				rtcSetGeometryUserData(pair.second->rtcVertGeom, pair.second);
				rtcSetGeometryBoundsFunction(pair.second->rtcVertGeom, vertBoundFunc, pair.first);
				rtcSetGeometryIntersectFunction(pair.second->rtcVertGeom, nullptr);
				rtcSetGeometryOccludedFunction(pair.second->rtcVertGeom, nullptr);

				// Attach
				rtcCommitGeometry(pair.second->rtcTriGeom);
				pair.second->triGeomID = rtcAttachGeometry(rtcTriScene, pair.second->rtcTriGeom);
				rtcReleaseGeometry(pair.second->rtcTriGeom);

				rtcCommitGeometry(pair.second->rtcEdgeGeom);
				pair.second->edgeGeomID = rtcAttachGeometry(rtcEdgeScene, pair.second->rtcEdgeGeom);
				rtcReleaseGeometry(pair.second->rtcEdgeGeom);

				rtcCommitGeometry(pair.second->rtcVertGeom);
				pair.second->vertGeomID = rtcAttachGeometry(rtcVertScene, pair.second->rtcVertGeom);
				rtcReleaseGeometry(pair.second->rtcVertGeom);

				pair.second->dirty = false;
			}

			if (pair.second->dirty) {
				pair.second->dirty = false;
				// Recommit
				rtcUpdateGeometryBuffer(pair.second->rtcTriGeom, RTC_BUFFER_TYPE_VERTEX, 0);
				rtcCommitGeometry(pair.second->rtcTriGeom);
				rtcUpdateGeometryBuffer(pair.second->rtcEdgeGeom, RTC_BUFFER_TYPE_VERTEX, 0);
				rtcCommitGeometry(pair.second->rtcEdgeGeom);
				rtcUpdateGeometryBuffer(pair.second->rtcVertGeom, RTC_BUFFER_TYPE_VERTEX, 0);
				rtcCommitGeometry(pair.second->rtcVertGeom);
			}
		}

		rtcCommitScene(rtcTriScene);
		rtcCommitScene(rtcEdgeScene);
		rtcCommitScene(rtcVertScene);
	}

	void MeshCCD::triVertCCD(struct RTCCollision* collisions, unsigned int num_collisions, std::function<void(TriVertCollision)> triVertColCallback) {
		for (int i = 0; i < (int)num_collisions; i++) {
			const RTCCollision& col = collisions[i];
			const RTCMesh& rtcTri = *(RTCMesh*)rtcGetGeometryUserData(rtcGetGeometry(rtcTriScene, col.geomID0));
			const RTCMesh& rtcVert = *(RTCMesh*)rtcGetGeometryUserData(rtcGetGeometry(rtcVertScene, col.geomID1));

			Mesh& triMesh = *rtcTri.mesh;
			Mesh& vertMesh = *rtcVert.mesh;

			ivec3 triIndices = ivec3(0, 1, 2) + 3 * (int)col.primID0;
			triIndices = ivec3(triMesh.indices[triIndices[0]], triMesh.indices[triIndices[1]], triMesh.indices[triIndices[2]]);

			if (&triMesh == &vertMesh && any(equal(triIndices, ivec3(col.primID1)))) continue;

			const mat4x3 pos(
				triMesh.vertices[triIndices[0]].pos,
				triMesh.vertices[triIndices[1]].pos,
				triMesh.vertices[triIndices[2]].pos,
				vertMesh.vertices[col.primID1].pos
			);

			mat4x3 delta(0);
			if (rtcTri.delta) {
				delta[0] = rtcTri.delta[triIndices[0]];
				delta[1] = rtcTri.delta[triIndices[1]];
				delta[2] = rtcTri.delta[triIndices[2]];
			}
			if (rtcVert.delta)
				delta[3] = rtcVert.delta[col.primID1];

			mat4x3 cur;
			TriVertCollision dat;
			if (intMovingTriPoint(pos, delta, cur, dat.t, dat.bary)) {
				if (!all(isfinite(dat.bary))) continue;

				// Found intersection
				dat.norm = cross(cur[1] - cur[0], cur[2] - cur[0]);
				vec3 relDelta = delta * vec4(dat.bary, -1);

				float s = dot(dat.norm, relDelta);
				if (s < 0) dat.norm *= -1;
				else if (s == 0) continue; // Discard in favor of edge-edge collision.

				float d = length2(dat.norm);
				if (d > 0) dat.norm *= inversesqrt(d);
				else continue;

				// call callback
				dat.triMesh = &triMesh;
				dat.vertMesh = &vertMesh;

				dat.triIndex = col.primID0;
				dat.vertIndex = col.primID1;

				triVertColCallback(dat);
			}
		}
	}

	void MeshCCD::edgeEdgeCCD(struct RTCCollision* collisions, unsigned int num_collisions, std::function<void(EdgeEdgeCollision)> edgeEdgeColCallback) {
		for (int i = 0; i < (int)num_collisions; i++) {
			const RTCCollision& col = collisions[i];
			if (col.geomID0 == col.geomID1 && col.primID0 == col.primID1) continue;
			const RTCMesh& artc = *(RTCMesh*)rtcGetGeometryUserData(rtcGetGeometry(rtcEdgeScene, col.geomID0));
			const RTCMesh& brtc = *(RTCMesh*)rtcGetGeometryUserData(rtcGetGeometry(rtcEdgeScene, col.geomID1));
			Mesh& a = *artc.mesh;
			Mesh& b = *brtc.mesh;

			const ivec2 aIndices = artc.edges[col.primID0];
			const ivec2 bIndices = brtc.edges[col.primID1];

			if (&a == &b)
				if (any(equal(aIndices, ivec2(bIndices.x))) ||
					any(equal(aIndices, ivec2(bIndices.y)))) continue;

			const mat4x3 pos(
				a.vertices[aIndices[0]].pos,
				a.vertices[aIndices[1]].pos,
				b.vertices[bIndices[0]].pos,
				b.vertices[bIndices[1]].pos
			);

			mat4x3 delta(0);
			if (artc.delta) {
				delta[0] = artc.delta[aIndices[0]];
				delta[1] = artc.delta[aIndices[1]];
			}
			if (brtc.delta) {
				delta[2] = brtc.delta[bIndices[0]];
				delta[3] = brtc.delta[bIndices[1]];
			}

			mat4x3 cur;
			EdgeEdgeCollision dat;
			if (intMovingEdgeEdge(pos, delta, cur, dat.t, dat.uv, dat.norm)) {
				// Found intersection
				vec3 relDelta = delta * vec4(1 - dat.uv.x, dat.uv.x, dat.uv.y - 1, -dat.uv.y);

				float s = dot(dat.norm, relDelta);
				if (s < 0) dat.norm *= -1;
				else if (s == 0) {
					dat.uv = clamp(lineClosestPoints(pos[0], pos[1], pos[2], pos[3]), vec2(0), vec2(1));
					dat.norm = delta * vec4(1 - dat.uv.x, dat.uv.x, dat.uv.y - 1, -dat.uv.y);
				}

				if (!all(isfinite(dat.uv))) continue;

				float d = length2(dat.norm);
				if (d > 0) dat.norm *= inversesqrt(d);
				else continue;

				// call callback
				dat.aMesh = &a;
				dat.bMesh = &b;

				dat.ai = aIndices;
				dat.bi = bIndices;

				edgeEdgeColCallback(dat);
			}
		}
	}

	void MeshCCD::collide(std::function<void(TriVertCollision)> triVertColCallback, std::function<void(EdgeEdgeCollision)> edgeEdgeColCallback) {
		struct params {
			MeshCCD* ccd;
			std::function<void(TriVertCollision)> tvf;
			std::function<void(EdgeEdgeCollision)> eef;
		};
		params p{ this, triVertColCallback, edgeEdgeColCallback };

		rtcCollide(rtcTriScene, rtcVertScene, [](void* userPtr, struct RTCCollision* collisions, unsigned int num_collisions) {
			params& p = *(params*)userPtr;
			p.ccd->triVertCCD(collisions, num_collisions, p.tvf);
			}, &p);

		rtcCollide(rtcEdgeScene, rtcEdgeScene, [](void* userPtr, struct RTCCollision* collisions, unsigned int num_collisions) {
			params& p = *(params*)userPtr;
			p.ccd->edgeEdgeCCD(collisions, num_collisions, p.eef);
			}, &p);
	}
}