#pragma once
// Jerry Hsu 2022

#include <unordered_map>
#include <embree3/rtcore.h>
#include <functional>

#include "../includes/modules/Mesh.h"

namespace Kitten {
	class MeshCCD {
	public:
		typedef struct {
			RTCGeometry rtcTriGeom;
			unsigned int triGeomID;

			RTCGeometry rtcEdgeGeom;
			unsigned int edgeGeomID;

			RTCGeometry rtcVertGeom;
			unsigned int vertGeomID;

			Kitten::Mesh* mesh;
			glm::vec3* delta;
			glm::ivec2* edges;
			int nEdges;
			bool dirty;
		} RTCMesh;

		typedef struct {
			// The mesh the triangle belongs to
			Kitten::Mesh* triMesh;
			// The mesh the vertex belogns to
			Kitten::Mesh* vertMesh;

			// The id of the triangle
			int triIndex;
			// The id of the vertex
			int vertIndex;

			// The time from 0 to 1 this collsion happens at
			float t;
			// The bary centric coordinates of the collision point on the triangle
			glm::vec3 bary;
			// The collision normal (pointing away from tri)
			glm::vec3 norm;
		} TriVertCollision;

		typedef struct {
			// The mesh the 'a' edge belongs to
			Kitten::Mesh* aMesh;
			// The mesh the 'b' edge belongs to
			Kitten::Mesh* bMesh;

			// The vertex ids of the end points of edge a.
			glm::ivec2 ai;
			// The vertex ids of the end points of edge b.
			glm::ivec2 bi;

			// The time from 0 to 1 this collsion happens at
			float t;
			// The position of the collision encoded as u in [0, 1] from a0 to a1 and v in [0, 1] from b0 to b1
			glm::vec2 uv;
			// The collision normal (pointing away from a)
			glm::vec3 norm;
		} EdgeEdgeCollision;

	private:
		RTCDevice rtcDevice;

		RTCScene rtcTriScene;
		RTCScene rtcEdgeScene;
		RTCScene rtcVertScene;

		std::unordered_map<Kitten::Mesh*, RTCMesh*> meshes;
		void triVertCCD(struct RTCCollision*, unsigned int, std::function<void(TriVertCollision)> triVertColCallback);
		void edgeEdgeCCD(struct RTCCollision*, unsigned int, std::function<void(EdgeEdgeCollision)> edgeEdgeColCallback);

	public:
		MeshCCD();
		~MeshCCD();

		glm::vec3*& operator[](Kitten::Mesh*);

		void attachStatic(Kitten::Mesh* mesh);
		void attach(Kitten::Mesh* mesh, glm::vec3* delta);
		void detach(Kitten::Mesh*);

		// Marks a mesh as dirty and update its bvh as such.
		void dirty(Kitten::Mesh* mesh);

		// Rebuilds the bvh. Needed if the bounding box of any triangles are changed.
		void rebuildBVH();

		// Performs actual collision detection using the bvh
		void collide(std::function<void(TriVertCollision)> triVertColCallback, std::function<void(EdgeEdgeCollision)> edgeEdgeColCallback);
	};
}