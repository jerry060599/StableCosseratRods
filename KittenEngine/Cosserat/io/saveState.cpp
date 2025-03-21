#include "../Cosserat.h"
#include "gzstream.h"

namespace Cosserat {
	constexpr uint32_t MAGIC_NUMBER = 0xc022e5a7;	// Spells out "Cosserat" in hex

	struct CosseratHeader {
		uint32_t magic = MAGIC_NUMBER;	// Spells out "Cosserat" in hex
		int numVertices;
		int numSegments;
		int numRestAngles;

		uint32_t integrationMethod;
		float maxH;
		int numItr;
		float time;

		float renderByDensity;
		float renderDensity;
		float renderRadius;
		float damping;

		vec3 gravity;
		float drag;

		int numStateArrays;				// Number of extra data arrays saved outside of geometry
	};

	void Sim::saveToBinary(const char* filename, bool saveOnlyGeomtry) {
		using namespace gzstream;
		ogzstream file(filename, std::ios::out | std::ios::binary); // Open compressed file in binary mode
		if (file.fail()) throw std::runtime_error("Failed to open file!");

		CosseratHeader header{};
		header.numVertices = verts.size();
		header.numSegments = segs.size();
		header.numRestAngles = restAngles.size();

		header.integrationMethod = (uint32_t)method;
		header.maxH = maxH;
		header.numItr = meta.numItr;
		header.time = meta.time;

		header.renderByDensity = renderByDensity;
		header.renderDensity = renderDensity;
		header.renderRadius = renderRadius;
		header.damping = meta.damping;

		header.gravity = meta.gravity;
		header.drag = meta.drag;

		header.numStateArrays = saveOnlyGeomtry ? 0 : 4;

		// Write main data
		file.write(reinterpret_cast<const char*>(&header), sizeof(CosseratHeader));
		file.write(reinterpret_cast<const char*>(verts.data()), sizeof(Vert) * verts.size());
		file.write(reinterpret_cast<const char*>(segs.data()), sizeof(Seg) * segs.size());
		file.write(reinterpret_cast<const char*>(restAngles.data()), sizeof(RestAngle) * restAngles.size());

		if (!saveOnlyGeomtry) {
			// Write lastVel, vels, xpbdSegWs, and lambdaTradGammas
			// The rest can be regenerated
			file.write(reinterpret_cast<const char*>(vels.data()), sizeof(vec3) * vels.size());
			file.write(reinterpret_cast<const char*>(lastVels.data()), sizeof(vec3) * lastVels.size());
			file.write(reinterpret_cast<const char*>(xpbdSegWs.data()), sizeof(vec3) * xpbdSegWs.size());
			file.write(reinterpret_cast<const char*>(lambdaTradGammas.data()), sizeof(float) * lambdaTradGammas.size());
		}

		file.close();
	}

	void Sim::loadFromBinary(const char* filename, bool loadOnlyGeometry) {
		using namespace gzstream;
		igzstream file(filename, std::ios::in | std::ios::binary); // Open compressed file in binary mode
		if (file.fail()) throw std::runtime_error("Failed to open file!");

		CosseratHeader header{};
		file.read(reinterpret_cast<char*>(&header), sizeof(CosseratHeader));

		if (header.magic != MAGIC_NUMBER) {
			throw std::runtime_error("Not a Cosserat save state file!");
			file.close();
			return;
		}

		// Restore the simulation
		method = (IntegrationMethod)header.integrationMethod;
		maxH = header.maxH;
		meta.numItr = header.numItr;
		meta.time = header.time;

		renderByDensity = header.renderByDensity;
		renderDensity = header.renderDensity;
		renderRadius = header.renderRadius;
		meta.damping = header.damping;

		meta.gravity = header.gravity;
		meta.drag = header.drag;

		// Restore sim data

		verts.resize(header.numVertices);
		segs.resize(header.numSegments);
		restAngles.resize(header.numRestAngles);

		file.read(reinterpret_cast<char*>(verts.data()), sizeof(Vert) * verts.size());
		file.read(reinterpret_cast<char*>(segs.data()), sizeof(Seg) * segs.size());
		file.read(reinterpret_cast<char*>(restAngles.data()), sizeof(RestAngle) * restAngles.size());

		vels.resize(verts.size(), vec3(0));

		if (!loadOnlyGeometry) {
			lastVels.resize(verts.size(), vec3(0));
			xpbdSegWs.resize(segs.size(), vec3(0));
			lambdaTradGammas.resize(segs.size(), 1.f);

			if (header.numStateArrays >= 1) file.read(reinterpret_cast<char*>(vels.data()), sizeof(vec3) * vels.size());
			if (header.numStateArrays >= 2) file.read(reinterpret_cast<char*>(lastVels.data()), sizeof(vec3) * lastVels.size());
			if (header.numStateArrays >= 3) file.read(reinterpret_cast<char*>(xpbdSegWs.data()), sizeof(vec3) * xpbdSegWs.size());
			if (header.numStateArrays >= 4) file.read(reinterpret_cast<char*>(lambdaTradGammas.data()), sizeof(float) * lambdaTradGammas.size());
		}

		file.close();
	}
}