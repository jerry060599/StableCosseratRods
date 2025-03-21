#include "../Cosserat.h"

namespace Cosserat {
	void Sim::addVertex(vec3 pos) {
		Vert v;
		v.pos = pos;
		v.invMass = numeric_limits<float>::infinity();
		vels.push_back(vec3(0));
		verts.push_back(v);
		vertexLastSegment.push_back(-1);
	}

	void Sim::addSegment(ivec2 i, float kss, float density) {
		if (i.x < 0 || i.x >= verts.size() ||
			i.y < 0 || i.y >= verts.size())
			throw invalid_argument("Segment vertex index out of bounds");

		Seg s;
		s.i = i;

		vec3 dir = verts[i.y].pos - verts[i.x].pos;
		s.l = length(dir);
		if (s.l <= 0) throw invalid_argument("Invalid segment length");

		s.q = Rotor::fromTo(vec3(1, 0, 0), dir / s.l);
		s.kss = kss * s.l;

		// Minimize twisting
		if (vertexLastSegment[i.x] >= 0) {
			auto lastSeg = segs[vertexLastSegment[i.x]];
			auto t = s.q.inverse() * lastSeg.q;
			if (t.x != 0 || t.w != 0) {
				t = normalize(vec4(t.x, 0, 0, t.w));
				s.q = s.q * t;
			}
		}
		else if (vertexLastSegment[i.y] >= 0) {
			auto lastSeg = segs[vertexLastSegment[i.y]];
			auto t = s.q.inverse() * lastSeg.q;
			if (t.x != 0 || t.w != 0) {
				t = normalize(vec4(t.x, 0, 0, t.w));
				s.q = s.q * t;
			}
		}

		vertexLastSegment[i.x] = segs.size();
		vertexLastSegment[i.y] = segs.size();

		segs.push_back(s);

		const float m = 0.5f * density * s.l;
		if (verts[i.x].invMass > 0) verts[i.x].invMass = 1 / (1 / verts[i.x].invMass + m);
		if (verts[i.y].invMass > 0) verts[i.y].invMass = 1 / (1 / verts[i.y].invMass + m);
	}

	void Sim::addRestAngle(ivec2 i, float kbt) {
		if (i.x < 0 || i.x >= segs.size() ||
			i.y < 0 || i.y >= segs.size())
			throw invalid_argument("Angle segment index out of bounds");
		if (i.x == i.y) throw invalid_argument("Angle segment indices are the same");

		RestAngle ra;
		ra.i = i;
		ra.qRest = segs[i.x].q.inverse() * segs[i.y].q;
		ra.kbt = 4 * kbt / (0.5f * (segs[i.x].l + segs[i.y].l));
		restAngles.push_back(ra);
	}

	void Sim::pinVertex(int i) {
		if (i < 0 || i >= verts.size())
			throw invalid_argument("Vertex index out of bounds");
		verts[i].invMass = 0;
	}

	void Sim::pinSegment(int i) {
		if (i < 0 || i >= segs.size())
			throw invalid_argument("Segment index out of bounds");
		segs[i].kss = -abs(segs[i].kss);
	}

	void Sim::init() {
		// Initialize vertex-segment indices
		vector<int> counts(verts.size() + 1, 0);
		for (auto& seg : segs) {
			counts[seg.i.x]++;
			counts[seg.i.y]++;
		}

		int runningSum = 0;
		vertSegIndices.resize(verts.size() + 1);
		for (size_t i = 0; i < vertSegIndices.size(); i++) {
			vertSegIndices[i] = runningSum;
			runningSum += counts[i];
			counts[i] = 0;
		}

		vertSegIds.resize(runningSum);
		for (size_t i = 0; i < segs.size(); i++) {
			ivec2 ind = segs[i].i;
			vertSegIds[vertSegIndices[ind.x] + counts[ind.x]++] = i;
			vertSegIds[vertSegIndices[ind.y] + counts[ind.y]++] = i;
		}

		// Initialize segment-rest angle indices
		counts.resize(segs.size() + 1);
		for (int& c : counts) c = 0;
		for (auto& ra : restAngles) {
			counts[ra.i.x]++;
			counts[ra.i.y]++;
		}

		runningSum = 0;
		segAngleIndices.resize(segs.size() + 1);
		for (size_t i = 0; i < segAngleIndices.size(); i++) {
			segAngleIndices[i] = runningSum;
			runningSum += counts[i];
			counts[i] = 0;
		}

		segAngleIds.resize(runningSum);
		for (size_t i = 0; i < restAngles.size(); i++) {
			ivec2 ind = restAngles[i].i;
			segAngleIds[segAngleIndices[ind.x] + counts[ind.x]++] = i;
			segAngleIds[segAngleIndices[ind.y] + counts[ind.y]++] = i;
		}

		// Init memory
		vertexLastSegment.clear();
		lastVels.resize(verts.size(), vec3(0));
		xpbdSegWs.resize(segs.size(), vec3(0));
		lambdaTradGammas.resize(segs.size(), 1.f);

		dxs.resize(verts.size(), vec3(0));
		xpbdSegLambdas.resize(segs.size(), vec3(0));
		xpbdAngleLambdas.resize(restAngles.size(), vec4(0));
		xpbdLastQs.resize(segs.size(), vec3(0));

		if (glad_glGetString != nullptr) {
			vertBuffer = new ComputeBuffer(sizeof(Vert), verts.size());
			segBuffer = new ComputeBuffer(sizeof(Seg), segs.size());
			renderMesh = Kit::genCylMesh(6, 1, false);
		}

		// Enable multithreading if there are enough segments
		if (segs.size() > COSSERAT_THREAD_BATCH * 2)
			multithreaded = true;
	}

	void Sim::uploadGPU() {
		vertBuffer->upload(verts.data());
		segBuffer->upload(segs.data(), segs.size());
	}

	Sim::Sim() {
		meta.gravity = vec3(0, -9.8f, 0);
	}

	Sim::~Sim() {
		if (vertBuffer) delete vertBuffer;
		if (segBuffer) delete segBuffer;
		if (renderMesh) delete renderMesh;
	}
}
