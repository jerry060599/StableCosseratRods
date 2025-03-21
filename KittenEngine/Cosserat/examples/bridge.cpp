#include "../CosseratExamples.h"

namespace Cosserat {
	Sim* generateLargeBridgeExample() {
		Sim* sim = new Cosserat::Sim();

		// 1:10 Scale version of the Margaret Hunt Hill Bridge in Dallas, Texas
		const int numArcSegs = 84;
		const int numSpanSegs = 72;
		const float mainSpanSegLen = 1.f;
		const int numCablesPerQuater = 14;
		const float maxCableSegLen = 0.7f;

		// Material parameters

		const float mainSpanDensity = 0.25f;
		const float arcDensity = 2.f;
		const float cableDensity = 0.01f;

		const float kssArc = 8e4;
		const float kbtArc = 4e4;
		const float kssSpan = 2e4;
		const float kbtSpan = 5;
		const float kbtTie = kbtArc;
		const float kssCable = 2e3;

		// Main span
		for (int i = 0; i <= numSpanSegs; i++)
			sim->addVertex(vec3(i * mainSpanSegLen, 2, 0));
		for (int i = 0; i < numSpanSegs; i++)
			sim->addSegment(ivec2(i, i + 1), kssSpan, mainSpanDensity);
		for (int i = 0; i < numSpanSegs - 1; i++)
			sim->addRestAngle(ivec2(i, i + 1), kbtSpan);

		// Arc
		const int arcVOffset = sim->verts.size();
		float totalArcLength = (exp(2.1f * 1.8f) - exp(-2.1f * 1.8f)) / 2.1f;
		float segArcLength = totalArcLength / numArcSegs;
		for (int i = 0; i <= numArcSegs; i++) {
			float t = 2.1f * (i * segArcLength - 0.5f * totalArcLength);
			float x = log(t + sqrt(1 + t * t)) / 2.1f;
			float y = 22.2f - 0.5f * (exp(2.1f * x) + exp(-2.1f * x));
			sim->addVertex(vec3(numSpanSegs / 2 * mainSpanSegLen, y, x));
		}

		for (int i = 0; i < numArcSegs; i++)
			sim->addSegment(ivec2(i + arcVOffset, i + 1 + arcVOffset), kssArc, arcDensity);
		for (int i = 0; i < numArcSegs - 1; i++)
			sim->addRestAngle(ivec2(arcVOffset - 1 + i, arcVOffset + i), kbtArc);

		sim->pinVertex(0);
		sim->pinVertex(numSpanSegs);

		sim->pinVertex(arcVOffset);
		sim->pinVertex(arcVOffset + numArcSegs);
		sim->pinSegment(arcVOffset - 1);
		sim->pinSegment(arcVOffset + numArcSegs - 2);

		// Side spans
		// Find the first arc vertex above the main span
		vec3 sideSpanPos(0, 0.3, 1);
		int arcSpanVertex = arcVOffset;
		for (; arcSpanVertex < sim->verts.size(); arcSpanVertex++)
			if (sim->verts[arcSpanVertex].pos.y > 2) {
				sideSpanPos = sim->verts[arcSpanVertex].pos;
				break;
			}

		// Build side span
		sideSpanPos.x = 0;
		int numHalfSpanSegs = numSpanSegs / 2 - 1;

		// Build side spans as quaters connecting at the arc
		// One side of the span
		int halfSpanVOffset = sim->verts.size();
		int halfSpanSOffset = sim->segs.size();
		for (int i = 0; i <= numSpanSegs; i++)
			if (i != numSpanSegs / 2)
				sim->addVertex(vec3(i * mainSpanSegLen, sideSpanPos.y, sideSpanPos.z));

		for (int i = 0; i < numSpanSegs; i++) {
			ivec2 seg(i, i + 1);
			if (seg.x == numSpanSegs / 2) seg.x = arcSpanVertex - halfSpanVOffset;
			else if (seg.x > numSpanSegs / 2) seg.x--;
			if (seg.y == numSpanSegs / 2) seg.y = arcSpanVertex - halfSpanVOffset;
			else if (seg.y > numSpanSegs / 2) seg.y--;
			sim->addSegment(seg + ivec2(halfSpanVOffset), kssSpan, mainSpanDensity);
		}
		for (int i = 0; i < numSpanSegs - 1; i++)
			sim->addRestAngle(ivec2(i + halfSpanSOffset, i + halfSpanSOffset + 1), kbtSpan);

		sim->pinVertex(halfSpanVOffset);
		sim->pinVertex(halfSpanVOffset + numSpanSegs - 1);

		// Other side of the span
		halfSpanVOffset = sim->verts.size();
		halfSpanSOffset = sim->segs.size();
		for (int i = 0; i <= numSpanSegs; i++)
			if (i != numSpanSegs / 2)
				sim->addVertex(vec3(i * mainSpanSegLen, sideSpanPos.y, -sideSpanPos.z));

		for (int i = 0; i < numSpanSegs; i++) {
			ivec2 seg(i, i + 1);
			if (seg.x == numSpanSegs / 2) seg.x = numArcSegs - arcSpanVertex + 2 * arcVOffset - halfSpanVOffset;
			else if (seg.x > numSpanSegs / 2) seg.x--;
			if (seg.y == numSpanSegs / 2) seg.y = numArcSegs - arcSpanVertex + 2 * arcVOffset - halfSpanVOffset;
			else if (seg.y > numSpanSegs / 2) seg.y--;
			sim->addSegment(seg + ivec2(halfSpanVOffset), kssSpan, mainSpanDensity);
		}
		for (int i = 0; i < numSpanSegs - 1; i++)
			sim->addRestAngle(ivec2(i + halfSpanSOffset, i + halfSpanSOffset + 1), kbtSpan);
		sim->pinVertex(halfSpanVOffset);
		sim->pinVertex(halfSpanVOffset + numSpanSegs - 1);

		// Connect the side spans with the center span
		for (int i = 4; i < numSpanSegs - 1; i += 4) {
			int left = arcVOffset + numArcSegs + i + 1;
			int right = arcVOffset + numArcSegs + numSpanSegs + i + 1;

			if (i == numSpanSegs / 2) {
				left = arcSpanVertex;
				right = numArcSegs - arcSpanVertex + 2 * arcVOffset;
			}
			else if (i > numSpanSegs / 2) {
				left--;
				right--;
			}

			sim->addSegment(ivec2(i, left), kssSpan, mainSpanDensity);
			sim->addSegment(ivec2(i, right), kssSpan, mainSpanDensity);
			int tie = sim->segs.size() - 2;
			sim->addRestAngle(ivec2(tie, tie + 1), kbtTie);
			sim->addRestAngle(ivec2(i, tie), kbtTie);
			sim->addRestAngle(ivec2(i, tie + 1), kbtTie);
		}

		// Connect the cables
		auto connectCable = [&](int i, int j) {
			vec3 a = sim->verts[i].pos;
			vec3 b = sim->verts[j].pos;
			float len = length(b - a);
			int numSegs = (int)ceil(len / maxCableSegLen);

			const int cableVOffset = sim->verts.size();
			for (int i = 1; i < numSegs; i++)
				sim->addVertex(a + (i / (float)numSegs) * (b - a));
			sim->addSegment(ivec2(i, cableVOffset), kssCable, cableDensity);
			for (int i = 0; i < numSegs - 2; i++)
				sim->addSegment(ivec2(cableVOffset + i, cableVOffset + i + 1), kssCable, cableDensity);
			sim->addSegment(ivec2(cableVOffset + numSegs - 2, j), kssCable, cableDensity);
			};

		for (size_t i = 0; i < numCablesPerQuater; i++) {
			connectCable(numSpanSegs / 2 - 2 - 2 * i, arcVOffset + numArcSegs / 2 - i - 1);
			connectCable(numSpanSegs / 2 + 3 + 2 * i, arcVOffset + numArcSegs / 2 - i - 1);
		}

		connectCable(numSpanSegs / 2 - 1, arcVOffset + numArcSegs / 2);
		connectCable(numSpanSegs / 2 + 1, arcVOffset + numArcSegs / 2);

		for (size_t i = 0; i < numCablesPerQuater; i++) {
			connectCable(numSpanSegs / 2 - 3 - 2 * i, arcVOffset + numArcSegs / 2 + i + 1);
			connectCable(numSpanSegs / 2 + 2 + 2 * i, arcVOffset + numArcSegs / 2 + i + 1);
		}

		sim->renderRadius = 0.05f;

		sim->maxH = 1e-3;
		sim->meta.numItr = 4;
		sim->meta.drag = 0.1f;
		sim->meta.damping = 5e-8;
		sim->meta.gravity = vec3(0, -0.981f, 0);	// This is 1:10 scale

		sim->renderByDensity = 1;
		sim->renderDensity = 10;
		sim->meshHeightLen = 0.1f;

		// 0.5% strain threshold
		// Real steel should have a strain threshold of 0.2% ish
		sim->meta.segStrainThreshold = 0.005f;
		sim->meshExcludeInlinePhong = 1;

		return sim;
	}

	Sim* generateSmallBridgeExample() {
		Sim* sim = new Cosserat::Sim();

		// A model of a simple foot bridge modeled after the sundial bridge in Redding, California
		const int numSpanSegs = 34;
		const int numPillarSegs = 20;
		const int numCables = 13;
		const float segLen = 0.5f;

		const float spanDensity = 0.1f;
		const float pillarDensity = 1.0f;
		const float cableDensity = 0.01f;

		const float kssPillar = 1e5;
		const float kbtPillar = 1e4;
		const float kssSpan = 1e3;
		const float kbtSpan = 1e1;
		const float kssCable = 1e3;

		const vec3 pillarPos(2, 0, 2);
		const vec3 pillarDir = normalize(vec3(2, 1.4f, -1));

		// Main span
		for (int i = 0; i <= numSpanSegs; i++)
			sim->addVertex(vec3(i * segLen, 1, 0));
		for (int i = 0; i < numSpanSegs; i++)
			sim->addSegment(ivec2(i, i + 1), kssSpan, spanDensity);
		for (int i = 0; i < numSpanSegs - 1; i++)
			sim->addRestAngle(ivec2(i, i + 1), kbtSpan);

		// Pillars
		const int pillarVOffset = sim->verts.size();
		const int pillarSOffset = sim->segs.size();
		for (int i = 0; i <= numPillarSegs; i++)
			sim->addVertex(pillarPos + i * segLen * pillarDir);
		for (int i = 0; i < numPillarSegs; i++)
			sim->addSegment(ivec2(i, i + 1) + ivec2(pillarVOffset), kssPillar, pillarDensity);
		for (int i = 0; i < numPillarSegs - 1; i++)
			sim->addRestAngle(ivec2(i, i + 1) + ivec2(pillarSOffset), kbtPillar);

		for (int i = 0; i <= numPillarSegs; i++) {
			vec3 p = pillarPos + i * segLen * pillarDir;
			p.x = +numSpanSegs * segLen - p.x;
			p.z *= -1;
			sim->addVertex(p);
		}
		for (int i = 0; i < numPillarSegs; i++)
			sim->addSegment(ivec2(i, i + 1) + ivec2(pillarVOffset + numPillarSegs + 1), kssPillar, pillarDensity);
		for (int i = 0; i < numPillarSegs - 1; i++)
			sim->addRestAngle(ivec2(i, i + 1) + ivec2(pillarSOffset + numPillarSegs), kbtPillar);

		// Connect the cables
		auto connectCable = [&](int i, int j) {
			vec3 a = sim->verts[i].pos;
			vec3 b = sim->verts[j].pos;
			float len = length(b - a);
			int numSegs = (int)ceil(len / segLen);
			numSegs = glm::max(numSegs, 3);

			const int cableVOffset = sim->verts.size();
			for (int i = 1; i < numSegs; i++)
				sim->addVertex(a + (i / (float)numSegs) * (b - a));
			sim->addSegment(ivec2(i, cableVOffset), kssCable, cableDensity);
			for (int i = 0; i < numSegs - 2; i++)
				sim->addSegment(ivec2(cableVOffset + i, cableVOffset + i + 1), kssCable, cableDensity);
			sim->addSegment(ivec2(cableVOffset + numSegs - 2, j), kssCable, cableDensity);
			};

		// Pin

		sim->pinVertex(0);
		sim->pinVertex(numSpanSegs);

		sim->pinVertex(pillarVOffset);
		sim->pinVertex(pillarVOffset + numPillarSegs + 1);
		sim->pinSegment(pillarSOffset);
		sim->pinSegment(pillarSOffset + numPillarSegs);

		for (int i = 0; i < numCables; i++)
			connectCable(i + 8, pillarVOffset + i + numPillarSegs - numCables - 1);
		for (int i = 0; i < numCables; i++)
			connectCable(numSpanSegs - (i + 8) - 1, pillarVOffset + i + 2 * numPillarSegs - numCables);

		sim->maxH = 1e-3;
		sim->meta.numItr = 4;
		sim->meta.drag = 0.1f;
		sim->meta.damping = 1e-7;
		sim->meta.gravity = vec3(0, -0.981f, 0);	// This is 1:10 scale

		sim->renderByDensity = 1;
		sim->renderDensity = 30;

		sim->meta.segStrainThreshold = 0.005f;

		return sim;
	}
}