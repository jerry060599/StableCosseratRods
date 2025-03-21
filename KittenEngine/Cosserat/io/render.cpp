#include "../Cosserat.h"

namespace Cosserat {
	void Sim::render() {
		uploadGPU();

		static auto base = Kit::get<Kit::Shader>("resources\\shaders\\rodBase.glsl");
		static auto forward = Kit::get<Kit::Shader>("resources\\shaders\\rodForward.glsl");

		vertBuffer->bind(5);
		segBuffer->bind(6);

		base->setFloat("radius", renderRadius);
		forward->setFloat("radius", renderRadius);
		forward->setFloat("density", renderDensity);
		forward->setFloat("density", renderDensity);
		forward->setFloat("renderByDensity", renderByDensity);
		forward->setFloat("renderByDensity", renderByDensity);

		Kit::renderInstancedForward(renderMesh, segs.size(), base, forward);
	}
}