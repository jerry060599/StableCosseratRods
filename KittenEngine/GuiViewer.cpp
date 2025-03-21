
#include <cassert>
#include <iostream>

#include "Cosserat/Cosserat.h"
#include "Cosserat/CosseratExamples.h"
#include "KittenEngine/includes/KittenEngine.h"
#include "KittenEngine/includes/modules/BasicCameraControl.h"
#include "KittenEngine/includes/modules/StopWatch.h"

using namespace glm;
using namespace std;

Kit::BasicCameraControl camera;
Cosserat::Sim* sim = nullptr;

bool simulate = false;
float timeScale = 1.;
Kit::Dist simSpeedDist;

const float EXPORT_DT = 1 / 30.f;
bool exportSim = false;
bool exportMesh = true;
int exportFrame = 0;

void setBridgeScenario(float cat);
void setTreeScenario();

void renderScene() {
	if (simulate) {
		// Dynamic dt
		float advTime = EXPORT_DT;
		const float realTime = ImGui::GetIO().DeltaTime * timeScale;
		if (!exportSim)
			advTime = glm::min(realTime, 1 / 40.f);

		Kit::StopWatch timer;
		sim->advance(advTime);
		float measuredTime = timer.time();

		float ss = advTime / measuredTime;
		simSpeedDist.accu(ss);

		if (exportSim) {
			char* filename = new char[100];

			if (true)
				if (exportMesh) {
					sprintf(filename, "./frames/frame_%03d.obj", exportFrame);
					sim->exportMesh(filename);
				}
				else {
					sprintf(filename, "./frames/frame_%03d.crs", exportFrame);
					sim->saveToBinary(filename, false);
				}

			if (exportFrame == 599) {
				exportSim = false;
				simulate = false;
			}

			exportFrame++;
		}
	}

	Kit::projMat = glm::perspective(45.0f, Kit::getAspect(), 0.005f, 512.f);
	Kit::viewMat = camera.getViewMatrix();

	// Render everything
	Kit::startRender();
	glClearColor(0.f, 0.f, 0.f, 1.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	sim->render();
}

void setBridgeScenario(float speed) {
	// Cat5 hurricane has wind speed of 70m/s.
	// At our scale, thats 7unit/s
	const float waveSpeed = 0.1f * speed;
	const float waveLen = 4.f;
	const float dragCoeff = 7e-3;

	const vec3 windDir = normalize(vec3(2, 0, 3));
	const float r = 0.7f;				// Noise ratio
	const float baseWindLoad = 1.f;		// Minimum wind load

	const float freq = 2 * glm::pi<float>() / waveLen;
	const float animationSpeed = 2 * glm::pi<float>() * waveSpeed / waveLen;
	const float maxSum = 1 / (1 - r);

	sim->setWindFieldExternalForce([=](vec3 pos, float time) -> vec3 {
		float phase = freq * dot(windDir, pos) - animationSpeed * time;
		float windMag = 0;
		float alpha = 1;
		for (int i = 1; i <= 3; i++) {
			windMag += alpha * sin(alpha * phase);
			alpha *= r;
		}

		windMag = (windMag + maxSum + baseWindLoad) / (baseWindLoad + 2 * maxSum);
		return (waveSpeed * windMag) * windDir;
		}, dragCoeff);
}

void setTreeScenario() {
	const float waveSpeed = 15;
	const float waveLen = 20.f;
	const float dragCoeff = 3;

	const vec3 windDir = normalize(vec3(2, 0, 3));
	const float r = 0.7f;				// Noise ratio
	const float baseWindLoad = 0.f;		// Minimum wind load

	const float freq = 2 * glm::pi<float>() / waveLen;
	const float animationSpeed = 2 * glm::pi<float>() * waveSpeed / waveLen;
	const float maxSum = 1 / (1 - r);

	sim->setWindFieldExternalForce([=](vec3 pos, float time) -> vec3 {
		float phase = freq * dot(windDir, pos) - animationSpeed * time;
		float windMag = 0;
		float alpha = 1;
		for (int i = 1; i <= 3; i++) {
			windMag += alpha * sin(alpha * phase);
			alpha *= r;
		}

		windMag = (windMag + maxSum + baseWindLoad) / (baseWindLoad + 2 * maxSum);
		return (waveSpeed * windMag) * windDir;
		}, dragCoeff);
}

void switchSim(Cosserat::Sim* newSim) {
	delete sim;
	sim = newSim;
	sim->init();
	sim->performSelfCheck();
	camera.pos = sim->verts[0].pos;
	exportFrame = 0;
	simSpeedDist = Kit::Dist();
}

void renderGui() {
	ImGui::Begin("Control Panel");
	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	ImGui::Text("Measured sim/wall time ratio avg %.3f, SD: %.3f, N=%d", simSpeedDist.mean(), simSpeedDist.sd(), simSpeedDist.num);
	ImGui::Text("Sim time: %.3f", sim->meta.time);
	ImGui::Text("%d verts, %d segs, and %d angles.", sim->verts.size(), sim->segs.size(), sim->restAngles.size());

	ImGui::Separator();
	if (ImGui::Button("Load Cantilever"))
		switchSim(Cosserat::generateCantileverExample());
	ImGui::SameLine();
	if (ImGui::Button("Load VBD Failure Case"))
		switchSim(Cosserat::generateVBDFailureExample());
	ImGui::SameLine();
	if (ImGui::Button("Load Slinky"))
		switchSim(Cosserat::generateSlinkyExample());

	if (ImGui::Button("Load 1 Tree"))
		switchSim(Cosserat::generateTreeExample(1));

	ImGui::SameLine();
	if (ImGui::Button("Load 3 Trees"))
		switchSim(Cosserat::generateTreeExample(3));

	ImGui::SameLine();
	if (ImGui::Button("Load 32 Trees"))
		switchSim(Cosserat::generateTreeExample(32));

	if (ImGui::Button("Load Big Bridge"))
		switchSim(Cosserat::generateLargeBridgeExample());

	ImGui::SameLine();
	if (ImGui::Button("Load Small Bridge"))
		switchSim(Cosserat::generateSmallBridgeExample());

	if (ImGui::Button("Load Slingshot"))
		switchSim(Cosserat::generateSlingshotExample());

	ImGui::Separator();

	const char* methodNames[]{ "Stable Cosserat Rods", "Stable Cosserat Exact",
			"VBD","VBD w/ Linearized Constraints", "VBD w/ LC & PSD Projection", "XPBD" };

	if (ImGui::TreeNode("Simulation")) {
		ImGui::Checkbox("Simulate", &simulate);
		ImGui::SliderFloat("Time Scale", &timeScale, 0.001, 2);

		ImGui::Separator();

		ImGui::Combo("Method", (int*)&sim->method, methodNames, 6);

		ImGui::DragFloat3("Gravity", (float*)&sim->meta.gravity);
		ImGui::InputFloat("Max dt", &sim->maxH, 1e-4, 1e-3, "%.5f");
		ImGui::InputInt("Itr", &sim->meta.numItr);

		ImGui::Separator();

		ImGui::SliderFloat("Drag", &sim->meta.drag, 0, 1, "%.3f");
		float damping = log10f(sim->meta.damping);
		if (ImGui::SliderFloat("Damping", &damping, -10, -0.1, "%.3f"))
			sim->meta.damping = powf(10, damping);

		ImGui::Separator();

		bool useSegStrain = glm::isfinite(sim->meta.segStrainThreshold);
		if (ImGui::Checkbox("Use Segment Strain Threshold", &useSegStrain))
			sim->meta.segStrainThreshold = useSegStrain ? 0.1f : INFINITY;
		if (useSegStrain)
			ImGui::SliderFloat("Seg Strain Threshold", &sim->meta.segStrainThreshold, 0, 1);
		else
			sim->meta.segStrainThreshold = INFINITY;

		bool useAngleStrain = glm::isfinite(sim->meta.angleStrainThreshold);
		if (ImGui::Checkbox("Use Angle Strain Threshold", &useAngleStrain))
			sim->meta.angleStrainThreshold = useAngleStrain ? 0.1f : INFINITY;
		if (useAngleStrain)
			ImGui::SliderFloat("Angle Strain Threshold", &sim->meta.angleStrainThreshold, 0, 1);
		else
			sim->meta.angleStrainThreshold = INFINITY;

		ImGui::Separator();

		ImGui::Checkbox("Multithread", &sim->multithreaded);

		ImGui::TreePop();
	}

	if (ImGui::TreeNode("Rendering")) {
		ImGui::SliderFloat("Render radius", &sim->renderRadius, 0.001, 1.f);
		ImGui::SliderFloat("Render density", &sim->renderDensity, 0.001, 100000.f);
		ImGui::SliderFloat("Render by density", &sim->renderByDensity, 0.f, 1.f);
		ImGui::TreePop();
	}

	if (ImGui::TreeNode("Export")) {
		ImGui::Checkbox("Export", &exportSim);
		ImGui::Checkbox("Export as Mesh", &exportMesh);
		static bool saveOnlyGeometry = false;
		ImGui::SameLine();
		ImGui::Checkbox("Save only geometry", &saveOnlyGeometry);

		if (ImGui::Button("Export mesh now")) {
			sim->exportMesh("./frame.obj");
			printf("Exported mesh to ./frame.obj\n");
		}
		ImGui::SameLine();
		if (ImGui::Button("Export sim now")) {
			sim->saveToBinary("./frame.crs", saveOnlyGeometry);
			printf("Exported sim to ./frame.crs\n");
		}

		if (ImGui::Button("Restore sim from frame")) {
			delete sim;
			sim = new Cosserat::Sim();
			sim->loadFromBinary("./frame.crs", saveOnlyGeometry);
			sim->init();
			printf("Restored sim from ./frame.crs\n");
		}
		ImGui::Separator();

		if (ImGui::Button("Set gravity scenario")) {
			// Start increasing gravity after 3 seconds
			sim->externalForce = [](Cosserat::Sim*, Cosserat::Vert& vert, float time) -> vec3 {
				float scale = glm::clamp(1 * (time - 3), 0.f, 5.f);
				return sim->meta.gravity * (scale / vert.invMass);
				};
		}

		static int hurricaneCat = 5;
		if (ImGui::Button("Set bridge scenario"))
			setBridgeScenario(70 + 15 * (hurricaneCat - 5));
		ImGui::SameLine();
		ImGui::InputInt("Bridge Hurricane Cat", &hurricaneCat);

		if (ImGui::Button("Set tree scenario"))
			setTreeScenario();

		if (ImGui::Button("Clear scenario"))
			sim->externalForce = nullptr;

		ImGui::TreePop();
	}

	if (ImGui::TreeNode("Debug")) {
		if (ImGui::Button("Zero velocity"))
			for (auto& v : sim->vels) v = vec3(0);
		if (ImGui::Button("Perform Self Check"))
			sim->performSelfCheck();

		if (ImGui::Button("Scramble")) {
			for (int i = 0; i < sim->verts.size(); i++)
				if (sim->verts[i].invMass > 0)
					sim->verts[i].pos += 2.f * (vec3(rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, rand() / (float)RAND_MAX) * 2.f - 1.f);
		}

		if (ImGui::Button("Test convergence")) {
			switchSim(Cosserat::generateLargeBridgeExample());
			sim->method = Cosserat::IntegrationMethod::StableCosserat;
			sim->checkSegConvergence();

			switchSim(Cosserat::generateLargeBridgeExample());
			sim->method = Cosserat::IntegrationMethod::StableCosseratExact;
			sim->checkSegConvergence();

			switchSim(Cosserat::generateLargeBridgeExample());
			sim->method = Cosserat::IntegrationMethod::VBD;
			sim->checkSegConvergence();

			switchSim(Cosserat::generateLargeBridgeExample());
			sim->method = Cosserat::IntegrationMethod::LCVBD;
			sim->checkSegConvergence();

			switchSim(Cosserat::generateLargeBridgeExample());
			sim->method = Cosserat::IntegrationMethod::PSDPLCVBD;
			sim->checkSegConvergence();

			printf("\n\nDone.\n");
		}

		if (ImGui::Button("Test convergence over time"))
			sim->checkSegConvergenceOverTime();
		if (ImGui::Button("Check gamma"))
			sim->checkExactLambdaConvergence();

		if (ImGui::Button("Perform Breaking Test")) {
			sim->deformDestroy(
				[&](Cosserat::Sim* sim, Cosserat::Seg& seg, vec3 strain)->bool {
					return true;
				},
				[&](Cosserat::Sim* sim, Cosserat::RestAngle& angle, vec4 strain)->bool {
					return true;
				}
			);
			sim->performSelfCheck();
		}


		if (ImGui::Button("Print stiffness range")) {
			float minK = INFINITY, maxK = -INFINITY;
			for (auto& seg : sim->segs) {
				minK = glm::min(minK, abs(seg.kss));
				maxK = glm::max(maxK, abs(seg.kss));
			}
			printf("Min kss: %f, Max kss: %f\n", minK, maxK);

			minK = INFINITY, maxK = -INFINITY;
			for (auto& angle : sim->restAngles) {
				minK = glm::min(minK, angle.kbt);
				maxK = glm::max(maxK, angle.kbt);
			}
			printf("Min kbt: %f, Max kbt: %f\n", minK, maxK);

			minK = INFINITY, maxK = -INFINITY;
			for (auto& vert : sim->verts)
				if (vert.invMass > 0) {
					minK = glm::min(minK, 1 / vert.invMass);
					maxK = glm::max(maxK, 1 / vert.invMass);
				}
			printf("Min mass: %f, Max mass: %f\n", minK, maxK);
		}


		if (ImGui::Button("Time simulation")) {
			const float warmupTime = 5;
			const float advTime = 15;
			printf("\nTimed simulation:\n");
			printf("%d verts, %d segs, and %d angles.\n", sim->verts.size(), sim->segs.size(), sim->restAngles.size());
			printf("Max DT: %f\n", sim->maxH);
			printf("Num itr: %d\n", sim->meta.numItr);
			printf("Method: %s\n", methodNames[(int)sim->method]);
			printf("Multithreaded: %s\n", sim->multithreaded ? "True" : "False");

			printf("Warming up for %.1f sec...\n", warmupTime);
			Kit::StopWatch warmupTimer;
			while (true) {
				sim->step(sim->maxH);
				warmupTimer.time();
				if (warmupTimer.totTime() > warmupTime) break;
			}

			printf("Advancing simulation by %.3f sec...\n", advTime);
			Kit::StopWatch timer;

			sim->advance(advTime);

			double realtime = timer.time();
			printf("Done. Sim took %.3f sec. Ratio %.4fx realtime. %.2f ms per 33.3ms (30fps) frame.\n", realtime, advTime / realtime, realtime / (30.f * advTime) * 1000);
		}

		if (ImGui::Button("Print gammas")) {
			for (size_t i = 0; i < sim->lambdaTradGammas.size(); i++)
				printf("Seg %d: %f\n", i, sim->lambdaTradGammas[i]);
		}

		ImGui::TreePop();
	}

	ImGui::End();
}

void initScene() {
	Kit::loadDirectory("resources");

	Kit::UBOLight light;
	light.col = vec4(1, 1, 1, 1);
	light.dir = normalize(vec3(1, -2, -1));
	light.hasShadow = false;
	light.type = (int)Kit::KittenLight::DIR;
	Kit::lights.push_back(light);
	light.shadowBias = 0.0001f;
	Kit::shadowDist = 0.5f;

	Kit::ambientLight.col = vec4(0);

	sim = Cosserat::generateTreeExample(1);

	sim->init();
	sim->performSelfCheck();

	if (!sim) exit(-1);

	camera.pos = sim->verts[0].pos;
	camera.minDistance = 0.01f;
}

void mouseButtonCallback(GLFWwindow* w, int button, int action, int mode) {
	camera.processMouseButton(button, action, mode);
}

void cursorPosCallback(GLFWwindow* w, double xp, double yp) {
	camera.processMousePos(xp, yp);
}

void scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
	camera.processMouseScroll(xoffset, yoffset);
}

void keyCallback(GLFWwindow* w, int key, int scancode, int action, int mode) {
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(w, true);

	if (key == GLFW_KEY_F && action == GLFW_PRESS)
		sim->step(sim->maxH);

	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
		simulate = !simulate;
}

int main(int argc, char** argv) {
	// Init window and OpenGL
	Kit::initWindow(ivec2(1920, 1080), "OpenGL Window");

	// Register callbacks
	Kit::getIO().mouseButtonCallback = mouseButtonCallback;
	Kit::getIO().cursorPosCallback = cursorPosCallback;
	Kit::getIO().scrollCallback = scrollCallback;
	Kit::getIO().keyCallback = keyCallback;

	// Init scene
	initScene();

	while (!Kit::shouldClose()) {
		Kit::startFrame();
		renderScene();		// Render
		renderGui();		// GUI Render
		Kit::endFrame();
	}

	Kit::terminate();
	return 0;
}