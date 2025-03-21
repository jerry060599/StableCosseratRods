#pragma once

/*

// KittenEngine Quick Start Skeleton Code

#include "KittenEngine/includes/KittenEngine.h"
#include "KittenEngine/includes/modules/BasicCameraControl.h"

using namespace glm;
using namespace std;

Kit::BasicCameraControl camera;

void renderScene() {
	Kit::startRender();
	Kit::projMat = glm::perspective(45.0f, Kit::getAspect(), 0.05f, 512.f);
	Kit::viewMat = camera.getViewMatrix();

	// Do rendering here
}

void renderGui() {
	ImGui::Begin("ImGui Window");
	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

	// Do ImGui here

	ImGui::End();
}

void initScene() {
	Kit::loadDirectory("resources");

	// Do scene preprocessing here
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

void framebufferSizeCallback(GLFWwindow* w, int width, int height) {
	// Frame buffer here
}

void keyCallback(GLFWwindow* w, int key, int scancode, int action, int mode) {
	// Key inputs here
}

int main(int argc, char** argv) {
	// Init window and OpenGL
	Kit::initWindow(ivec2(800, 600), "OpenGL Window");

	// Register callbacks
	Kit::getIO().mouseButtonCallback = mouseButtonCallback;
	Kit::getIO().cursorPosCallback = cursorPosCallback;
	Kit::getIO().scrollCallback = scrollCallback;
	Kit::getIO().framebufferSizeCallback = framebufferSizeCallback;
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

*/

#include "modules/Common.h"
#include "modules/KittenAssets.h"
#include "modules/KittenRendering.h"

#include "modules/KittenInit.h"

namespace Kit = Kitten;
