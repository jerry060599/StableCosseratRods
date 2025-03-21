#include "../includes/modules/KittenInit.h"
#include "../includes/modules/KittenRendering.h"
#include "../includes/modules/Gizmos.h"

#include <iostream>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

namespace Kitten {
	GLFWwindow* window = nullptr;
	GLFWCallbacks glfwCallbacks = {};
	GLFWCallbacks imGuiCallbacks = {};

	GLFWCallbacks& getIO() {
		return glfwCallbacks;
	}

	void mouseButtonCallback(GLFWwindow* w, int button, int action, int mode) {
		if (imGuiCallbacks.mouseButtonCallback)
			imGuiCallbacks.mouseButtonCallback(w, button, action, mode);
		gizmoMouseButtonCallback(button, action, mode);
		if (ImGui::GetIO().WantCaptureMouse) return;
		if (glfwCallbacks.mouseButtonCallback)
			glfwCallbacks.mouseButtonCallback(w, button, action, mode);
	}

	void cursorPosCallback(GLFWwindow* w, double xp, double yp) {
		if (imGuiCallbacks.cursorPosCallback)
			imGuiCallbacks.cursorPosCallback(w, xp, yp);
		gizmoCursorPosCallback(xp, yp);
		if (ImGui::GetIO().WantCaptureMouse) return;
		if (glfwCallbacks.cursorPosCallback)
			glfwCallbacks.cursorPosCallback(w, xp, yp);
	}

	void scrollCallback(GLFWwindow* w, double xoffset, double yoffset) {
		if (imGuiCallbacks.scrollCallback)
			imGuiCallbacks.scrollCallback(w, xoffset, yoffset);
		if (ImGui::GetIO().WantCaptureMouse) return;
		if (glfwCallbacks.scrollCallback)
			glfwCallbacks.scrollCallback(w, xoffset, yoffset);
	}

	void framebufferSizeCallback(GLFWwindow* w, int width, int height) {
		windowRes.x = width;
		windowRes.y = height;
		glViewport(0, 0, width, height);

		if (imGuiCallbacks.framebufferSizeCallback)
			imGuiCallbacks.framebufferSizeCallback(w, width, height);
		if (glfwCallbacks.framebufferSizeCallback)
			glfwCallbacks.framebufferSizeCallback(w, width, height);
	}

	void keyCallback(GLFWwindow* w, int key, int scancode, int action, int mode) {
		if (imGuiCallbacks.keyCallback)
			imGuiCallbacks.keyCallback(w, key, scancode, action, mode);
		if (glfwCallbacks.keyCallback)
			glfwCallbacks.keyCallback(w, key, scancode, action, mode);
	}

	void initWindow(glm::ivec2 res, const char* title) {
		windowRes = res;
		glfwInit();
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
		glfwWindowHint(GLFW_SAMPLES, 16);
		window = glfwCreateWindow(res.x, res.y, title, nullptr, nullptr);
		if (!window) {
			std::cerr << "Cannot create window";
			std::exit(1);
		}
		glfwMakeContextCurrent(window);

		assert(window);
		if (gladLoadGLLoader((GLADloadproc)(glfwGetProcAddress)) == 0) {
			std::cerr << "Failed to intialize OpenGL loader" << std::endl;
			std::exit(1);
		}
		assert(glGetError() == GL_NO_ERROR);

		glEnable(GL_MULTISAMPLE);

		IMGUI_CHECKVERSION();
		ImGui::CreateContext();
		ImGui::StyleColorsDark();

		if (!ImGui_ImplGlfw_InitForOpenGL(window, true))
			printf("error: ImGui GLFW init failure.");

		const char* glsl_version = "#version 130";
		if (!ImGui_ImplOpenGL3_Init(glsl_version))
			printf("error: ImGui OpenGL init failure.");

		imGuiCallbacks.mouseButtonCallback = glfwSetMouseButtonCallback(window, mouseButtonCallback);
		imGuiCallbacks.cursorPosCallback = glfwSetCursorPosCallback(window, cursorPosCallback);
		imGuiCallbacks.framebufferSizeCallback = glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);
		imGuiCallbacks.keyCallback = glfwSetKeyCallback(window, keyCallback);
		imGuiCallbacks.scrollCallback = glfwSetScrollCallback(window, scrollCallback);

		initFreetype();
		initRender();
		initGizmos();
	}

	void terminate() {
		ImGui_ImplOpenGL3_Shutdown();
		ImGui_ImplGlfw_Shutdown();
		ImGui::DestroyContext();

		glfwDestroyWindow(window);
		glfwTerminate();
	}
}