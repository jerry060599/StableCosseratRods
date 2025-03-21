#pragma once
// Jerry Hsu, 2021

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <imgui.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "Font.h"

namespace Kitten {
	typedef struct {
		GLFWcursorposfun cursorPosCallback;
		GLFWframebuffersizefun framebufferSizeCallback;
		GLFWmousebuttonfun mouseButtonCallback;
		GLFWscrollfun scrollCallback;
		GLFWkeyfun keyCallback;
	} GLFWCallbacks;

	extern GLFWwindow* window;

	GLFWCallbacks& getIO();
	void initWindow(glm::ivec2 res, const char* title = "OpenGL Window");
	void terminate();
}