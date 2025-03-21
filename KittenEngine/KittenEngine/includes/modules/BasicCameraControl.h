#pragma once
// Jerry Hsu 2022

#include <iostream>
#include <fstream>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

namespace Kitten {
	/// <summary>
	/// An implementaiton of basic camera controls. 
	/// Call processMousePos(), processMouseButton(), and processMouseScroll() to update.
	/// 
	/// Left mouse drag - pan/move
	/// Right mouse draw - rotate
	/// Scroll wheel - zoom
	/// 
	/// </summary>
	class BasicCameraControl {
	public:
		// Current camera position
		vec3 pos = vec3(0);

		// Current camera angle
		vec2 angle = vec2(0);

		// Current camera distance
		float distance = 2;

		// General mouse sensitivity for everything
		vec2 mouseSensitivity = vec2(0.75, 1);

		// Per axis move speed (x, y, z)
		vec3 moveSpeed = vec3(0.002f);

		// Per axis rotation speed (yaw, pitch)
		vec2 rotSpeed = vec3(0.8f);

		// Zoom speed
		float zoomSpeed = 0.2f;

		// Maximum angle allowed, in degrees
		float maxPitch = 80.f;

		// Minimum camera distance allowed.
		float minDistance = 0.1f;

		// The last mouse position recorded from 0 to 1
		vec2 lastMousePos = vec2(0);

	private:
		bvec2 buttonDown = bvec2(false);

	public:
		mat4 getViewMatrix() {
			mat4 mat(1);
			mat = glm::rotate(mat, glm::radians(angle.y), { 1.0f, 0.0f, 0.0f });
			mat = glm::rotate(mat, glm::radians(angle.x), { 0.0f, 1.0f, 0.0f });
			mat = glm::translate(mat, -pos);
			mat[3][2] -= distance;
			return mat;
		}

		void processMousePos(double xp, double yp) {
			vec2 mousePos = vec2(xp, yp);
			vec2 mouseDelta = mousePos - lastMousePos;

			if (buttonDown[0]) {
				vec3 dir = vec3(-mouseSensitivity.x * mouseDelta.x, mouseSensitivity.y * mouseDelta.y, 0) * mat3(getViewMatrix());
				pos += distance * moveSpeed * dir;
			}
			if (buttonDown[1]) {
				angle += rotSpeed * mouseSensitivity * mouseDelta;
				angle.y = glm::clamp(angle.y, -maxPitch, maxPitch);
			}

			lastMousePos = mousePos;
		}

		void processMouseButton(int button, int action, int mode) {
			if (button == GLFW_MOUSE_BUTTON_LEFT) buttonDown[0] = action == GLFW_PRESS;
			if (button == GLFW_MOUSE_BUTTON_RIGHT) buttonDown[1] = action == GLFW_PRESS;
		}

		void processMouseScroll(double xoffset, double yoffset) {
			distance = glm::max(minDistance, distance * powf(1 + zoomSpeed, (float)(yoffset + xoffset)));
		}

		void loadSettings(const char* fileName) {
			if (FILE* file = fopen(fileName, "rb")) {
				fread(this, sizeof(BasicCameraControl), 1, file);
				fclose(file);
			}
		}

		void saveSettings(const char* fileName) {
			FILE* file = fopen(fileName, "wb");
			fwrite(this, sizeof(BasicCameraControl), 1, file);
			fclose(file);
		}
	};

	/// <summary>
	/// An implementaiton of basic camera controls for 2D. 
	/// Call processMousePos(), processMouseButton(), and processMouseScroll() to update.
	/// Remember to multiply xp by aspect ratio
	/// 
	/// Left mouse drag - pan/move
	/// Scroll wheel - zoom
	/// 
	/// </summary>
	class BasicCameraControl2D {
	public:
		// Current camera position
		vec2 pos = vec3(0);

		// Current camera distance
		float distance = 2;

		// Zoom speed
		float zoomSpeed = 0.2f;

		// Minimum camera distance allowed.
		float minDistance = 0.1f;

		// The last mouse position recorded from 0 to 1
		vec2 lastMousePos = vec2(0);

	private:
		bool buttonDown = false;

	public:
		mat4 getViewMatrix() {
			return glm::translate(mat4(1), vec3(-pos, 0));
		}

		void processMousePos(double xp, double yp) {
			vec2 mousePos = vec2(xp, yp);
			vec2 mouseDelta = mousePos - lastMousePos;

			if (buttonDown) {
				vec2 dir = vec2(-mouseDelta.x, mouseDelta.y);
				pos += dir * distance * 2.f;
			}

			lastMousePos = mousePos;
		}

		void processMouseButton(int button, int action, int mode) {
			if (button == GLFW_MOUSE_BUTTON_RIGHT) buttonDown = action == GLFW_PRESS;
		}

		void processMouseScroll(double xoffset, double yoffset) {
			distance = glm::max(minDistance, distance * powf(1 + zoomSpeed, -(float)(yoffset + xoffset)));
		}

		void loadSettings(const char* fileName) {
			if (FILE* file = fopen(fileName, "rb")) {
				fread(this, sizeof(BasicCameraControl2D), 1, file);
				fclose(file);
			}
		}

		void saveSettings(const char* fileName) {
			FILE* file = fopen(fileName, "wb");
			fwrite(this, sizeof(BasicCameraControl2D), 1, file);
			fclose(file);
		}
	};
}