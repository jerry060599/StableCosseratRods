#pragma once
// Jerry Hsu, 2024

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <functional>

namespace Kitten {
	using namespace glm;

	// Gizmos can be queued with queueGizmo___ functions.
	// Gizmos provide a quick way to visualize/manipulate data in the scene.
	// Think of it as ImGui for the 3D renderer.

	void queueGizmoTranslate(const char* label, vec3* ptr, bool interactive = true);

	void queueGizmoLine(vec3 a, vec3 b, float thickness, vec4 color);
	void queueGizmoLineScreenspace(vec2 a, vec2 b, float thickness, vec4 color);

	void queueGizmoSquare(vec3 a, float size, vec4 color);
	void queueGizmoSquare(mat4 modelMat, vec4 color);
	void queueGizmoSquareScreenspace(vec2 a, float size, vec4 color);

	// Do not call these functions below directly.

	typedef struct GizmoData {
		mat4 model;
		vec4 color;
	};

	struct GizmoCallback {
		std::function<void(vec2, float)> onDrag = nullptr;
		std::function<void(float)> onMouseOver = nullptr;
		std::function<void(float)> onClick = nullptr;
	};

	void addGizmoLineCallback(vec3 a, vec3 b, float thickness, void* ptr, GizmoCallback callback);
	void addGizmoCircleCallback(vec3 a, float radius, void* ptr, GizmoCallback callback);
	void addGizmoSquareCallback(mat4 mat, void* ptr, GizmoCallback callback);

	void recordCallback(float dist, void* ptr, GizmoCallback callback);
	void processCallback(float dist, void* ptr, GizmoCallback callback);

	void queueGizmo(GizmoData mat);

	void initGizmos();
	void renderGizmos();
	void gizmoMouseButtonCallback(int button, int action, int mode);
	void gizmoCursorPosCallback(double xp, double yp);
}