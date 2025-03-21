#include "../includes/modules/Gizmos.h"
#include "../includes/modules/ComputeBuffer.h"
#include "../includes/modules/KittenAssets.h"
#include "../includes/modules/KittenRendering.h"
#include "../includes/modules/Shader.h"
#include <vector>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

namespace Kitten {
	std::vector<GizmoData> gizmoRenderQueue;
	ComputeBuffer* gizmoBuffer;
	constexpr int MAX_GIZMO_PER_DRAW = 1024;

	vec3 mouseRayOrigin = vec3(0);
	vec3 mouseRayDir = vec3(0, 0, -1);
	vec2 mouseNDC = vec2(0);

	vec2 ndcDragAmount = vec2(0);
	bool isDragging = false;
	bool isClicking = false;
	bool isLifting = false;

	float closestDist = INFINITY;
	void* closestPtr = nullptr;
	GizmoCallback closestCallback = {};

	void* pendingPtr = nullptr;
	GizmoCallback pendingCallback = {};
	bool onDragPending = false;
	bool onMouseOverPending = false;
	bool onClickPending = false;

	inline vec3 cameraPos() {
		return -transpose((mat3)viewMat) * (vec3)viewMat[3];
	}

	inline void translateGizmoAxis(vec3* ptr, vec3 axis, vec4 color, float scale) {
		float thickness = 0.0025f * scale;
		queueGizmoLine(*ptr, *ptr + (0.2f * scale) * axis, thickness, color);
		addGizmoLineCallback(*ptr, *ptr + (0.2f * scale) * axis, thickness * 2, ptr, {
			[=](vec2 drag, float dist) {
				vec3 ori, dir;
				ndcToWorldRay(mouseNDC + drag, ori, dir);
				float uvd = lineClosestPoints(*ptr, *ptr + axis, ori, ori + dir).x
					- lineClosestPoints(*ptr, *ptr + axis, mouseRayOrigin, mouseRayOrigin + mouseRayDir).x;
				*ptr += (2.f * uvd) * axis;
			},
			[=](float dist) {
				queueGizmoLine(*ptr, *ptr + (0.2f * scale) * axis, thickness * 1.3f, color);
			}, nullptr
			});
	}

	inline void translateGizmoPlane(vec3* ptr, mat3 basis, int axis, vec4 color, float scale) {
		float size = 0.04f * scale;

		mat4 model = mat4(size * basis);
		std::swap(model[axis], model[2]);
		model[3] = vec4(*ptr + 0.5f * (vec3)model[0] + 0.5f * (vec3)model[1], 1);
		queueGizmoSquare(model, color);

		addGizmoSquareCallback(model, ptr, {
			[=](vec2 drag, float dist) {
				mat4 inv = inverse(model);

				vec3 ori, dir;
				ndcToWorldRay(mouseNDC + drag, ori, dir);

				ori = inv * vec4(ori, 1);
				dir = inv * vec4(dir, 0);
				vec2 a = (vec2)ori - (vec2)vec2(dir) * (ori.z / dir.z);

				ori = inv * vec4(mouseRayOrigin, 1);
				dir = inv * vec4(mouseRayDir, 0);
				vec2 b = (vec2)ori - (vec2)vec2(dir) * (ori.z / dir.z);

				vec2 diff = 2.f * (a - b);

				*ptr += (vec3)(model * vec4(diff, 0, 0));
			},
			[=](float dist) {
				mat4 m = model;
				m[3] = vec4(*ptr, 1);
				queueGizmoSquare(m, color);
			}, nullptr
			});
	}

	void queueGizmoTranslate(const char* label, vec3* ptr, bool interactive) {
		float scale = 0.5f * length(viewMat[3]);

		if (ptr == pendingPtr) {
			translateGizmoAxis(ptr, vec3(1, 0, 0), vec4(0.8, 0.1, 0.1, 1), scale);
			translateGizmoAxis(ptr, vec3(0, 1, 0), vec4(0.1, 0.8, 0.1, 1), scale);
			translateGizmoAxis(ptr, vec3(0, 0, 1), vec4(0.1, 0.1, 0.8, 1), scale);

			mat3 basis = mat3(1);
			translateGizmoPlane(ptr, basis, 0, vec4(0.8, 0.1, 0.1, 0.2), scale);
			translateGizmoPlane(ptr, basis, 1, vec4(0.1, 0.8, 0.1, 0.2), scale);
			translateGizmoPlane(ptr, basis, 2, vec4(0.1, 0.1, 0.8, 0.2), scale);
		}

		if (label) {
			mat4 oldModel = modelMat;
			modelMat = transpose((mat3)viewMat);
			modelMat[3] = vec4(*ptr, 1);
			defFont->render(label, scale * 3.f, vec4(1, 1, 1, 1), Bound<2>(vec2(scale * 0.01f, 0))
				, TextWrap::NONE, TextWrap::NONE, TextJustification::LEFT, TextJustification::CENTER);
			modelMat = oldModel;
		}

		queueGizmoSquare(*ptr, scale * 0.004f, vec4(0.8, 0.8, 0.8, 1));
		addGizmoCircleCallback(*ptr, scale * 0.01f, ptr, {
			[=](vec2 drag, float dist) {
				vec3 ori, dir;
				ndcToWorldRay(mouseNDC + drag, ori, dir);
				*ptr += 2.f * ((ori + dist * dir) - (mouseRayOrigin + dist * mouseRayDir));
			},
			[=](float dist) {
				queueGizmoSquare(*ptr, scale * 0.006f, vec4(0.8, 0.8, 0.8, 1));
			}, nullptr
			});
	}

	void addGizmoLineCallback(vec3 a, vec3 b, float thickness, void* ptr, GizmoCallback callback) {
		// Test intersection and register callback
		vec2 uv = lineClosestPoints(a, b, mouseRayOrigin, mouseRayOrigin + mouseRayDir);
		uv.x = glm::clamp(uv.x, 0.f, 1.f);
		float dist = length(mix(a, b, uv.x) - mix(mouseRayOrigin, mouseRayOrigin + mouseRayDir, uv.y));
		processCallback(uv.y, ptr, callback);
		if (dist < thickness)
			recordCallback(uv.y, ptr, callback);
	}

	void addGizmoCircleCallback(vec3 a, float radius, void* ptr, GizmoCallback callback) {
		float dist = dot(a - mouseRayOrigin, mouseRayDir);
		processCallback(dist, ptr, callback);
		if (length2(mouseRayOrigin + dist * mouseRayDir - a) < radius * radius)
			recordCallback(dist - radius, ptr, callback);
	}

	void addGizmoSquareCallback(mat4 modelMat, void* ptr, GizmoCallback callback) {
		mat4 inv = inverse(modelMat);
		vec3 ori = inv * vec4(mouseRayOrigin, 1);
		vec3 dir = inv * vec4(mouseRayDir, 0);

		float dist = -ori.z / dir.z;
		vec2 uv = (vec2)ori + (vec2)vec2(dir) * dist;
		processCallback(dist, ptr, callback);
		if (0 <= uv.x && uv.x <= 1 && 0 <= uv.y && uv.y <= 1)
			recordCallback(dist, ptr, callback);
	}

	void queueGizmoLine(vec3 a, vec3 b, float thickness, vec4 color) {
		vec3 cam = cameraPos();
		vec3 dir = normalize(mix(a, b, lineClosestPoints(a, b, cam)) - cameraPos());

		GizmoData g;
		g.model[1] = vec4(b - a, 0);
		g.model[0] = vec4(normalize(cross((vec3)g.model[1], dir)), 0);
		g.model[2] = vec4(dir, 0);
		g.model[0] *= thickness;
		g.model[3] = vec4(a - 0.5f * (vec3)g.model[0], 1);
		g.color = color;
		g.model = projMat * viewMat * g.model;
		gizmoRenderQueue.push_back(g);
	}

	void queueGizmoSquare(mat4 model, vec4 color) {
		vec3 dir = normalize((vec3)model[3] - cameraPos());
		mat3 basis = orthoBasisZ(dir);
		GizmoData g;
		g.model = model;
		g.color = color;
		g.model = projMat * viewMat * g.model;
		gizmoRenderQueue.push_back(g);
	}

	void queueGizmoSquare(vec3 a, float size, vec4 color) {
		vec3 dir = normalize(a - cameraPos());
		mat3 basis = orthoBasisZ(dir);
		GizmoData g;
		g.model = (mat4)(basis * mat3(size));
		g.model[3] = vec4(a - 0.5f * (vec3)g.model[0] - 0.5f * (vec3)g.model[1], 1);
		g.color = color;
		g.model = projMat * viewMat * g.model;
		gizmoRenderQueue.push_back(g);
	}

	void queueGizmoLineScreenspace(vec2 a, vec2 b, float thickness, vec4 color) {
		GizmoData g;
		g.model[0] = vec4(b - a, 0, 0);
		g.model[1] = vec4(-g.model[0].y * thickness, g.model[0].x * thickness, 0, 0);
		g.model[2] = vec4(0);
		float aspect = getAspect();
		a.x *= aspect;
		g.model[3] = vec4(vec3(a, 0) - 0.5f * (vec3)g.model[1], 1);
		g.model = glm::ortho(-aspect, aspect, -1.f, 1.f) * g.model;
		g.color = color;
		gizmoRenderQueue.push_back(g);
	}

	void queueGizmoSquareScreenspace(vec2 a, float size, vec4 color) {
		GizmoData g;
		g.model = (mat4)mat3(size);
		float aspect = getAspect();
		a.x *= aspect;
		g.model[3] = vec4(vec3(a, 0) - 0.5f * (vec3)g.model[0] - 0.5f * (vec3)g.model[1], 1);
		g.model = glm::ortho(-aspect, aspect, -1.f, 1.f) * g.model;
		g.color = color;
		gizmoRenderQueue.push_back(g);
	}

	void initGizmos() {
		gizmoBuffer = new ComputeBuffer(sizeof(GizmoData), MAX_GIZMO_PER_DRAW, GL_DYNAMIC_DRAW);
	}

	void recordCallback(float dist, void* ptr, GizmoCallback callback) {
		if (closestDist > dist && dist >= 0) {
			closestDist = dist;
			closestPtr = ptr;
			closestCallback = callback;
		}
	}

	void processCallback(float dist, void* ptr, GizmoCallback callback) {
		// Call callbacks and reset if needed
		if (ptr == pendingPtr) {
			if (onDragPending && pendingCallback.onDrag) {
				pendingCallback.onDrag(ndcDragAmount, dist);
				ndcDragAmount = vec2(0);
				onDragPending = false;
			}
			if (onMouseOverPending && pendingCallback.onMouseOver) {
				pendingCallback.onMouseOver(dist);
				onMouseOverPending = false;
			}
			if (onClickPending && pendingCallback.onClick) {
				pendingCallback.onClick(dist);
				onClickPending = false;
			}
		}
	}

	void gizmoMouseButtonCallback(int button, int action, int mode) {
		if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_PRESS)
			isClicking = true;
		else if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_RELEASE)
			isLifting = true;
	}

	void gizmoCursorPosCallback(double xp, double yp) {
		vec2 ndc = vec2(xp / windowRes.x, yp / windowRes.y) * 2.f - 1.f;
		ndc.y = -ndc.y;

		ndcToWorldRay(ndc, mouseRayOrigin, mouseRayDir);

		ndcDragAmount = ndc - mouseNDC;
		mouseNDC = ndc;
	}

	void queueGizmo(GizmoData g) {
		gizmoRenderQueue.push_back(g);
	}

	void renderGizmos() {
		// Process events
		onDragPending = onMouseOverPending = onClickPending = false;
		if (isDragging) {	// Dragging
			onDragPending = ndcDragAmount.x != 0 || ndcDragAmount.y != 0;
			onMouseOverPending = true;
		}
		else if (closestPtr) {				// Not draging but about to start draging
			isDragging = isClicking;
			pendingPtr = closestPtr;
			pendingCallback = closestCallback;
			onMouseOverPending = true;
		}
		else if (isClicking) {
			pendingPtr = nullptr;
			pendingCallback = {};
		}

		if (isLifting) {					// Process click and end drag
			isDragging = false;
			onClickPending = true;
		}

		// Reset
		closestDist = INFINITY;
		closestPtr = nullptr;
		closestCallback = {};
		isLifting = isClicking = false;

		// Render all gizmos
		if (gizmoRenderQueue.size()) {
			glTempVar<GL_BLEND_ALPHA> blend(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glTempVar<GL_CULL_FACE> cull(false);
			glTempVar<GL_DEPTH_WRITEMASK> zwrite(false);

			static auto shader = get<Shader>("KittenEngine\\shaders\\gizmo.glsl");
			shader->use();
			gizmoBuffer->bind(3);
			glBindVertexArray(defMesh->VAO);

			int i = 0;
			for (; i < (int)gizmoRenderQueue.size() - MAX_GIZMO_PER_DRAW; i += MAX_GIZMO_PER_DRAW) {
				gizmoBuffer->upload(&gizmoRenderQueue[i]);
				glDrawElementsInstanced(shader->drawMode(), (GLsizei)defMesh->indices.size(), GL_UNSIGNED_INT, 0, MAX_GIZMO_PER_DRAW);
			}
			int leftOver = (int)gizmoRenderQueue.size() - i;
			gizmoBuffer->upload(&gizmoRenderQueue[i], leftOver);
			glDrawElementsInstanced(shader->drawMode(), (GLsizei)defMesh->indices.size(), GL_UNSIGNED_INT, 0, leftOver);

			shader->unuse();
			gizmoRenderQueue.clear();
		}
	}
}