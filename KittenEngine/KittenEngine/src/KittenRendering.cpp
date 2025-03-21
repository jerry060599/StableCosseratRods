#include "../includes/modules/KittenInit.h"

#include "../includes/modules/KittenRendering.h"
#include "../includes/modules/KittenPreprocessor.h"
#include "../includes/modules/UniformBuffer.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

using namespace glm;

namespace Kitten {
	mat4 projMat;
	mat4 invProjMat;
	mat4 viewMat;
	mat4 invViewMat;
	mat4 modelMat(1);

	Material defMaterial{ vec4(1), vec4(1), vec4(1), vec4(1) };

	vector<UBOLight> lights;
	UBOLight ambientLight = { vec4(0.f, 0.f, 0.f, 1.f), vec3(0), 0, 0, 0, 0, 0, vec3(0), (int)KittenLight::AMBIENT };

	UniformBuffer<UBOCommon>* uboCommon;
	UniformBuffer<UBOModel>* uboModel;
	UniformBuffer<UBOMat>* uboMat;
	UniformBuffer<UBOLight>* uboLight;

	int shadowRes = 2048;
	float shadowDist = 50.f;

	Texture* defTexture;
	Texture* defCubemap;
	Mesh* defMesh, * defMeshPoly;
	Font* defFont;

	Shader* defBaseShader;
	Shader* defForwardShader;
	Shader* defUnlitShader;
	Shader* defEnvShader;
	Shader* defBlitShader;

	bool initialized = false;
	vector<FrameBuffer*> shadowMaps;

	ivec2 windowRes;

	void allocShadowMaps() {
		while (lights.size() > shadowMaps.size())
			shadowMaps.push_back(new FrameBuffer(shadowRes, shadowRes, 0));
		if (lights.size() < shadowMaps.size()) {
			for (size_t i = lights.size(); i < shadowMaps.size(); i++)
				delete shadowMaps[i];
			shadowMaps.resize(lights.size());
		}
		for (size_t i = 0; i < shadowMaps.size(); i++) {
			shadowMaps[i]->resize(shadowRes, shadowRes);
			mat4 proj(1);
			if (lights[i].type == (int)KittenLight::DIR)
				proj = glm::ortho(-shadowDist, shadowDist, -shadowDist, shadowDist, -shadowDist, shadowDist)
				* rotateView(lights[i].dir) * glm::translate(mat4(1), -lights[i].pos);
			else if (lights[i].type == (int)KittenLight::SPOT)
				proj = glm::perspective(2 * (180 / glm::pi<float>()) * glm::acos(lights[i].spread), 1.f, lights[i].radius, shadowDist)
				* rotateView(lights[i].dir) * glm::translate(mat4(1), -lights[i].pos);
			lights[i].shadowProj = proj;
		}
	}

	void clearShadowMaps() {
		for (auto fb : shadowMaps) {
			fb->bind();
			glClear(GL_DEPTH_BUFFER_BIT);
			fb->unbind();
		}
	}

	double getTime() {
		return glfwGetTime();
	}

	float getAspect() {
		return windowRes.x / (float)windowRes.y;
	}

	bool shouldClose() {
		return glfwWindowShouldClose(window);
	}

	void checkErr(const char* tag) {
		unsigned int error = glGetError();
		if (error != GL_NO_ERROR)
			if (tag)
				printf("GL error %d at %s\n", error, tag);
			else
				printf("GL error %d\n", error);
	}

	void initRender() {
		if (initialized) return;
		initialized = true;

		includePaths.push_back("KittenEngine\\shaders");

		loadDirectory("KittenEngine\\shaders");
		defBaseShader = get<Kitten::Shader>("KittenEngine\\shaders\\blingBase.glsl");
		defForwardShader = get<Kitten::Shader>("KittenEngine\\shaders\\blingForward.glsl");
		defUnlitShader = get<Kitten::Shader>("KittenEngine\\shaders\\unlit.glsl");
		defEnvShader = get<Kitten::Shader>("KittenEngine\\shaders\\env.glsl");
		defBlitShader = get<Kitten::Shader>("KittenEngine\\shaders\\blit.glsl");

		loadDirectory("KittenEngine\\fonts");
		defFont = get<Kitten::Font>("KittenEngine\\fonts\\Quicksand_Regular.ttf");

		uboCommon = new UniformBuffer<UBOCommon>;
		uboModel = new UniformBuffer<UBOModel>;
		uboMat = new UniformBuffer<UBOMat>;
		uboLight = new UniformBuffer<UBOLight>;

		glEnable(GL_CULL_FACE);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);

		glDepthFunc(GL_LEQUAL);
		glDepthMask(GL_TRUE);
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

		const float patchLvl[]{ 16, 16, 16, 16 };
		glPatchParameteri(GL_PATCH_VERTICES, 3);
		glPatchParameterfv(GL_PATCH_DEFAULT_OUTER_LEVEL, patchLvl);
		glPatchParameterfv(GL_PATCH_DEFAULT_INNER_LEVEL, patchLvl);

		{// Gen white texture for binding null textures
			unsigned int handle;
			glGenTextures(1, &handle);
			glBindTexture(GL_TEXTURE_2D, handle);

			unsigned int data = ~0;
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, &data);
			glGenerateMipmap(GL_TEXTURE_2D);

			glBindTexture(GL_TEXTURE_2D, 0);

			Texture* img = new Texture;
			img->glHandle = handle;
			img->width = 1;
			img->height = 1;
			img->deviceFormat = GL_RGBA;
			img->hostFormat = GL_RGBA8;
			img->hostDataType = GL_UNSIGNED_BYTE;
			img->ratio = 1;
			img->borders = ivec4(0);
			img->rawData = new unsigned char[4];
			img->rawData[0] = ~0;
			img->rawData[1] = ~0;
			img->rawData[2] = ~0;
			img->rawData[3] = ~0;
			resources["white.tex"] = img;
			defTexture = img;
		}

		// Gen standard quad for FBO bliting
		defMesh = genQuadMesh(1, 1);
		defMeshPoly = genQuadMesh(1, 1);
		defMeshPoly->polygonize();

		defCubemap = new Texture(defTexture, defTexture, defTexture, defTexture, defTexture, defTexture);
	}

	inline void uploadUniformBuff(unsigned int handle, void* data, size_t size) {
		glBindBuffer(GL_UNIFORM_BUFFER, handle);
		glBufferSubData(GL_UNIFORM_BUFFER, 0, size, data);
		glBindBuffer(GL_UNIFORM_BUFFER, 0);
	}

	void uploadUBOCommonBuff() {
		UBOCommon dat;

		dat.projMat = projMat;
		invProjMat = dat.projMatInv = inverse(projMat);
		dat.viewMat = viewMat;
		invViewMat = dat.viewMatInv = inverse(viewMat);

		dat.viewMat_n = normalTransform(dat.viewMat);
		dat.vpMat = dat.projMat * dat.viewMat;
		dat.vpMatInv = inverse(dat.vpMat);

		uboCommon->upload(dat);
	}

	void startRender() {
		initRender();

		uboCommon->bind(0);
		uboModel->bind(1);
		uboMat->bind(2);
		uboLight->bind(3);

		ambientLight.type = (int)KittenLight::AMBIENT;
		uploadUBOCommonBuff();

		allocShadowMaps();
		clearShadowMaps();

		modelMat = mat4(1);
		checkErr("render_start");
	}

	void startFrame() {
		glfwPollEvents();
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		checkErr("frame_start");
	}

	void endFrame() {
		renderGizmos();
		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		checkErr("frame_end");
		glfwSwapBuffers(window);
	}

	void startRenderMesh(mat4 transform) {
		UBOModel dat;
		dat.modelMat = modelMat * transform;
		dat.modelMatInv = inverse(dat.modelMat);
		dat.modelMat_n = normalTransform(dat.modelMat);
		uboModel->upload(dat);
	}

	void startRenderMaterial(Material* mat) {
		if (!mat) mat = &defMaterial;

		uboMat->upload(mat->props);

		for (size_t i = 0; i <= MAT_TEX5; i++) {
			glActiveTexture((GLenum)(GL_TEXTURE0 + i));
			glBindTexture(GL_TEXTURE_2D, mat->texs[i] ? mat->texs[i]->glHandle : defTexture->glHandle);
		}
		glActiveTexture((GLenum)(GL_TEXTURE0 + MAT_CUBEMAP));
		glBindTexture(GL_TEXTURE_CUBE_MAP, mat->texs[MAT_CUBEMAP] ? mat->texs[MAT_CUBEMAP]->glHandle : defCubemap->glHandle);
	}

	void render(Mesh* mesh, Shader* base) {
		renderForward(mesh, base);
	}

	void uploadAmbientLight() {
		uboLight->upload(ambientLight);
	}

	void renderAdditive(Mesh* mesh, Shader* base) {
		startRenderMesh(mesh->defTransform);
		startRenderMaterial(mesh->defMaterial);
		uploadAmbientLight();

		glTempVar<GL_BLEND_ALPHA> blend(GL_ONE, GL_ONE);
		if (!base) base = defUnlitShader;
		base->use();
		glBindVertexArray(mesh->VAO);
		glDrawElements(base->drawMode(), (GLsizei)mesh->indices.size(), GL_UNSIGNED_INT, 0);
	}

	void renderLine(Mesh* mesh, Shader* base) {
		startRenderMesh(mesh->defTransform);
		startRenderMaterial(mesh->defMaterial);
		uploadAmbientLight();

		glTempVar<GL_BLEND_ALPHA> blend(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		if (!base) base = defUnlitShader;
		base->use();
		glBindVertexArray(mesh->VAO);
		glDrawElements(GL_LINES, (GLsizei)mesh->indices.size(), GL_UNSIGNED_INT, 0);
	}

	void renderInstanced(Mesh* mesh, int count, Shader* base) {
		startRenderMesh(mesh->defTransform);
		startRenderMaterial(mesh->defMaterial);
		uploadAmbientLight();

		glTempVar<GL_BLEND_ALPHA> blend(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		if (!base) base = defUnlitShader;
		base->use();
		glBindVertexArray(mesh->VAO);
		glDrawElementsInstanced(base->drawMode(), (GLsizei)mesh->indices.size(), GL_UNSIGNED_INT, 0, count);
	}

	inline bool canSkipBase() {
		return lights.size() == 1 && trace(vec3(ambientLight.col)) * ambientLight.col.w <= 0;
	}

	void renderForward(Mesh* mesh, Shader* base, Shader* light) {
		startRenderMesh(mesh->defTransform);
		startRenderMaterial(mesh->defMaterial);
		uploadAmbientLight();

		bool skipBase = canSkipBase() && light;

		glTempVar<GL_BLEND_ALPHA> blend(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glBindVertexArray(mesh->VAO);
		if (!skipBase) {
			if (!base) base = defUnlitShader;
			base->use();
			glDrawElements(base->drawMode(), (GLsizei)mesh->indices.size(), GL_UNSIGNED_INT, 0);
		}

		if (light) {
			glTempVar<GL_BLEND_ALPHA> blend(GL_ONE, skipBase ? GL_ZERO : GL_ONE);
			glTempVar<GL_DEPTH_WRITEMASK> zwrite(skipBase);

			light->use();
			for (size_t i = 0; i < lights.size(); i++) {
				if (lights[i].type == (int)KittenLight::SPOT
					|| lights[i].type == (int)KittenLight::DIR) {
					glActiveTexture((GLenum)(GL_TEXTURE0 + MAT_SHADOW));
					glBindTexture(GL_TEXTURE_2D, shadowMaps[i]->depthStencil);
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);
				}

				uboLight->upload(lights[i]);
				glDrawElements(light->drawMode(), (GLsizei)mesh->indices.size(), GL_UNSIGNED_INT, 0);
			}
		}

		glBindVertexArray(0);
		glUseProgram(0);
	}

	void renderInstancedForward(Mesh* mesh, int count, Shader* base, Shader* light) {
		startRenderMesh(mesh->defTransform);
		startRenderMaterial(mesh->defMaterial);
		uploadAmbientLight();

		bool skipBase = canSkipBase() && light;

		glTempVar<GL_BLEND_ALPHA> blend(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glBindVertexArray(mesh->VAO);

		if (!skipBase) {
			if (!base) base = defUnlitShader;
			base->use();
			glDrawElementsInstanced(base->drawMode(), (GLsizei)mesh->indices.size(), GL_UNSIGNED_INT, 0, count);
		}

		if (light) {
			glTempVar<GL_BLEND_ALPHA> blend(GL_ONE, skipBase ? GL_ZERO : GL_ONE);
			glTempVar<GL_DEPTH_WRITEMASK> zwrite(skipBase);

			light->use();
			for (size_t i = 0; i < lights.size(); i++) {
				if (lights[i].type == (int)KittenLight::SPOT
					|| lights[i].type == (int)KittenLight::DIR) {
					glActiveTexture((GLenum)(GL_TEXTURE0 + MAT_SHADOW));
					glBindTexture(GL_TEXTURE_2D, shadowMaps[i]->depthStencil);
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);
				}

				uboLight->upload(lights[i]);
				glDrawElementsInstanced(light->drawMode(), (GLsizei)mesh->indices.size(), GL_UNSIGNED_INT, 0, count);
			}
		}

		glBindVertexArray(0);
		glUseProgram(0);
	}

	void renderShadows(Mesh* mesh, Shader* base) {
		startRenderMesh(mesh->defTransform);
		startRenderMaterial(mesh->defMaterial);

		mat4 oldView = viewMat;
		mat4 oldProj = projMat;
		viewMat = mat4(1);
		glTempVar<GL_CULL_FACE> cull(false);
		glTempVar<GL_COLOR_WRITEMASK> colwrite(bvec4(false));

		if (base == nullptr) base = defUnlitShader;
		base->use();
		glBindVertexArray(mesh->VAO);
		for (size_t i = 0; i < lights.size(); i++)
			if (lights[i].hasShadow
				&& (lights[i].type == (int)KittenLight::SPOT
					|| lights[i].type == (int)KittenLight::DIR)) {
				projMat = lights[i].shadowProj;
				uploadUBOCommonBuff();

				shadowMaps[i]->bind();
				glDrawElements(base->drawMode(), (GLsizei)mesh->indices.size(), GL_UNSIGNED_INT, 0);
				shadowMaps[i]->unbind();
			}

		viewMat = oldView;
		projMat = oldProj;
		uploadUBOCommonBuff();
		glBindVertexArray(0);
		glUseProgram(0);
	}

	void renderInstancedShadows(Mesh* mesh, int count, Shader* base) {
		startRenderMesh(mesh->defTransform);
		startRenderMaterial(mesh->defMaterial);

		mat4 oldView = viewMat;
		mat4 oldProj = projMat;
		viewMat = mat4(1);
		glTempVar<GL_CULL_FACE> cull(false);
		glTempVar<GL_COLOR_WRITEMASK> colwrite(bvec4(false));

		if (base == nullptr) base = defUnlitShader;
		base->use();
		glBindVertexArray(mesh->VAO);
		for (size_t i = 0; i < lights.size(); i++)
			if (lights[i].hasShadow
				&& (lights[i].type == (int)KittenLight::SPOT
					|| lights[i].type == (int)KittenLight::DIR)) {
				projMat = lights[i].shadowProj;
				uploadUBOCommonBuff();

				shadowMaps[i]->bind();
				glDrawElementsInstanced(base->drawMode(), (GLsizei)mesh->indices.size(), GL_UNSIGNED_INT, 0, count);
				shadowMaps[i]->unbind();
			}

		viewMat = oldView;
		projMat = oldProj;
		uploadUBOCommonBuff();
		glBindVertexArray(0);
		glUseProgram(0);
	}

	void renderEnv(Texture* cubemap) {
		glActiveTexture((GLenum)(GL_TEXTURE0 + MAT_CUBEMAP));
		glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap->glHandle);
		defEnvShader->use();
		glBindVertexArray(defMesh->VAO);
		glDrawElements(GL_TRIANGLES, (GLsizei)defMesh->indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
		glUseProgram(0);
	}

	void fixedUpdateAdapter(std::function<void(double)> dynamicUpdate, std::function<void(double)> fixedUpdate,
		double dt, double fixedDT, double& timeSinceFixed) {
		double timeLeft = dt;

		while (true) {
			float timeTillNextFixed = std::max(fixedDT - timeSinceFixed, 0.);
			if (timeLeft > timeTillNextFixed) {
				dynamicUpdate(timeTillNextFixed);
				timeLeft -= timeTillNextFixed;
				timeSinceFixed = 0;
				fixedUpdate(fixedDT);
			}
			else break;
		}

		dynamicUpdate(timeLeft);
		timeSinceFixed += timeLeft;
	}
}