#pragma once
// Jerry Hsu, 2021

#include <string>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glad/glad.h>
#include "Shader.h"

namespace Kitten {
	class Texture {
	public:
		bool isCubeMap = false;
		unsigned int glHandle;
		int width;
		int height;
		GLenum deviceFormat;	// Format in GPU memory
		GLenum hostFormat;		// Internal data format
		GLenum hostDataType;

		ivec4 borders = ivec4(0);
		float ratio;
		unsigned char* rawData = nullptr;

		Texture();
		Texture(int width, int height, GLenum deviceFormat = GL_RGBA, GLenum hostFormat = GL_RGBA8, GLenum hostDataType = GL_UNSIGNED_BYTE);
		Texture(Texture* xpos, Texture* xneg, Texture* ypos, Texture* yneg, Texture* zpos, Texture* zneg);
		~Texture();

		void bind(int index = 0);
		void debugBlit(Kitten::Shader* shader = nullptr);
		vec4 sample(vec2 uv);
		void genMipmap();
		void resize(int w, int h);
		void setAniso(float v);
		void setFilter(GLenum mode);
		void setWrap(GLenum mode);
		void save(const char* path);
	};

	void savePNG(const char* path, unsigned char* data, int width, int height);
}