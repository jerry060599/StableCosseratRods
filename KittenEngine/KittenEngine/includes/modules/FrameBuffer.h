#pragma once
// Jerry Hsu, 2021

#include <glm/glm.hpp>
#include "Texture.h"

namespace Kitten {
	class FrameBuffer {
	public:
		bool managed = false;
		int width, height;
		unsigned int glHandle = 0;
		unsigned int depthStencil = 0;
		Texture* buffs[8] = {};
		vec4 lastViewport;
		GLint lastHandle;

		FrameBuffer();
		FrameBuffer(int width, int height, int numGBuffers = 1);
		FrameBuffer(Texture** textures, int numGBuffers = 1);
		~FrameBuffer();

		void clear(vec4 color = vec4(0, 0, 0, 0), bool clearDepth = true);
		void clearDepth();

		void resize(glm::ivec2 res);
		void resize(int width, int height);
		void bind();
		void unbind();
	};
}