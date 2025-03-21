#include <glad/glad.h>
#include "../includes/modules/FrameBuffer.h"

namespace Kitten {
	FrameBuffer::FrameBuffer() {}
	FrameBuffer::FrameBuffer(int width, int height, int numGBuffers)
		:managed(true), width(width), height(height) {
		glGenFramebuffers(1, &glHandle);
		glBindFramebuffer(GL_FRAMEBUFFER, glHandle);

		glGenTextures(1, &depthStencil);
		glBindTexture(GL_TEXTURE_2D, depthStencil);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glBindTexture(GL_TEXTURE_2D, 0);

		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthStencil, 0);

		for (size_t i = 0; i < numGBuffers; i++) {
			buffs[i] = new Texture(width, height);
			glFramebufferTexture2D(GL_FRAMEBUFFER, GLenum(GL_COLOR_ATTACHMENT0 + i), GL_TEXTURE_2D, buffs[i]->glHandle, 0);
		}
		if (numGBuffers < 8)
			buffs[numGBuffers] = nullptr;

		if (numGBuffers == 0)
			glDrawBuffer(GL_NONE);

		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
			printf("error: fbo creation failure!!\n");

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}

	FrameBuffer::FrameBuffer(Texture** textures, int numGBuffers) {
		width = textures[0]->width;
		height = textures[0]->height;
		managed = true;

		glGenFramebuffers(1, &glHandle);
		glBindFramebuffer(GL_FRAMEBUFFER, glHandle);

		glGenTextures(1, &depthStencil);
		glBindTexture(GL_TEXTURE_2D, depthStencil);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glBindTexture(GL_TEXTURE_2D, 0);

		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthStencil, 0);

		for (size_t i = 0; i < numGBuffers; i++) {
			buffs[i] = textures[i];
			glFramebufferTexture2D(GL_FRAMEBUFFER, GLenum(GL_COLOR_ATTACHMENT0 + i), GL_TEXTURE_2D, buffs[i]->glHandle, 0);
		}
		if (numGBuffers < 8)
			buffs[numGBuffers] = nullptr;

		if (numGBuffers == 0)
			glDrawBuffer(GL_NONE);

		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
			printf("error: fbo creation failure!!\n");

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}

	FrameBuffer::~FrameBuffer() {
		if (managed) {
			for (int i = 0; buffs[i] && i < 8; i++)
				delete buffs[i];
			glDeleteTextures(1, &depthStencil);
			glDeleteFramebuffers(1, &glHandle);
		}
	}

	void FrameBuffer::clear(vec4 color, bool clearDepth) {
		glClearColor(color.x, color.y, color.z, color.w);
		if (clearDepth) glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		else glClear(GL_COLOR_BUFFER_BIT);
	}

	void FrameBuffer::clearDepth() {
		glClear(GL_DEPTH_BUFFER_BIT);
	}

	void FrameBuffer::resize(glm::ivec2 res) {
		resize(res.x, res.y);
	}

	void FrameBuffer::resize(int nw, int nh) {
		if (nw == width && nh == height) return;
		width = nw;
		height = nh;
		bind();

		int i = 0;
		while (buffs[i] && i < 8) buffs[i++]->resize(nw, nh);

		glBindTexture(GL_TEXTURE_2D, depthStencil);
		glTexImage2D(
			GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8, width, height, 0,
			GL_DEPTH_STENCIL, GL_UNSIGNED_INT_24_8, NULL
		);
		glBindTexture(GL_TEXTURE_2D, 0);
		unbind();
	}

	void FrameBuffer::bind() {
		glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &lastHandle);
		glGetFloatv(GL_VIEWPORT, (GLfloat*)&lastViewport);
		glViewport(0, 0, width, height);
		glBindFramebuffer(GL_FRAMEBUFFER, glHandle);
	}

	void FrameBuffer::unbind() {
		glBindFramebuffer(GL_FRAMEBUFFER, lastHandle);
		glViewport((GLsizei)lastViewport.x, (GLsizei)lastViewport.y, (GLsizei)lastViewport.z, (GLsizei)lastViewport.w);
	}
}