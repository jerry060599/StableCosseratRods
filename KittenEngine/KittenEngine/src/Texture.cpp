#include "../includes/modules/Texture.h"
#include "../includes/modules/KittenRendering.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

namespace Kitten {
	Texture::Texture() {}
	Texture::Texture(int width, int height, GLenum deviceFormat, GLenum hostFormat, GLenum hostDataType)
		:width(width), height(height), deviceFormat(deviceFormat), ratio(float(width) / height), hostFormat(hostFormat), hostDataType(hostDataType) {
		glGenTextures(1, &glHandle);
		glBindTexture(GL_TEXTURE_2D, glHandle);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glTexImage2D(GL_TEXTURE_2D, 0, hostFormat, width, height, 0, deviceFormat, hostDataType, nullptr);
		glBindTexture(GL_TEXTURE_2D, 0);
	}

	Texture::Texture(Texture* xpos, Texture* xneg, Texture* ypos, Texture* yneg, Texture* zpos, Texture* zneg) {
		isCubeMap = true;
		glGenTextures(1, &glHandle);
		glBindTexture(GL_TEXTURE_CUBE_MAP, glHandle);

		glTexImage2D(
			GL_TEXTURE_CUBE_MAP_POSITIVE_X,
			0, GL_RGB, xpos->width, xpos->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, xpos->rawData
		);
		glTexImage2D(
			GL_TEXTURE_CUBE_MAP_NEGATIVE_X,
			0, GL_RGB, xneg->width, xneg->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, xneg->rawData
		);
		glTexImage2D(
			GL_TEXTURE_CUBE_MAP_POSITIVE_Y,
			0, GL_RGB, ypos->width, ypos->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, ypos->rawData
		);
		glTexImage2D(
			GL_TEXTURE_CUBE_MAP_NEGATIVE_Y,
			0, GL_RGB, yneg->width, yneg->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, yneg->rawData
		);
		glTexImage2D(
			GL_TEXTURE_CUBE_MAP_POSITIVE_Z,
			0, GL_RGB, zpos->width, zpos->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, zpos->rawData
		);
		glTexImage2D(
			GL_TEXTURE_CUBE_MAP_NEGATIVE_Z,
			0, GL_RGB, zneg->width, zneg->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, zneg->rawData
		);

		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	}

	Texture::~Texture() {
		glDeleteTextures(1, &glHandle);
		if (rawData) delete[] rawData;
	}

	void Texture::bind(int index) {
		glActiveTexture((GLenum)(GL_TEXTURE0 + index));
		glBindTexture(isCubeMap ? GL_TEXTURE_CUBE_MAP : GL_TEXTURE_2D, glHandle);
	}

	void Texture::debugBlit(Kitten::Shader* shader) {
		if (!shader) shader = Kitten::defBlitShader;
		Kitten::Texture*& matTex = Kitten::defMaterial.texs[0];
		auto tex = matTex;
		matTex = this;
		Kitten::render(Kitten::defMesh, shader);
		matTex = tex;
	}

	vec4 Texture::sample(vec2 uv) {
		if (uv.x < 0) return vec4(0, 0, 0, 0);
		if (uv.x > 1) return vec4(0, 0, 0, 0);
		if (uv.y < 0) return vec4(0, 0, 0, 0);
		if (uv.y > 1) return vec4(0, 0, 0, 0);
		if (!rawData) return vec4(1, 0, 1, 0);
		//return vec4(mix(mix(vec3(1, 0, 0), vec3(1, 0, 2), uv.x), mix(vec3(0, 0, 1), vec3(0, 0, 0), uv.x), uv.y), 1.f);
		//float w = uv.x < 0.75f ? uv.x / 0.75f : 1 - (uv.x - 0.75f) / 0.25f;
		//w *= uv.y < 0.75f ? uv.y / 0.75f : 1 - (uv.y - 0.75f) / 0.25f;
		//return mix(vec4(1, 1, 1, 1), vec4(1, 0, 0, 1), w);
		size_t x = std::min((size_t)floor(uv.x * width), size_t(width - 1));
		size_t y = std::min((size_t)floor(uv.y * height), size_t(height - 1));
		size_t idx = 4 * (y * width + x);
		switch (hostDataType) {
		case GL_UNSIGNED_BYTE:
			return vec4((float)rawData[idx + 0], (float)rawData[idx + 1], (float)rawData[idx + 2], (float)rawData[idx + 3]) / 255.f;
		}
		return vec4(0);
	}

	void Texture::genMipmap() {
		glBindTexture(GL_TEXTURE_2D, glHandle);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glGenerateMipmap(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 0);
	}

	void Texture::setAniso(float v) {
		glBindTexture(GL_TEXTURE_2D, glHandle);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY, v);
		glBindTexture(GL_TEXTURE_2D, 0);
	}

	void Texture::setFilter(GLenum mode) {
		glBindTexture(GL_TEXTURE_2D, glHandle);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, mode);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mode);
		glBindTexture(GL_TEXTURE_2D, 0);
	}

	void Texture::setWrap(GLenum mode) {
		glBindTexture(GL_TEXTURE_2D, glHandle);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, mode);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, mode);
		glBindTexture(GL_TEXTURE_2D, 0);
	}

	void Texture::save(const char* path) {
		unsigned char* data = new unsigned char[(size_t)(width * height * 4)];
		glMemoryBarrier(GL_FRAMEBUFFER_BARRIER_BIT);

		bind(0);
		glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);

		stbi_flip_vertically_on_write(true);
		stbi_write_png(path, width, height, 4, data, 0);
		Kitten::checkErr("texture_save");

		delete[] data;
	}

	void Texture::resize(int nw, int nh) {
		if (nw == width && nh == height) return;
		width = nw;
		height = nh;
		glBindTexture(GL_TEXTURE_2D, glHandle);
		glTexImage2D(GL_TEXTURE_2D, 0, deviceFormat, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
		glBindTexture(GL_TEXTURE_2D, 0);
	}

	void savePNG(const char* path, unsigned char* data, int width, int height) {
		stbi_write_png(path, width, height, 4, data, 0);
	}
}