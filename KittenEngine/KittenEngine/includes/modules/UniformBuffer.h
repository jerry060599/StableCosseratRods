#pragma once
// Jerry Hsu, 2021

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glad/glad.h>

namespace Kitten {
	template<typename T>
	class UniformBuffer {
	public:
		unsigned int glHandle;

		UniformBuffer() {
			glGenBuffers(1, &glHandle);
			glBindBuffer(GL_UNIFORM_BUFFER, glHandle);
			glBufferData(GL_UNIFORM_BUFFER, sizeof(T), NULL, GL_STATIC_DRAW);
			glBindBuffer(GL_UNIFORM_BUFFER, 0);
		}
		~UniformBuffer() {
			glDeleteBuffers(1, &glHandle);
		}

		void bind(int loc) {
			glBindBufferRange(GL_UNIFORM_BUFFER, loc, glHandle, 0, sizeof(T));
		}

		void upload(T& data) {
			glBindBuffer(GL_UNIFORM_BUFFER, glHandle);
			glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(T), &data);
			glBindBuffer(GL_UNIFORM_BUFFER, 0);
		}
	};
}