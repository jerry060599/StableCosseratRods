#pragma once
// Jerry Hsu, 2021

#include <string>
#include <vector>
#include <glad/glad.h> 
#include <glm/glm.hpp>
#include "ComputeBuffer.h"
#include "Rotor.h"

namespace Kitten {
	using namespace std;

	enum class ShaderType {
		VERT = 1,
		FRAG = 2,
		GEOM = 4,
		TESS = 8,
		COMP = 16
	};

	class Shader {
	public:
		~Shader();

		unsigned int glHandle = 0;
		int type = 0;
		bool tryLinked = false;
		vector<unsigned int> unlinkedHandles;

		bool link();
		void use();
		void unuse();
		GLenum drawMode();

		void dispatchCompute(int numThreads);
		void dispatchCompute(glm::ivec2 numThreads);
		void dispatchCompute(glm::ivec3 numThreads);

		void setBuffer(const char* name, ComputeBuffer* buffer);
		void setBool(const char* name, bool v);
		void setInt(const char* name, int v);
		void setFloat(const char* name, float v);
		void setFloat2(const char* name, glm::vec2 v);
		void setFloat3(const char* name, glm::vec3 v);
		void setFloat4(const char* name, glm::vec4 v);
		void setRotor(const char* name, Rotor v);
		void setMat2(const char* name, glm::mat2 v);
		void setMat3(const char* name, glm::mat3 v);
		void setMat4(const char* name, glm::mat4 v);
	};

	bool compileShader(string path, GLenum type, unsigned int* handle);
}