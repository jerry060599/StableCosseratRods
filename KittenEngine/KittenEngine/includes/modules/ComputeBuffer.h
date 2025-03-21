#pragma once
// Jerry Hsu, 2021

#include <string>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glad/glad.h>

#if __has_include("cuda_runtime.h")
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#endif

namespace Kitten {
	class ComputeBuffer {
	public:
		unsigned int glHandle;
		size_t elementSize, size;
		GLenum usage;

		ComputeBuffer() = delete;
		ComputeBuffer(size_t elementSize, size_t count, GLenum usage = GL_DYNAMIC_READ);
		~ComputeBuffer();

		void bind(int loc);
		void resize(size_t newSize);
		void upload(void* src);
		void upload(void* src, size_t count);
		void download(void* dst);
		void download(void* dst, size_t count);
	};

#ifdef __CUDA_RUNTIME_H__
	class CudaComputeBuffer : public ComputeBuffer {
	public:
		cudaGraphicsResource* cudaRes;
		void* cudaPtr;

		CudaComputeBuffer() = delete;
		CudaComputeBuffer(size_t elementSize, size_t count, GLenum usage = GL_DYNAMIC_READ) :
			ComputeBuffer(elementSize, count, usage) {
			cudaGraphicsGLRegisterBuffer(&cudaRes, glHandle, cudaGraphicsRegisterFlagsNone);
			cudaGraphicsMapResources(1, &cudaRes);

			size_t size;
			cudaGraphicsResourceGetMappedPointer(&cudaPtr, &size, cudaRes);
		}

		~CudaComputeBuffer() {
			cudaGraphicsUnmapResources(1, &cudaRes, 0);
			cudaGraphicsUnregisterResource(cudaRes);
		}
	};
#endif

}