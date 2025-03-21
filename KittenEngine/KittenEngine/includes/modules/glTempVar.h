#pragma once
// Jerry Hsu, 2022

#include "Shader.h"
#include <glad/glad.h>
#include <glm/glm.hpp>

#define GEN_BOOL_SPEC(key)											\
template <>															\
class glTempVar<key> {												\
public:																\
	const GLboolean oldVal;											\
	glTempVar(bool val) : oldVal(0) {								\
		glGetBooleanv(key, (GLboolean*)&oldVal);					\
		if(val) glEnable(key); else glDisable(key);					\
	}																\
	~glTempVar() {													\
		if(oldVal) glEnable(key); else glDisable(key);				\
	}																\
};

#define GEN_INT_SPEC(key, setfunc)									\
template <>															\
class glTempVar<key> {												\
public:																\
	const int oldVal;												\
	glTempVar(int val) : oldVal(0) {								\
		glGetIntegerv(key, (int*)&oldVal);							\
		setfunc(val);												\
	}																\
	~glTempVar() {													\
		setfunc(oldVal);											\
	}																\
};

#define GEN_INT_SPEC_PARAM(key, setfunc, param)						\
template <>															\
class glTempVar<key> {												\
public:																\
	const int oldVal;												\
	glTempVar(int val) : oldVal(0) {								\
		glGetIntegerv(key, (int*)&oldVal);							\
		setfunc(param, val);										\
	}																\
	~glTempVar() {													\
		setfunc(param, oldVal);										\
	}																\
};

#define GEN_FLOAT_SPEC(key, setfunc)								\
template <>															\
class glTempVar<key> {												\
public:																\
	const float oldVal;												\
	glTempVar(float val) : oldVal(0) {								\
		glGetFloatv(key, (float*)&oldVal);							\
		setfunc(val);												\
	}																\
	~glTempVar() {													\
		setfunc(oldVal);											\
	}																\
};

#define GEN_VEC2_SPEC(key, setfunc)									\
template <>															\
class glTempVar<key> {												\
public:																\
	const vec2 oldVal;												\
	glTempVar(vec2 val) : oldVal(0) {								\
		glGetFloatv(key, (float*)&oldVal);							\
		setfunc(val.x, val.y);										\
	}																\
	glTempVar(float a, float b) :									\
		glTempVar(vec2(a, b)) {}									\
	~glTempVar() {													\
		setfunc(oldVal.x, oldVal.y);								\
	}																\
};

#define GEN_VEC4_SPEC(key, setfunc)									\
template <>															\
class glTempVar<key> {												\
public:																\
	const vec4 oldVal;												\
	glTempVar(vec4 val) : oldVal(0) {								\
		glGetFloatv(key, (float*)&oldVal);							\
		setfunc(val.x, val.y, val.z, val.w);						\
	}																\
	glTempVar(float a, float b, float c, float d) :					\
		glTempVar(vec4(a, b, c, d)) {}								\
	~glTempVar() {													\
		setfunc(oldVal.x, oldVal.y, oldVal.z, oldVal.w);			\
	}																\
};

#define GEN_IVEC4_SPEC(key, setfunc)								\
template <>															\
class glTempVar<key> {												\
public:																\
	const ivec4 oldVal;												\
	glTempVar(ivec4 val) : oldVal(0) {								\
		glGetIntegerv(key, (int*)&oldVal);							\
		setfunc(val.x, val.y, val.z, val.w);						\
	}																\
	glTempVar(int a, int b, int c, int d) :							\
		glTempVar(ivec4(a, b, c, d)) {}								\
	~glTempVar() {													\
		setfunc(oldVal.x, oldVal.y, oldVal.z, oldVal.w);			\
	}																\
};

#define GEN_BVEC4_SPEC(key, setfunc)								\
template <>															\
class glTempVar<key> {												\
public:																\
	const GLboolean oldVal[4]{};									\
	glTempVar(bvec4 val) {											\
		glGetBooleanv(key, (GLboolean*)oldVal);						\
		setfunc(val.x, val.y, val.z, val.w);						\
	}																\
	glTempVar(bool a, bool b, bool c, bool d) :						\
		glTempVar(bvec4(a, b, c, d)) {}								\
	~glTempVar() {													\
		setfunc(oldVal[0], oldVal[1], oldVal[2], oldVal[3]);		\
	}																\
};

#define GEN_BLEND_SPEC(key)											\
template <>															\
class glTempVar<key> {												\
public:																\
	const ivec2 oldVal;												\
	glTempVar(ivec2 val) : oldVal(0) {								\
		glGetIntegerv(GL_BLEND_SRC_ALPHA, (int*)&oldVal.x);			\
		glGetIntegerv(GL_BLEND_DST_ALPHA, (int*)&oldVal.y);			\
		glBlendFunc(val.x, val.y);									\
	}																\
	glTempVar(GLenum sfactor, GLenum  dfactor) :					\
		glTempVar(ivec2(sfactor, dfactor)) {}						\
	~glTempVar() {													\
		glBlendFunc(oldVal.x, oldVal.y);							\
	}																\
};

#define GEN_STENCIL_OP_SPEC(key)									\
template <>															\
class glTempVar<key> {												\
public:																\
	const ivec3 oldVal;												\
	glTempVar(ivec3 val) : oldVal(0) {								\
		glGetIntegerv(GL_STENCIL_FAIL, (int*)&oldVal.x);			\
		glGetIntegerv(GL_STENCIL_PASS_DEPTH_FAIL, (int*)&oldVal.y);	\
		glGetIntegerv(GL_STENCIL_PASS_DEPTH_PASS, (int*)&oldVal.y);	\
		glStencilOp(val.x, val.y, val.z);							\
	}																\
	glTempVar(GLenum sfail, GLenum dpfail, GLenum dppass) :			\
		glTempVar(ivec3(sfail, dpfail, dppass)) {}					\
	~glTempVar() {													\
		glStencilOp(oldVal.x, oldVal.y, oldVal.z);					\
	}																\
};

#define GEN_STENCIL_FUNC_SPEC(key)									\
template <>															\
class glTempVar<key> {												\
public:																\
	const GLenum oldVal0;											\
	const GLint oldVal1;											\
	const GLuint oldVal2;											\
	glTempVar(GLenum func, GLint ref, GLuint mask) :				\
			oldVal0(0), oldVal1(0), oldVal2(0) {					\
		glGetIntegerv(GL_STENCIL_FUNC, (int*)&oldVal0);				\
		glGetIntegerv(GL_STENCIL_REF, (int*)&oldVal1);				\
		glGetIntegerv(GL_STENCIL_VALUE_MASK, (int*)&oldVal2);		\
		glStencilFunc(func, ref, mask);								\
	}																\
	~glTempVar() {													\
		glStencilFunc(oldVal0, oldVal1, oldVal2);					\
	}																\
};

namespace Kitten {
	using namespace glm;

	template<GLenum key>
	class glTempVar {};

	GEN_BOOL_SPEC(GL_BLEND);
	GEN_BOOL_SPEC(GL_CULL_FACE);
	GEN_BOOL_SPEC(GL_DEPTH_TEST);
	GEN_BOOL_SPEC(GL_DITHER);
	GEN_BOOL_SPEC(GL_LINE_SMOOTH);
	GEN_BOOL_SPEC(GL_PROGRAM_POINT_SIZE);
	GEN_BOOL_SPEC(GL_POLYGON_SMOOTH);
	GEN_BOOL_SPEC(GL_SCISSOR_TEST);
	GEN_BOOL_SPEC(GL_STENCIL_TEST);

	// For some reason GL_DEPTH_WRITEMASK is special
	template <>
	class glTempVar<GL_DEPTH_WRITEMASK> {
	public:
		const GLboolean oldVal;
		glTempVar(bool val) : oldVal(0) {
			glGetBooleanv(GL_DEPTH_WRITEMASK, (GLboolean*)&oldVal);
			glDepthMask(val);
		}
		~glTempVar() { glDepthMask(oldVal); }
	};

	GEN_INT_SPEC(GL_ACTIVE_TEXTURE, glActiveTexture);
	GEN_INT_SPEC(GL_COLOR_LOGIC_OP, glLogicOp);
	GEN_INT_SPEC(GL_CULL_FACE_MODE, glCullFace);
	// GEN_INT_SPEC(GL_CURRENT_PROGRAM, glUseProgram);
	GEN_INT_SPEC(GL_DEPTH_FUNC, glDepthFunc);
	GEN_INT_SPEC(GL_DRAW_BUFFER, glDrawBuffer);
	GEN_INT_SPEC(GL_LOGIC_OP_MODE, glLogicOp);
	GEN_INT_SPEC(GL_PRIMITIVE_RESTART_INDEX, glPrimitiveRestartIndex);
	GEN_INT_SPEC(GL_PROGRAM_PIPELINE_BINDING, glBindProgramPipeline);
	GEN_INT_SPEC(GL_PROVOKING_VERTEX, glProvokingVertex);
	GEN_INT_SPEC(GL_READ_BUFFER, glReadBuffer);
	GEN_INT_SPEC(GL_STENCIL_CLEAR_VALUE, glClearStencil);
	GEN_INT_SPEC(GL_VERTEX_ARRAY_BINDING, glBindVertexArray);
	GEN_INT_SPEC(GL_STENCIL_WRITEMASK, glStencilMask);

	template <>
	class glTempVar<GL_CURRENT_PROGRAM> {
	public:
		const int oldVal;

		glTempVar(int val) : oldVal(0) {
			glGetIntegerv(GL_CURRENT_PROGRAM, (int*)&oldVal);
			glUseProgram(val);
		}

		glTempVar(Kitten::Shader* shader) : oldVal(0) {
			glGetIntegerv(GL_CURRENT_PROGRAM, (int*)&oldVal);
			shader->use();
		}

		~glTempVar() {
			glUseProgram(oldVal);
		}
	};

	GEN_INT_SPEC_PARAM(GL_ARRAY_BUFFER_BINDING, glBindBuffer, GL_ARRAY_BUFFER);
	GEN_INT_SPEC_PARAM(GL_DRAW_FRAMEBUFFER_BINDING, glBindFramebuffer, GL_DRAW_FRAMEBUFFER);
	GEN_INT_SPEC_PARAM(GL_READ_FRAMEBUFFER_BINDING, glBindFramebuffer, GL_READ_FRAMEBUFFER);
	GEN_INT_SPEC_PARAM(GL_ELEMENT_ARRAY_BUFFER_BINDING, glBindBuffer, GL_ELEMENT_ARRAY_BUFFER);
	GEN_INT_SPEC_PARAM(GL_FRAGMENT_SHADER_DERIVATIVE_HINT, glHint, GL_FRAGMENT_SHADER_DERIVATIVE_HINT);
	GEN_INT_SPEC_PARAM(GL_LINE_SMOOTH_HINT, glHint, GL_LINE_SMOOTH_HINT);
	GEN_INT_SPEC_PARAM(GL_PACK_ALIGNMENT, glPixelStorei, GL_PACK_ALIGNMENT);
	GEN_INT_SPEC_PARAM(GL_PACK_IMAGE_HEIGHT, glPixelStorei, GL_PACK_IMAGE_HEIGHT);
	GEN_INT_SPEC_PARAM(GL_PACK_LSB_FIRST, glPixelStorei, GL_PACK_LSB_FIRST);
	GEN_INT_SPEC_PARAM(GL_PACK_ROW_LENGTH, glPixelStorei, GL_PACK_ROW_LENGTH);
	GEN_INT_SPEC_PARAM(GL_PACK_SKIP_IMAGES, glPixelStorei, GL_PACK_SKIP_IMAGES);
	GEN_INT_SPEC_PARAM(GL_PACK_SKIP_PIXELS, glPixelStorei, GL_PACK_SKIP_PIXELS);
	GEN_INT_SPEC_PARAM(GL_PACK_SKIP_ROWS, glPixelStorei, GL_PACK_SKIP_ROWS);
	GEN_INT_SPEC_PARAM(GL_PACK_SWAP_BYTES, glPixelStorei, GL_PACK_SWAP_BYTES);
	GEN_INT_SPEC_PARAM(GL_PIXEL_PACK_BUFFER_BINDING, glBindBuffer, GL_PIXEL_PACK_BUFFER);
	GEN_INT_SPEC_PARAM(GL_PIXEL_UNPACK_BUFFER_BINDING, glBindBuffer, GL_PIXEL_UNPACK_BUFFER);
	GEN_INT_SPEC_PARAM(GL_POLYGON_SMOOTH_HINT, glHint, GL_POLYGON_SMOOTH_HINT);
	GEN_INT_SPEC_PARAM(GL_RENDERBUFFER_BINDING, glBindRenderbuffer, GL_RENDERBUFFER);
	GEN_INT_SPEC_PARAM(GL_SHADER_STORAGE_BUFFER_BINDING, glBindBuffer, GL_SHADER_STORAGE_BUFFER);
	GEN_INT_SPEC_PARAM(GL_STENCIL_BACK_WRITEMASK, glStencilMaskSeparate, GL_BACK);
	GEN_INT_SPEC_PARAM(GL_TEXTURE_BINDING_1D, glBindTexture, GL_TEXTURE_1D);
	GEN_INT_SPEC_PARAM(GL_TEXTURE_BINDING_1D_ARRAY, glBindTexture, GL_TEXTURE_1D_ARRAY);
	GEN_INT_SPEC_PARAM(GL_TEXTURE_BINDING_2D, glBindTexture, GL_TEXTURE_2D);
	GEN_INT_SPEC_PARAM(GL_TEXTURE_BINDING_2D_ARRAY, glBindTexture, GL_TEXTURE_2D_ARRAY);
	GEN_INT_SPEC_PARAM(GL_TEXTURE_BINDING_2D_MULTISAMPLE, glBindTexture, GL_TEXTURE_2D_MULTISAMPLE);
	GEN_INT_SPEC_PARAM(GL_TEXTURE_BINDING_2D_MULTISAMPLE_ARRAY, glBindTexture, GL_TEXTURE_2D_MULTISAMPLE_ARRAY);
	GEN_INT_SPEC_PARAM(GL_TEXTURE_BINDING_3D, glBindTexture, GL_TEXTURE_3D);
	GEN_INT_SPEC_PARAM(GL_TEXTURE_BINDING_BUFFER, glBindTexture, GL_TEXTURE_BUFFER);
	GEN_INT_SPEC_PARAM(GL_TEXTURE_BINDING_CUBE_MAP, glBindTexture, GL_TEXTURE_CUBE_MAP);
	GEN_INT_SPEC_PARAM(GL_TEXTURE_BINDING_RECTANGLE, glBindTexture, GL_TEXTURE_RECTANGLE);
	GEN_INT_SPEC_PARAM(GL_TEXTURE_COMPRESSION_HINT, glHint, GL_TEXTURE_COMPRESSION_HINT);
	GEN_INT_SPEC_PARAM(GL_TRANSFORM_FEEDBACK_BUFFER_BINDING, glBindBuffer, GL_TRANSFORM_FEEDBACK_BUFFER);
	GEN_INT_SPEC_PARAM(GL_UNIFORM_BUFFER_BINDING, glBindBuffer, GL_UNIFORM_BUFFER);
	GEN_INT_SPEC_PARAM(GL_UNPACK_ALIGNMENT, glPixelStorei, GL_UNPACK_ALIGNMENT);
	GEN_INT_SPEC_PARAM(GL_UNPACK_IMAGE_HEIGHT, glPixelStorei, GL_UNPACK_IMAGE_HEIGHT);
	GEN_INT_SPEC_PARAM(GL_UNPACK_LSB_FIRST, glPixelStorei, GL_UNPACK_LSB_FIRST);
	GEN_INT_SPEC_PARAM(GL_UNPACK_ROW_LENGTH, glPixelStorei, GL_UNPACK_ROW_LENGTH);
	GEN_INT_SPEC_PARAM(GL_UNPACK_SKIP_IMAGES, glPixelStorei, GL_UNPACK_SKIP_IMAGES);
	GEN_INT_SPEC_PARAM(GL_UNPACK_SKIP_PIXELS, glPixelStorei, GL_UNPACK_SKIP_PIXELS);
	GEN_INT_SPEC_PARAM(GL_UNPACK_SKIP_ROWS, glPixelStorei, GL_UNPACK_SKIP_ROWS);
	GEN_INT_SPEC_PARAM(GL_UNPACK_SWAP_BYTES, glPixelStorei, GL_UNPACK_SWAP_BYTES);
	GEN_INT_SPEC_PARAM(GL_POLYGON_MODE, glPolygonMode, GL_FRONT_AND_BACK);

	GEN_FLOAT_SPEC(GL_DEPTH_CLEAR_VALUE, glClearDepthf);
	GEN_FLOAT_SPEC(GL_LINE_WIDTH, glLineWidth);
	GEN_FLOAT_SPEC(GL_POINT_SIZE, glPointSize);

	GEN_VEC2_SPEC(GL_DEPTH_RANGE, glDepthRangef);

	GEN_VEC4_SPEC(GL_BLEND_COLOR, glBlendColor);
	GEN_VEC4_SPEC(GL_COLOR_CLEAR_VALUE, glClearColor);

	GEN_IVEC4_SPEC(GL_SCISSOR_BOX, glScissor);
	GEN_IVEC4_SPEC(GL_VIEWPORT, glViewport);

	GEN_BVEC4_SPEC(GL_COLOR_WRITEMASK, glColorMask);

	inline constexpr GLenum GL_BLEND_ALPHA = GL_BLEND_SRC_ALPHA;

	GEN_BLEND_SPEC(GL_BLEND_ALPHA);
	GEN_BLEND_SPEC(GL_BLEND_DST_ALPHA);

	inline constexpr GLenum GL_STENCIL_OP = GL_STENCIL_PASS_DEPTH_PASS;

	GEN_STENCIL_OP_SPEC(GL_STENCIL_OP);
	GEN_STENCIL_OP_SPEC(GL_STENCIL_FAIL);
	GEN_STENCIL_OP_SPEC(GL_STENCIL_PASS_DEPTH_FAIL);

	GEN_STENCIL_FUNC_SPEC(GL_STENCIL_FUNC);
	GEN_STENCIL_FUNC_SPEC(GL_STENCIL_REF);
	GEN_STENCIL_FUNC_SPEC(GL_STENCIL_VALUE_MASK);
}