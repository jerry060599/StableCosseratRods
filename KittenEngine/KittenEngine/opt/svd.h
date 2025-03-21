#pragma once

#include "../includes/modules/Common.h"

namespace Kitten {
	void svd(const mat3& m, mat3& u, vec3& s, mat3& v);
	void svd(mat2 m, mat2& U, vec2& sig, mat2& V);
	void polarDecomp(mat2 m, mat2& R, mat2& S);
}