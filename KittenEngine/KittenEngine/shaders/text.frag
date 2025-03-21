#version 430 core

#include "kittenCommonFrag.glsl"

out vec4 fragColor;

uniform vec4 textColor;
layout(binding = 7) uniform sampler2D atlas;

in vec2 uv;

void main() {
	fragColor = vec4(textColor.xyz, texture(atlas, uv).r);
	if (fragColor.w < 0.5)
		discard;
	fragColor.w *= textColor.w;
}