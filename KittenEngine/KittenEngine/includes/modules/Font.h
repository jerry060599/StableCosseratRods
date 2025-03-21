#pragma once
// Jerry Hsu, 2024

#include <filesystem>
#include "Texture.h"
#include "Shader.h"
#include "UniformBuffer.h"
#include "Bound.h"

namespace Kitten {
	using namespace std::filesystem;
	extern float defaultFontLoadRes;

	enum class TextWrap {
		NONE,
		TRUNCATE,
		WRAP
	};

	enum class TextJustification {
		CENTER,
		LEFT,
		RIGHT,
		TOP,
		BOTTOM
	};

	struct CharInfo {
		ivec2 size;			// Size of the glyph image
		ivec2 bearing;		// Glyph bearing
		int adv;			// How much to advance horizontally after this character.
		int atlasPos = -1;	// Horizontal position in the atlas
	};

	class Font {
	public:
		float renderedSize = 0;
		Texture* atlas = nullptr;
		UniformBuffer<vec4[129]>* uvBuffer = nullptr;
		ivec2 atlasSize;
		CharInfo matrix[129];

		~Font();

		vec2 renderScreenspace(const char* buff, float fontSize, vec4 color, Kitten::Bound<2> b,
			TextWrap horWrap = TextWrap::NONE, TextWrap vertWrap = TextWrap::NONE,
			TextJustification horJust = TextJustification::LEFT, TextJustification vertJust = TextJustification::TOP);
		vec2 render(const char* buff, float fontSize, vec4 color, Kitten::Bound<2> b,
			TextWrap horWrap = TextWrap::NONE, TextWrap vertWrap = TextWrap::NONE,
			TextJustification horJust = TextJustification::LEFT, TextJustification vertJust = TextJustification::TOP);
	private:
		vec2 renderInternal(const char* buff, float fontSize, vec4 color, Shader* shader, Kitten::Bound<2> b,
			TextWrap horWrap = TextWrap::NONE, TextWrap vertWrap = TextWrap::NONE,
			TextJustification horJust = TextJustification::LEFT, TextJustification vertJust = TextJustification::TOP);
	};

	void initFreetype();
	Font* loadFont(path path);
}