
#include "../includes/modules/Font.h"
#include <iostream>
#include "../includes/modules/KittenAssets.h"
#include "../includes/modules/KittenRendering.h"
#include <vector>
#include <freetype/freetype.h>

namespace Kitten {
	float defaultFontLoadRes = 256;
	FT_Library ftLib;
	ComputeBuffer* meshBuffer;
	constexpr int MAX_CHAR_PER_DRAW = 1024;
	vector<vec4> meshDataBuffer;

	Font::~Font() {
		if (atlas) delete atlas;
		if (uvBuffer) delete uvBuffer;
	}

	void initFreetype() {
		if (FT_Init_FreeType(&ftLib)) {
			std::cout << "error: Error initializing freetype" << std::endl;
			return;
		}
		meshBuffer = new ComputeBuffer(sizeof(vec4), MAX_CHAR_PER_DRAW, GL_DYNAMIC_DRAW);
	}

	vec2 Font::render(const char* buff, float fontSize, vec4 color, Kitten::Bound<2> b,
		TextWrap horWrap, TextWrap vertWrap, TextJustification horJust, TextJustification vertJust) {
		glTempVar<GL_BLEND_ALPHA> blend(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glTempVar<GL_CULL_FACE> cull(false);
		auto ret = renderInternal(buff, fontSize, color,
			Kitten::get<Shader>("KittenEngine\\shaders\\text.glsl"), b, horWrap, vertWrap, horJust, vertJust);
		return ret;
	}

	vec2 Font::renderScreenspace(const char* buff, float fontSize, vec4 color, Kitten::Bound<2> b,
		TextWrap horWrap, TextWrap vertWrap, TextJustification horJust, TextJustification vertJust) {
		glTempVar<GL_DEPTH_WRITEMASK> zwrite(false);
		glTempVar<GL_BLEND_ALPHA> blend(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glTempVar<GL_DEPTH_TEST> dtest(false);
		glTempVar<GL_CULL_FACE> cull(false);
		auto oldM = modelMat;
		float aspect = Kitten::getAspect();
		b.min.x *= aspect;
		b.max.x *= aspect;
		modelMat = glm::ortho(-aspect, aspect, -1.f, 1.f) * modelMat;	// Hacky way to get screenspace
		auto ret = renderInternal(buff, fontSize, color,
			Kitten::get<Shader>("KittenEngine\\shaders\\textScreenspace.glsl"), b, horWrap, vertWrap, horJust, vertJust);
		modelMat = oldM;
		return ret;
	}

	void horJustify(int start, Bound<2> b, TextJustification horJust, Font* font, float scale) {
		if (!meshDataBuffer.size()) return;
		vec4 last = meshDataBuffer.back();
		auto lastChar = font->matrix[(int)last.x];
		float lastX = last.z - lastChar.bearing.x * scale + lastChar.adv * scale;
		if (horJust == TextJustification::CENTER) {
			float diff = (b.max.x - lastX) * 0.5f;
			for (int i = start; i < meshDataBuffer.size(); i++)
				meshDataBuffer[i].z += diff;
		}
		else if (horJust == TextJustification::RIGHT) {
			float diff = b.max.x - lastX;
			for (int i = start; i < meshDataBuffer.size(); i++)
				meshDataBuffer[i].z += diff;
		}
	}

	vec2 Font::renderInternal(const char* buff, float fontSize, vec4 color, Shader* shader, Kitten::Bound<2> b,
		TextWrap horWrap, TextWrap vertWrap, TextJustification horJust, TextJustification vertJust) {
		fontSize /= renderedSize;
		float scale = fontSize / atlas->height;
		vec2 pos = vec2(b.min.x, b.max.y);
		pos.y -= fontSize;
		if (vertWrap == TextWrap::TRUNCATE && pos.y < b.min.y)
			return pos;

		// Very primitive type setting algorithm
		meshDataBuffer.clear();
		int buffIndex = 0;
		int rollback = -1;
		int dataRollback = 0;
		int nlRollback = 0;
		int numLines = 1;
		for (; buff[buffIndex]; buffIndex++) {
			char c = buff[buffIndex];

			if (c == ' ') { 	// Handle wrapping points
				rollback = buffIndex;
				dataRollback = meshDataBuffer.size();
			}
			else if (c == '\n') {				// Handle new lines
				horJustify(nlRollback, b, horJust, this, scale);
				pos.x = b.min.x;
				pos.y -= fontSize;
				rollback = -1;
				nlRollback = dataRollback = meshDataBuffer.size();
				numLines++;
				if (vertWrap == TextWrap::TRUNCATE && pos.y < b.min.y)
					break;
			}
			else if (c < ' ' || c > 127)		// Ignore special chars
				continue;

			auto cinfo = matrix[c];

			// Unknown char
			if (cinfo.atlasPos == -1)
				continue;

			vec4 textInfo;
			textInfo.x = c;
			textInfo.y = fontSize;
			textInfo.z = pos.x + cinfo.bearing.x * scale;
			textInfo.w = pos.y + cinfo.bearing.y * scale;
			pos.x += cinfo.adv * scale;

			if (c != ' ')
				if (pos.x > b.max.x) {								// Out of bounds
					if (horWrap == TextWrap::WRAP) {
						if (meshDataBuffer.size() == nlRollback)	// At the minimum we want one char per line
							meshDataBuffer.push_back(textInfo);
						else if (rollback < 0) 						// No wrap point we just go back one char
							buffIndex--;
						else if (rollback >= 0) {					// We have a wrap point	
							buffIndex = rollback;
							meshDataBuffer.resize(dataRollback);
						}

						// New line
						horJustify(nlRollback, b, horJust, this, scale);
						pos.x = b.min.x;
						pos.y -= fontSize;
						rollback = -1;
						nlRollback = dataRollback = meshDataBuffer.size();
						numLines++;
						if (vertWrap == TextWrap::TRUNCATE && pos.y < b.min.y)
							break;
					}
					else if (horWrap == TextWrap::NONE)
						meshDataBuffer.push_back(textInfo);
				}
				else meshDataBuffer.push_back(textInfo);
		}
		horJustify(nlRollback, b, horJust, this, scale);

		// Do vertical justification
		if (vertJust == TextJustification::CENTER) {
			float diff = (b.max.y - b.min.y - numLines * fontSize) * 0.5f;
			for (auto& m : meshDataBuffer)
				m.w -= diff;
		}
		else if (vertJust == TextJustification::BOTTOM) {
			float diff = b.max.y - b.min.y - numLines * fontSize;
			for (auto& m : meshDataBuffer)
				m.w -= diff;
		}

		//Draw in batches
		if (meshDataBuffer.size()) {
			startRenderMesh(mat4(1));
			shader->use();
			uvBuffer->bind(3);
			meshBuffer->bind(4);
			atlas->bind(7);
			shader->setFloat4("textColor", color);
			glBindVertexArray(defMesh->VAO);

			int i = 0;
			for (; i < (int)meshDataBuffer.size() - MAX_CHAR_PER_DRAW; i += MAX_CHAR_PER_DRAW) {
				meshBuffer->upload(&meshDataBuffer[i]);
				glDrawElementsInstanced(shader->drawMode(), (GLsizei)defMesh->indices.size(), GL_UNSIGNED_INT, 0, MAX_CHAR_PER_DRAW);
			}
			int leftOver = (int)meshDataBuffer.size() - i;
			meshBuffer->upload(&meshDataBuffer[i], leftOver);
			glDrawElementsInstanced(shader->drawMode(), (GLsizei)defMesh->indices.size(), GL_UNSIGNED_INT, 0, leftOver);

			shader->unuse();
		}

		return pos;
	}

	Font* loadFont(path path) {
		if (resources.count(path.string())) return (Font*)resources[path.string()];
		cout << "asset: loading font " << path.string().c_str() << endl;

		FT_Face face;
		int err;
		if (err = FT_New_Face(ftLib, path.string().c_str(), 0, &face)) {
			std::cout << "error loading: Unable to load font " << path << err << endl;
			return nullptr;
		}

		FT_Set_Pixel_Sizes(face, 0, defaultFontLoadRes);

		Font* font = new Font;
		font->renderedSize = defaultFontLoadRes;
		resources[path.string()] = font;

		int atlasWidth = 0;
		int maxHeight = 0;

		for (size_t i = 32; i < 128; i++) {
			if (FT_Load_Char(face, i, FT_LOAD_RENDER)) {
				printf("error: Unable to load ascii \"%c\"\n", i);
				continue;
			}
			font->matrix[i] = CharInfo{
				ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
				ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
				int(face->glyph->advance.x) >> 6,
				atlasWidth
			};
			maxHeight = maxHeight < int(face->glyph->bitmap.rows) ? face->glyph->bitmap.rows : maxHeight;
			atlasWidth += face->glyph->bitmap.width;
		}

		//Define tab as four spaces
		font->matrix['\t'] = CharInfo{
			ivec2(0, 0),
			ivec2(0, 0),
			font->matrix[' '].adv * 4,
			-1
		};

		font->atlas = new Texture(atlasWidth, maxHeight, GL_RED, GL_RED, GL_UNSIGNED_BYTE);

		vec4 buff[129];
		glBindTexture(GL_TEXTURE_2D, font->atlas->glHandle);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		for (size_t i = 32; i < 128; i++) {
			if (font->matrix[i].atlasPos == -1)
				continue;
			FT_Load_Char(face, i, FT_LOAD_RENDER);
			glTexSubImage2D(GL_TEXTURE_2D, 0, font->matrix[i].atlasPos, 0, face->glyph->bitmap.width, face->glyph->bitmap.rows, GL_RED, GL_UNSIGNED_BYTE, face->glyph->bitmap.buffer);
			buff[i] = vec4(
				float(font->matrix[i].atlasPos) / atlasWidth,					// Position on atlas
				float(face->glyph->bitmap.width) / face->glyph->bitmap.rows,	// Aspect ratio
				face->glyph->bitmap.width / float(atlasWidth),					// Widtht in uv
				face->glyph->bitmap.rows / float(maxHeight));					// Height in uv
		}
		font->atlas->genMipmap();
		glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

		FT_Done_Face(face);

		buff[128].x = float(atlasWidth);
		buff[128].y = float(maxHeight);
		//Define a tab as four spaces
		font->matrix['\t'].atlasPos = 0;
		buff['\t'] = vec4(0, 0, 0, 0);

		font->uvBuffer = new UniformBuffer<vec4[129]>;
		font->uvBuffer->upload(buff);

		return font;
	}
}