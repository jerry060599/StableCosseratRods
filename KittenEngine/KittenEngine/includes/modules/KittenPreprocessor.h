#pragma once
// Jerry Hsu, 2021

#include <map>
#include <set>
#include <vector>
#include <string>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

namespace Kitten {
	using namespace std;
	using namespace glm;

	extern vector<string> includePaths;
	typedef map<string, ivec4> Tags;
	void parseAssetTag(string& ori, string& name, Tags& tags);
	bool endsWith(string const& str, string const& ending);
	string loadText(string name);
	string loadTextWithIncludes(string path);
	void printWithLineNumber(string str);
}