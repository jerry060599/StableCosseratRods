#include "../includes/modules/KittenAssets.h"
#include "../includes/modules/KittenPreprocessor.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <ratio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <regex>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

namespace Kitten {
	map<string, void*> resources;

	void loadDirectory(path root) {
		if (root.is_relative())
			for (auto& p : recursive_directory_iterator(root)) {
				path f(p);
				if (is_regular_file(f))
					loadAsset(f);
			}
		else
			cout << "err: resource path must be relative" << endl;
	}

	void loadAsset(path path) {
		glGetError();
		string ext = path.extension().string();
		if (ext == ".png" || ext == ".jpg" || ext == ".jpeg")
			loadTexture(path);
		else if (ext == ".vert" || ext == ".frag" || ext == ".comp" || ext == ".geom" || ext == ".tesc" || ext == ".tese")
			loadShader(path);
		else if (ext == ".obj" || ext == ".fbx" || ext == ".ply")
			loadMesh(path);
		else if (ext == ".node" || ext == ".face" || ext == ".ele")
			loadTetgenMesh(path);
		else if (ext == ".ttf")
			loadFont(path);
		else if (ext == ".glsl" || ext == ".include" || ext == ".mtl") {
		}
		else if (ext == ".csv" || ext == ".txt" || ext == ".cfg" || ext == ".json") {
			cout << "asset: loading text " << path.string().c_str() << endl;
			resources[path.string()] = new string(loadText(path.string()));
		}
		else
			cout << "err: unknown asset type " << path << endl;
		unsigned int error = glGetError();
		if (error != GL_NO_ERROR)
			printf("err: internal GL error: %d\n", error);
	}

	void loadMesh(path path) {
		cout << "asset: loading model " << path.string().c_str() << endl;
		loadMeshFrom(path);
	}

	void loadTetgenMesh(path path) {
		// cout << "asset: loading tetgen " << path.string().c_str() << endl;
		loadTetgenFrom(path);
	}

	void loadTexture(path path) {
		if (resources.count(path.string()))
			return;
		string filename = path.filename().string();
		string parsedName;
		Tags tags;
		parseAssetTag(filename, parsedName, tags);
		printf("asset: loading image %s (%s)\n", path.string().c_str(), parsedName.c_str());

		parsedName = path.parent_path().string().append("\\").append(parsedName);

		int width, height, nrChannels;
		stbi_set_flip_vertically_on_load(1);
		unsigned char* data = stbi_load(path.string().c_str(), &width, &height, &nrChannels, 4);
		if (data) {
			unsigned int handle;
			glGenTextures(1, &handle);
			glBindTexture(GL_TEXTURE_2D, handle);

			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			if ((tags["point"]).x) {
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			}
			else {
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			}

			Texture* img = new Texture;
			img->glHandle = handle;
			img->width = width;
			img->deviceFormat = GL_RGBA;
			img->hostFormat = GL_RGBA8;
			img->hostDataType = GL_UNSIGNED_BYTE;
			img->height = height;
			img->ratio = float(width) / height;
			img->borders = tags["border"];
			img->rawData = data;

			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
			glGenerateMipmap(GL_TEXTURE_2D);

			resources[parsedName] = img;
			resources[path.string()] = img;
		}
		else
			cout << "err: failed to load " << path << endl;
	}

	void loadShader(path path) {
		string ext = path.extension().string();
		string name = path.filename().string();

		unsigned int handle;
		int type = 0;
		Shader* shader;
		string sName = path.string().substr(0, path.string().length() - 5).append(".glsl");

		printf("asset: loading shader %s (%s)\n", sName.c_str(), ext.c_str());

		if (resources.count(sName))
			shader = (Shader*)resources[sName];
		else {
			shader = new Shader;
			resources[sName] = shader;
		}

		if (ext == ".vert") {
			type = (int)ShaderType::VERT;
			compileShader(path.string(), GL_VERTEX_SHADER, &handle);
		}
		else if (ext == ".frag") {
			type = (int)ShaderType::FRAG;
			compileShader(path.string(), GL_FRAGMENT_SHADER, &handle);
		}
		else if (ext == ".geom") {
			type = (int)ShaderType::GEOM;
			compileShader(path.string(), GL_GEOMETRY_SHADER, &handle);
		}
		else if (ext == ".tesc") {
			type = (int)ShaderType::TESS;
			compileShader(path.string(), GL_TESS_CONTROL_SHADER, &handle);
		}
		else if (ext == ".tese") {
			type = (int)ShaderType::TESS;
			compileShader(path.string(), GL_TESS_EVALUATION_SHADER, &handle);
		}
		else if (ext == ".comp") {
			type = (int)ShaderType::COMP;
			compileShader(path.string(), GL_COMPUTE_SHADER, &handle);
		}
		if (type == 0) {
			printf("err: unknown shader type %s.", name.c_str());
			return;
		}

		shader->unlinkedHandles.push_back(handle);
		shader->type |= type;
	}

	// Returns all the caches in a directory
	void getCaches(path p, vector<path>& paths) {
		auto name = p.filename().string();
		// Iterate over all files in the directory
		for (const auto& entry : std::filesystem::directory_iterator(p.parent_path())) {
			if (entry.is_regular_file()) {
				auto ep = entry.path();
				auto cname = ep.filename().string();

				// Check for .tmp extension and matching name
				if (ep.has_extension() && ep.extension() == ".tmp" &&
					cname.compare(0, name.length(), name) == 0) {

					// Check for '_'
					if (cname.length() > name.length() && cname[name.length()] == '_') {
						// printf("Found %s\n", ep.string().c_str());
						paths.push_back(ep);
					}
				}
			}
		}
	}

	std::size_t getLastModifiedEpoch(const std::filesystem::path& path) {
		auto ftime = std::filesystem::last_write_time(path);
		auto sctp = std::chrono::time_point_cast<std::chrono::seconds>(ftime);
		auto epoch = sctp.time_since_epoch();
		return std::chrono::duration_cast<std::chrono::seconds>(epoch).count();
	}

	std::string size_t2Hex(std::size_t num) {
		std::stringstream stream;
		stream << std::hex << num;
		return stream.str();
	}

	std::filesystem::path getCache(std::string key, size_t hash, int numCache) {
		auto targetPath = path(key + "_" + size_t2Hex(hash) + ".tmp");
		auto target = targetPath.filename().string();

		// Search if the number of hash
		vector<path> paths;
		getCaches(targetPath.has_parent_path() ? key : "./" + key, paths);

		// Check if we have a matching hash
		for (auto& p : paths)
			if (p.filename().string() == target) {
				// Touch file
				auto now = std::filesystem::file_time_type::clock::now();
				std::filesystem::last_write_time(p, now);

				return p;
			}

		// Cache not found, delete the oldest one if we have too many
		if (paths.size() && paths.size() >= numCache) {
			// Find the oldest one and delete it
			path oldest = paths[0];
			std::size_t oldestTime = std::numeric_limits<size_t>::max();

			for (auto& p : paths) {
				auto time = getLastModifiedEpoch(p);
				if (time < oldestTime) {
					oldestTime = time;
					oldest = p;
				}
			}

			std::filesystem::remove(oldest);
		}

		return target;
	}
}