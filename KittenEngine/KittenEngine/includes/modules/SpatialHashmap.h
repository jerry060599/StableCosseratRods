#pragma once
// Jerry Hsu, 2022

#include <atomic>
#include <glm/glm.hpp>

namespace Kitten {
	using namespace glm;

	/// <summary>
	/// An atomically synchronized hashmap for spatial hashing.
	/// add() and getNeighbors() cannot be called concurrently. 
	/// However, they are each thread-safe when not mixed.
	/// </summary>
	/// <typeparam name="T"></typeparam>
	template <typename T>
	class SpatialHashmap {
	public:
		const size_t maxSize;
		const float cellSize;

	private:
		typedef struct {
			int order;
			std::atomic<uint32_t> hash;
		} key;

		const float invCellSize;
		key* keys;
		T* vals;

		ivec3 getCell(vec3 pos) {
			return glm::ceil(pos * invCellSize);
		}

		uint32_t getCellHash(ivec3 cell) {
			uint32_t hash = cell.x + cell.y * 257 + cell.z * 65537;
			hash ^= (uint32_t)cell.y + 0x9e3779b9 + (hash << 3) + (hash >> 1);
			hash ^= (uint32_t)cell.x + 0x9e3779b9 + (hash << 3) + (hash >> 1);
			return hash + (!hash);
		}

	public:
		void clear() {
#pragma omp parallel for schedule(static, 4096)
			for (long long i = 0; i < (long long)maxSize; i++)
				keys[i].hash = 0;
		}

		SpatialHashmap(const size_t maxCol, const float maxDiameter) :
			maxSize(2llu << (int)ceil(log2((double)maxCol))), cellSize(maxDiameter), invCellSize(1.f / maxDiameter) {
			keys = new key[maxSize];
			vals = new T[maxSize];
			clear();
		}

		~SpatialHashmap() {
			delete[] keys;
			delete[] vals;
		}

		void add(vec3 pos, T val) {
			uint32_t hash = getCellHash(getCell(pos));
			uint32_t cell = hash % maxSize;
			uint32_t stride = 1;

			int order = 0;
			uint32_t expected = 0;

			while (true) {
				if (keys[cell].hash.compare_exchange_strong(expected, hash)) {
					keys[cell].order = order;
					vals[cell] = val;
					break;
				}

				order++;
				if (expected != hash) stride += 2;
				cell = (cell + stride) % maxSize;
				expected = 0;
			}
		}

		struct neighborhoodIter {
			using iterator_category = std::forward_iterator_tag;
			using difference_type = std::ptrdiff_t;
			using value_type = T;
			using pointer = value_type*;
			using reference = value_type&;

		private:
			inline const static ivec3 offsets[27]{
				ivec3(-1, -1, -1), ivec3(-1, -1, 0), ivec3(-1, -1, 1),
				ivec3(-1, 0, -1), ivec3(-1, 0, 0), ivec3(-1, 0, 1),
				ivec3(-1, 1, -1), ivec3(-1, 1, 0), ivec3(-1, 1, 1),

				ivec3(0, -1, -1), ivec3(0, -1, 0), ivec3(0, -1, 1),
				ivec3(0, 0, -1), ivec3(0, 0, 0), ivec3(0, 0, 1),
				ivec3(0, 1, -1), ivec3(0, 1, 0), ivec3(0, 1, 1),

				ivec3(1, -1, -1), ivec3(1, -1, 0), ivec3(1, -1, 1),
				ivec3(1, 0, -1), ivec3(1, 0, 0), ivec3(1, 0, 1),
				ivec3(1, 1, -1), ivec3(1, 1, 0), ivec3(1, 1, 1)
			};

			SpatialHashmap* map;
			ivec3 center;

			uint32_t hash;
			uint32_t cell;
			uint32_t stride;
			int order;
			int index;

			neighborhoodIter(SpatialHashmap* map, ivec3 center) : map(map), center(center) {
				if (!map) return;

				index = order = 0;
				stride = 1;
				hash = map->getCellHash(center + offsets[index]);
				cell = hash - 1;
				++(*this);
			}
			friend class SpatialHashmap;

		public:
			reference operator*() const {
				return map->vals[cell];
			}

			pointer operator->() {
				return &map->vals[cell];
			}

			neighborhoodIter& operator++() {
				if (!map) return *this;
				do {
					cell = (cell + stride) % map->maxSize;

					// End of this cell
					while (map->keys[cell].hash == 0) {
						// End of all cells
						if (++index == 27) {
							map = nullptr;
							break;
						}

						// reset to this cell
						hash = map->getCellHash(center + offsets[index]);
						cell = hash % map->maxSize;
						stride = 1;
						order = 0;
					}
					if (!map) break;

					if (map->keys[cell].hash == hash) {
						if (map->keys[cell].order == order) break;
					}
					else stride += 2;

					order++;
				} while (true);
				order++;

				return *this;
			}

			neighborhoodIter operator++(int) { neighborhoodIter tmp = *this; ++(*this); return tmp; }
			friend bool operator== (const neighborhoodIter& a, const neighborhoodIter& b) { return a.map == b.map; };
			friend bool operator!= (const neighborhoodIter& a, const neighborhoodIter& b) { return a.map != b.map; };
		};

		neighborhoodIter getNeighbors(vec3 pos) {
			return neighborhoodIter(this, getCell(pos));
		}

		neighborhoodIter end() {
			return neighborhoodIter(nullptr, ivec3(0));
		}
	};

	inline void testSpatialHashmap() {
		const int N = 3;
		SpatialHashmap<int> map(256 * 256 * N, 1);

#pragma omp parallel for schedule(dynamic, 1024)
		for (int i = 0; i < 256 * 256; i++) {
			for (size_t j = 0; j < N; j++)
				map.add(vec3(i / 256, i % 256, 0) + 0.5f, 10 * i + j);
		}

		for (auto itr = map.getNeighbors(vec3(34, 65, 0)); itr != map.end(); ++itr) {
			int n = *itr / 10;
			int a = n / 256;
			int b = n % 256;
			printf("%d %d %d\n", a, b, *itr % 10);
		}
	}
}