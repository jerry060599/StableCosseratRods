#pragma once
// Jerry Hsu, 2022

#include <atomic>
#include <stdio.h>

namespace Kitten {
	template <typename T>
	class UniqueList {
	private:
		std::atomic<size_t> count;

		typedef struct {
			std::atomic<uint64_t> key;
			size_t index;
		} entry;

		const size_t tableSize;
		entry* table;
		T* data;

	public:
		UniqueList(const size_t capacity) :
			tableSize(2llu << (int)ceil(log2((double)capacity))) {
			data = new T[capacity];
			table = new entry[tableSize];
			clear();
		}

		~UniqueList() {
			delete[] data;
			delete[] table;
		}

		uint64_t scramble(uint64_t x) {
			uint64_t a = (x + 11 * (x >> 32)) ^ 0x9e3779b9llu;
			if (a == 0xFFFFFFFFllu) return (0xFFFFFFFFllu + 11 * (0xFFFFFFFFllu >> 32)) ^ 0x9e3779b9llu;
			return a;
		}

		void clear() {
			count = 0;
#pragma omp parallel for schedule(static, 4096)
			for (long long i = 0; i < (long long)tableSize; i++)
				table[i].key = 0xFFFFFFFFllu;
		}

		bool add(uint64_t key, T val) {
			key = scramble(key);
			uint64_t cell = key % tableSize;
			uint64_t stride = 1;

			uint64_t expected = 0xFFFFFFFFllu;

			while (true) {
				if (table[cell].key.compare_exchange_strong(expected, key)) {
					size_t n = count.fetch_add(1);
					table[cell].index = n;
					data[n] = val;
					return true;
				}
				if (expected == key)
					return false;

				cell = (cell + (stride++)) % tableSize;
				expected = 0xFFFFFFFFllu;
			}
		}

		T* find(uint64_t key) {
			key = scramble(key);
			uint64_t cell = key % tableSize;
			uint64_t stride = 1;

			uint64_t expected = 0;

			while (true) {
				size_t v = table[cell].key;
				if (v == key)
					return &data[table[cell].index];
				else if (v == 0xFFFFFFFFllu)
					return nullptr;

				cell = (cell + (stride++)) % tableSize;
				expected = 0xFFFFFFFFllu;
			}
		}

		bool contains(uint64_t key) {
			return find(key) != nullptr;
		}

		T& operator[](size_t i) {
			return data[i];
		}

		size_t size() {
			return count;
		}
	};

	inline void testUniqueList() {
		UniqueList<int> list(1024);

		list.add(0, 1);
		list.add(3, 2);
		list.add(15, 3);
		list.add(3, 4);
		list.add(16, 5);

		for (size_t i = 0; i < list.size(); i++)
			printf("%d\n", list[i]);
	}
}