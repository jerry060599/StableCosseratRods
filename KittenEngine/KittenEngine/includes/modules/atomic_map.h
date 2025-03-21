#pragma once

// Jerry Hsu, 5/4/2021
// A simple atomic hash table implementation.

#include <atomic>
#include <utility>
#include <functional>
#include <algorithm>
#include <thread>

template<typename TKey, typename TVal>
class atomic_map {
private:
	enum class cell_state { FREE = 0, WRITE, READ, ABORTED };

	struct entry {
		std::atomic<size_t> hash;
		size_t idx;
		std::atomic<int> rev;
		std::atomic<int> lock;
	};

	size_t maxSize;
	std::atomic<size_t> c;
	std::pair<TKey, TVal>* list;
	entry* table;
	size_t tableSize;

	void clearTable() {
#pragma omp parallel for schedule(static, 1024)
		for (int64_t i = 0; i < (int64_t)tableSize; i++) {
			table[i].hash = 0;
			table[i].rev = 0;
			table[i].lock = (int)cell_state::FREE; // 0 = unallocated, 1 = writing, 2 = good data, 3 = aborted due to capcity limits.
			table[i].idx = 0;
		}
	}

	inline static const double MAX_LOAD_FACTOR = 0.5;
	inline static const size_t XOR_MASK = 0x3b6c307e19239b62;

public:
	size_t internalHash(TKey key) const {
		size_t hash = std::hash<TKey>()(key);
		hash ^= XOR_MASK;
		return hash ? hash : 1;
	}

	typedef void* handle;

	TVal read(handle handle) const {
		auto& ent = *(entry*)handle;
		TVal val;

		while (true) {
			auto prevRev = ent.rev.load();
			val = list[ent.idx].second;

			if (ent.lock.load() != (int)cell_state::READ) { // Make sure some write didnt start.
				std::this_thread::yield();
				continue;
			}
			if (prevRev == ent.rev.load()) break;
		}
		return val;
	}

	void write(handle handle, TVal& val) {
		auto& ent = *(entry*)handle;
		// Get write lock
		int lock = (int)cell_state::READ;

		// Get write lock
		while (!ent.lock.compare_exchange_strong(lock, (int)cell_state::WRITE)) lock = (int)cell_state::READ;
		list[ent.idx].second = val;
		ent.rev.fetch_add(1); // Invalidate all concurrent reads
		ent.lock.store((int)cell_state::READ); // release lock
	}

	void map(handle handle, std::function<bool(TVal&, bool)> mapping) {
		auto& ent = *(entry*)handle;
		// Get write lock
		int lock = (int)cell_state::READ;

		// Get write lock
		while (!ent.lock.compare_exchange_strong(lock, (int)cell_state::WRITE)) lock = (int)cell_state::READ;

		if (mapping(list[ent.idx].second, true)) ent.rev.fetch_add(1); // Invalidate all concurrent reads

		ent.lock.store((int)cell_state::READ); // release lock
	}

	handle getHandleWithMap(TKey& key, std::function<bool(TVal&, bool)> mapping, size_t hash = 0) {
		if (!hash) hash = internalHash(key);

		size_t idx = hash % tableSize;
		while (true) {
			size_t old = 0;
			auto& ent = table[idx];
			bool res = ent.hash.compare_exchange_strong(old, hash);
			if (res) { // We sucsessfully allocated a new entry
				// Write data
				size_t listIndex = c.fetch_add(1);
				if (listIndex >= maxSize) { // Abort the write
					c.fetch_sub(1);
					ent.hash.store(0);
					ent.lock.store((int)cell_state::ABORTED);
					return nullptr;
				}

				TVal val = {};
				mapping(val, false);
				list[listIndex] = std::make_pair(key, val);

				// Write table
				ent.idx = listIndex;
				ent.lock.store((int)cell_state::READ);
				return &table[idx];
			}
			else if (old == hash) {  // We hit a possible existing entry
				int lVal;
				while (!(lVal = ent.lock.load())); // Spinlock until this entry is consistent (lock == 0 means that its still allocating)

				if (lVal == (int)cell_state::ABORTED) return nullptr;
				if (list[ent.idx].first == key) { // Make sure both the hash and key match
					handle handle = &table[idx];
					map(handle, mapping);
					return handle;
				}
			}
			idx = (idx + 1) % tableSize;
		}

		return nullptr;
	}

	handle getHandle(TKey& key, size_t hash = 0) const {
		if (!hash) hash = internalHash(key);

		size_t idx = hash % tableSize;
		while (true) {
			auto& ent = table[idx];
			size_t cell = ent.hash.load();
			if (cell == 0) return nullptr;
			if (cell == hash) {  // We hit a possible existing entry
				int lVal;
				while (!(lVal = ent.lock.load())); // Spinlock until this entry is consistent (lock == 0 means that its still allocating)

				if (lVal == (int)cell_state::ABORTED) return nullptr;
				if (list[ent.idx].first == key)  // Make sure both the hash and key match
					return &table[idx];
			}
			idx = (idx + 1) % tableSize;
		}

		return nullptr;
	}

	/// <summary>
	/// The number of entries in this table
	/// </summary>
	/// <returns></returns>
	size_t size() const {
		return std::min(c.load(), maxSize);
	}

	/// <summary>
	/// The max number of entries in this table
	/// </summary>
	/// <returns></returns>
	size_t max_size() const {
		return maxSize;
	}

	/// <summary>
	/// Atomically maps an old value to a new one with a custom mapping. 
	/// </summary>
	/// <param name="key">key</param>
	/// <param name="mapping">a function that takes a value reference and whether it is an existing value. should return true if ref is modified.</param>
	/// <param name="hash">optional hash where the hash can be reused</param>
	/// <returns>whether the value existed/was added to this table</returns>
	bool map(TKey& key, std::function<bool(TVal&, bool)> mapping, size_t hash = 0) {
		return getHandleWithMap(key, mapping, hash);
	}

	/// <summary>
	/// Inserts an element into this table
	/// </summary>
	/// <param name="key">key</param>
	/// <param name="val">value</param>
	/// <param name="hash">optional hash where the hash can be reused</param>
	/// <returns>whether the value existed/was added to this table</returns>
	bool insert(TKey& key, TVal& val, size_t hash = 0) {
		return getHandleWithMap(key, [&val](TVal& v, bool init)->bool { v = val; return true; }, hash);
	}

	/// <summary>
	/// Gets an entry from this table
	/// </summary>
	/// <param name="key">key</param>
	/// <param name="val">value reference to be stored in</param>
	/// <param name="hash">optional hash where the hash can be reused</param>
	/// <returns>whether the entry was found</returns>
	bool get(TKey& key, TVal& val, size_t hash = 0) const {
		handle handle = getHandle(key, hash);
		if (handle) {
			val = read(handle);
			return true;
		}
		return false;
	}

	/// <summary>
	/// Culls the entries of this table in place. NOT ATOMIC/THREAD SAFE!!
	/// </summary>
	/// <param name="predicate">a predicate determining whether to keep this value or not</param>
	void cull(std::function<bool(const TKey&, TVal&)> predicate) {
		int64_t s = (int64_t)size();
		c.store(0);
		clearTable();
		for (int64_t i = 0; i < s; i++)
			if (predicate(list[i].first, list[i].second))
				insert(list[i].first, list[i].second);
	}

	/// <summary>
	/// Resizes this table to a new max size. NOT ATOMIC/THREAD SAFE!!
	/// </summary>
	/// <param name="newSize">the new size</param>
	void resize(size_t newSize) {
		if (newSize == maxSize) return;

		auto oldList = list;
		int64_t s = (int64_t)size();

		delete[] table;
		c.store(0);
		maxSize = newSize;
		tableSize = (size_t)std::ceil(maxSize / MAX_LOAD_FACTOR);
		list = new std::pair<TKey, TVal>[maxSize];
		table = new entry[tableSize];
		clearTable();

#pragma omp parallel for schedule(static, 512)
		for (int64_t i = 0, s = (int64_t)size(); i < s; i++)
			insert(oldList[i].first, oldList[i].second);

		delete[] oldList;
	}

	/// <summary>
	/// Creates a table
	/// </summary>
	/// <param name="size">maximum number of entries</param>
	atomic_map(size_t size) : maxSize(size), tableSize((size_t)std::ceil(size / MAX_LOAD_FACTOR)) {
		c = 0;
		list = new std::pair<TKey, TVal>[maxSize];
		table = new entry[tableSize];
		clearTable();
	}

	~atomic_map() {
		delete[] table;
		delete[] list;
	}
};

inline bool AtomicTableValidate() {
	bool good = true;

	{ // Simple parallel insertion test
		const int N = 1024;
		atomic_map<int, int> table(N);
#pragma omp parallel for schedule(static, 1)
		for (int i = 0; i < N; i++) {
			int key = 3 * i + 1;
			int val = key * key - key;
			table.insert(key, val);
		}

		for (int i = 0; i < N; i++) {
			int key = 3 * i + 1;
			int val = key * key - key;
			int aval;
			bool ret = table.get(key, aval);
			if (ret) {
				if (aval != val) {
					printf("err: wrong value found. 1\n");
					good = false;
				}
			}
			else {
				printf("err: expected value. found none.\n");
				good = false;
			}
		}

		for (int i = 0; i < N; i++) {
			int key = 3 * i + 2;
			int val = key * key - key;
			int aval;
			bool ret = table.get(key, aval);
			if (ret) {
				printf("err: expected no value but value exists.\n");
				good = false;
			}
		}
	}

	{ // Write self consistent test
		// Makes sure that large writes get correctly written atomically (i.e. all parts of TVal agree)
		struct testStruct {
			int data[64];
		};

		atomic_map<int, testStruct> table(16);

		const int N = 1024 * 512;
#pragma omp parallel for schedule(static, 1)
		for (int i = 0; i < N; i++) {
			int key = i % 2;
			int r = rand() % 100000;
			testStruct val;
			for (int j = 0; j < 64; j++)
				val.data[j] = j * r + j;
			table.insert(key, val);
		}

		for (int i = 0; i < 2; i++) {
			int key = i;
			testStruct val;

			bool ret = table.get(key, val);
			if (ret) {
				int r = val.data[1] - 1;
				for (int j = 0; j < 64; j++)
					if (val.data[j] != j * r + j) {
						printf("err: wrong value found. 2\n");
						good = false;
					}
			}
			else {
				printf("err: expected value. found none.\n");
				good = false;
			}
		}
	}

	{ // replace atomicity test
		const int N = 1024 * 512;

		atomic_map<int, int> table(16);

#pragma omp parallel for schedule(static, 1)
		for (int i = 0; i < N; i++) {
			int key = i % 2;
			table.map(key, [&i](int& v, bool init)->bool {
				if (init) {
					if (v < i) { // A value already exists
						v = i;
						return true;
					}
				}
				else  // No value found.
					v = i;
				return false;
				});
		}

		for (int i = 0; i < 2; i++) {
			int val;
			bool ret = table.get(i, val);
			if (ret) {
				if (val != N - 2 + i) {
					printf("err: wrong value found %d. 3\n", val);
					good = false;
				}
			}
			else {
				printf("err: expected value. found none.\n");
				good = false;
			}
		}
	}

	{ // concurrent read/write test
		struct testStruct {
			int data[64];
		};

		atomic_map<int, testStruct> table(16);

		const int N = 1024 * 512;
#pragma omp parallel for schedule(static, 1)
		for (int i = 0; i < N; i++) {
			bool testWrite = rand() % 2;
			int key = i % 2;
			if (testWrite) {
				int r = rand() % 100000;
				testStruct val;
				for (int j = 0; j < 64; j++)
					val.data[j] = j * r + j;
				table.insert(key, val);
			}
			else {
				testStruct val;
				bool ret = table.get(key, val);
				if (ret) {
					int r = val.data[1] - 1;
					for (int j = 0; j < 64; j++)
						if (val.data[j] != j * r + j) {
							printf("err: wrong value found. 4\n");
							good = false;
						}
				}
			}
		}
	}

	printf("table valid: %s\n", good ? "true" : "false");
	return good;
}