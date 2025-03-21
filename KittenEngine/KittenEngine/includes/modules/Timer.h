#pragma once
// Jerry Hsu, 2021

#include <string>
#include <vector>
#include <chrono>
#include <unordered_map>
#include "Dist.h"

namespace Kitten {
	using namespace std::chrono;

	/// <summary>
	/// A timer that keeps track of the distribution instead of the timeline
	/// </summary>
	class Timer {
	public:
		struct entry {
			Dist dist;
			steady_clock::time_point lastPoint;
			bool inFence = false;
			bool gpuSync = true;
		};

	private:
		double totSecs = 0;
		std::unordered_map<const char*, entry> entries;

	public:
		Timer();
		void printTimes();
		void start(const char* tag = nullptr, bool gpuSync = true);
		double end(const char* tag = nullptr);
		void reset();
	};
}