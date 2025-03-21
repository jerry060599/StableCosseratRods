#pragma once
// Jerry Hsu, 2021

#include <string>
#include <vector>
#include <chrono>

namespace Kitten {
	using namespace std::chrono;

	/// <summary>
	/// A timer that keeps track of the timeline
	/// </summary>
	class StopWatch {
	public:
		bool gpuSync = true;

		steady_clock::time_point lastPoint;
		std::vector<double> deltaTimes;
		std::vector<double> times;
		std::vector<const char*> tags;
	private:
		double totSecs = 0;

	public:
		StopWatch();
		double totTime();
		void printTimes();
		double time(const char* tag = nullptr);
		void reset();
	};
}