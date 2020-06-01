#include "utils.hpp"

#include <chrono>
#include <cassert>

namespace triangulator {
int Log::log_level_ = 10000;
void Log::SetLogLevel(int lvl) {
  log_level_ = lvl;
}

Timer::Timer() {
  timing = false;
  elapsedTime = std::chrono::duration<double>(std::chrono::duration_values<double>::zero());
}

void Timer::start() {
  if (timing) return;
  timing = true;
  startTime = std::chrono::steady_clock::now();
}

void Timer::stop() {
  if (!timing) return;
  timing = false;
  std::chrono::time_point<std::chrono::steady_clock> endTime = std::chrono::steady_clock::now();
  elapsedTime += (endTime - startTime);
}

double Timer::get() {
  if (timing) {
    stop();
    std::chrono::duration<double> ret = elapsedTime;
    start();
    return ret.count();
  }
  else {
    return elapsedTime.count();
  }
}

void Timer::clear() {
  elapsedTime = std::chrono::duration<double>(std::chrono::duration_values<double>::zero());
}
} // namespace triangulator
