#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <queue>
#include <random>

#include "bitset.hpp"

namespace triangulator {
// Interface

namespace utils {
template<typename T>
void SortAndDedup(std::vector<T>& vec);

template<typename T>
void InitZero(std::vector<T>& vec, size_t size);

template<typename T>
std::vector<T> PermInverse(const std::vector<T>& perm);

template<typename T>
T GetRand(T a, T b, std::mt19937& gen);
} // namespace utils

class Log {
 public:
  template<typename T, typename... Args>
  static void P(T first_message, Args... message);
  
  template<typename T>
  static void Write(int lvl, T message);

  template<typename T, typename... Args>
  static void Write(int lvl, T first_message, Args... message);

  static void SetLogLevel(int lvl);
 private:
  template<typename T>
  static void WriteImpl(std::vector<T> message);

  template<typename T>
  static void WriteImpl(std::pair<T, T> message);

  template<size_t chunks>
  static void WriteImpl(FBitset<chunks> message);
  
  template<typename T>
  static void WriteImpl(T message);
  static int log_level_;
};

class Timer {
 private:
  bool timing;
  std::chrono::duration<double> elapsedTime;
  std::chrono::time_point<std::chrono::steady_clock> startTime;
 public:
  Timer();
  void start();
  void stop();
  void clear();
  double get();
};

class PolyHash {
 private:
  uint64_t val_;
 public:
  PolyHash() {
    val_ = 0;  
  }
  void Add(uint64_t n) {
    n %= 1000000007;
    val_ *= 65599;
    val_ += n + 59;
    val_ %= 1000000007;
  }
  uint64_t Value() const {
    return val_;
  }
};

// Implementation

namespace utils {

template<typename T>
T GetRand(T a, T b, std::mt19937& gen) {
  return std::uniform_int_distribution<T>(a,b)(gen);
}

template<typename T>
void SortAndDedup(std::vector<T>& vec) {
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

template<typename T>
void InitZero(std::vector<T>& vec, size_t size) {
  vec.resize(size);
  std::fill(vec.begin(), vec.begin() + size, 0);
}

template<typename T>
std::vector<T> PermInverse(const std::vector<T>& perm) {
  std::vector<T> ret(perm.size());
  for (int i = 0; i < perm.size(); i++) {
    ret[perm[i]] = i;
  }
  return ret;
}

inline int GetU(int x, std::vector<int>& un) {
  if (un[x] == x) return x;
  else {
    un[x] = GetU(un[x], un);
    return un[x];
  }
}
} // namespace utils
template<typename T>
void Log::WriteImpl(T message) {
  std::cerr<<message;
}
template<>
inline void Log::WriteImpl(std::vector<char> message) {
  std::cerr<<"{ ";
  for (char e : message) {
    std::cerr<<(int)e<<" ";
  }
  std::cerr<<"}";
}
template<>
inline void Log::WriteImpl(std::vector<int> message) {
  std::cerr<<"{";
  for (int i=0;i<(int)message.size();i++){
    std::cerr<<message[i];
    if (i+1<(int)message.size()) {
      std::cerr<<", ";
    }
  }
  std::cerr<<"}";
}
template<typename T>
inline void Log::WriteImpl(std::pair<T, T> message) {
  std::cerr<<"{";
  WriteImpl(message.first);
  std::cerr<<", ";
  WriteImpl(message.second);
  std::cerr<<"}";
}
template<typename T>
inline void Log::WriteImpl(std::vector<T> message) {
  std::cerr<<"{";
  for (int i=0;i<(int)message.size();i++){
    WriteImpl(message[i]);
    if (i+1<(int)message.size()) {
      std::cerr<<", ";
    }
  }
  std::cerr<<"}";
}
template<size_t chunks>
inline void Log::WriteImpl(FBitset<chunks> b) {
  std::cerr<<"{";
  for (size_t i=0;i<chunks*BITS;i++){
    std::cerr<<b.Get(i);
  }
  std::cerr<<"}";
}
template<typename T>
void Log::Write(int lvl, T message) {
  if (lvl > log_level_) return;
  WriteImpl(message);
  std::cerr<<std::endl;
}
template<typename T, typename... Args>
void Log::Write(int lvl, T first_message, Args... message) {
  if (lvl > log_level_) return;
  WriteImpl(first_message);
  Write(lvl, message...);
}
template<typename T, typename... Args>
void Log::P(T first_message, Args... message) {
  Write(0, first_message, message...);
}
} // namespace triangulator
