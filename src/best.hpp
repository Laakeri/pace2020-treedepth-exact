#pragma once

#include "graph.hpp"

namespace sms {
namespace best {
void InitBest(const SparseGraph& graph);
int SetBest(const std::vector<int>& par, bool is_best);
void PrintBest();
} // namespace best
} // namespace sms
