#pragma once

#include <vector>
#include <string>
#include <istream>
#include <map>

#include "graph.hpp"
#include "utils.hpp"

namespace sms {
class Io {
 public:
  SparseGraph ReadGraph(std::istream& in);
};
} // namespace sms
