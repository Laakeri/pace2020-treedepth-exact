#include "io.hpp"

#include <vector>
#include <string>
#include <istream>
#include <sstream>
#include <algorithm>
#include <cassert>

#include "graph.hpp"
#include "utils.hpp"

namespace triangulator {

std::vector<std::string> GetTokens(std::string s) {
  std::stringstream ss;
  ss<<s;
  std::vector<std::string> tokens;
  while (ss>>s) {
    tokens.push_back(s);
  }
  return tokens;
}

SparseGraph Io::ReadGraph(std::istream& in) {
  std::ios_base::sync_with_stdio(0);
  std::cin.tie(0);
  std::vector<std::pair<int, int> > edges;
  std::string input_line;
  bool pace = false;
  int n = -1;
  int m = -1;
  while (std::getline(in, input_line)) {
    assert(input_line.size() > 0);
    auto tokens = GetTokens(input_line);
    if (tokens.size() > 0 && tokens[0] == "c") continue;
    if (tokens.size() == 4 && tokens[0] == "p" && tokens[1] == "tdp") {
      assert(!pace);
      pace = true;
      n = std::stoi(tokens[2]);
      m = std::stoi(tokens[3]);
    } else if (tokens.size() == 2 && pace) {
      assert(tokens[0] != tokens[1]);
      edges.push_back({std::stoi(tokens[0]), std::stoi(tokens[1])});
    }
  }
  SparseGraph graph(edges);
  assert(graph.n() == n && graph.m() == m);
  return graph;
}
} // namespace triangulator
