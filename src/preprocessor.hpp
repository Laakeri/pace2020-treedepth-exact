#pragma once

#include <vector>

#include "graph.hpp"

namespace triangulator {

std::vector<int> ColToPar(const SparseGraph& graph, const std::vector<int>& col);

class Preprocessor {
 public:
 	Preprocessor() : org_graph(0), org_graph2(0) {}
	SparseGraph Preprocess(SparseGraph graph);
	SparseGraph TamakiRules(SparseGraph graph, int k);
	std::vector<int> Reconstruct(std::vector<int> colors) const;
	SparseGraph org_graph;
	SparseGraph org_graph2;
 private:
	void ParseTrees(SparseGraph& graph);
	std::vector<int> SolveTree(int v, const SparseGraph& graph, const std::vector<int>& parent);
	void GetRvs(int v, int ml, const SparseGraph& graph, const std::vector<int>& parent, std::vector<std::vector<int>>& rvs);
	std::pair<int, int> CTPDFS(int v, int p, int de, const SparseGraph& graph, const std::vector<int>& parent);
	int ConstructTreePar(int v, int p, int de, const SparseGraph& graph, const std::vector<int>& parent);
	struct GTree {
		std::vector<int> dummys;
		std::vector<int> vertices;
		std::vector<std::vector<int>> rvs;
		std::vector<std::vector<int>> clqs;
	};
	int n_;
  std::vector<int> tree_lab_;
  std::vector<int> tree_par_;
  std::vector<int> vertex_map_;
  std::vector<GTree> trees_;

  std::vector<std::pair<int, std::vector<int>>> tamaki_elim_;
  std::vector<int> vertex_map2_;
};

} // namespace triangulator