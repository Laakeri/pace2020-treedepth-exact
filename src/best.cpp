 #include "best.hpp"

#include <iostream>
#include <vector>
#include <cassert>

#include "graph.hpp"

using std::vector;
using std::cout;
using std::flush;
using std::max;

namespace triangulator {
namespace best {
namespace {
int GetDepth(int v, const vector<vector<int>>& tr) {
	int d = 1;
	for (int nv : tr[v]) {
		d = max(d, GetDepth(nv, tr)+1);
	}
	return d;
}
} // namespace
int best_comp[10101010];
int v_map[10101010];
int best_td;
int n;

void InitBest(const SparseGraph& graph) {
  best_td = graph.n();
  n = graph.n();
  for (int i=0;i<n;i++){
    best_comp[i] = i;
    v_map[i] = graph.MapBack(i);
    assert(v_map[i] == i+1);
  }
}

int SetBest(const vector<int>& par, bool is_best) {
	assert((int)par.size() == n);
	int root = -1;
	vector<vector<int>> tr(n);
	for (int i=0;i<n;i++){
		if (par[i] == -1){
			assert(root == -1);
			root = i;
		} else {
			assert(par[i] >= 0 && par[i] < n);
			tr[par[i]].push_back(i);
		}
	}
	assert(root >= 0);
	int td = GetDepth(root, tr);
	Log::Write(3, "setbest ", td);
	if (td < best_td) {
		best_td = td;
		for (int i=0;i<n;i++) {
			if (par[i] == -1) {
				best_comp[i] = 0;
			} else {
				best_comp[i] = v_map[par[i]];
			}
		}
	} else {
		assert(!is_best);
	}
	return td;
}

void PrintBest() {
	cout<<best_td<<'\n';
	for (int i=0;i<n;i++){
		cout<<best_comp[i]<<'\n';
	}
	cout<<flush;
}

} // namespace best
} // namespace triangulator