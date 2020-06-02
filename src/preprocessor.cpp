#include "preprocessor.hpp"

#include "graph.hpp"
#include "utils.hpp"

#define F first
#define S second

using std::vector;
using std::queue;
using std::max;
using std::pair;

namespace sms {
namespace {
int SubtreeSize(int v, const SparseGraph& graph, const vector<int>& parent) {
  int sz = 1;
  for (int nv : graph.Neighbors(v)) {
    if (parent[nv] == v) {
      sz += SubtreeSize(nv, graph, parent);
    }
  }
  return sz;
}

void DelTree(int v, SparseGraph& graph, const vector<int>& parent, vector<int>& isols) {
  auto nbs = graph.Neighbors(v);
  for (int nv : nbs) {
    if (parent[nv] == v) {
      DelTree(nv, graph, parent, isols);
      graph.RemoveEdge(v, nv);
      assert(graph.Degree(nv) == 0);
      isols.push_back(nv);
    }
  }
}
} // namespace

vector<int> Preprocessor::SolveTree(int v, const SparseGraph& graph, const vector<int>& parent) {
  vector<int> rank = {1};
  int root = 0;
  for (int nv : graph.Neighbors(v)) {
    if (parent[nv] == v) {
      auto sr = SolveTree(nv, graph, parent);
      rank.resize(max(sr.size(), rank.size())+1);
      for (int i=root;i<(int)sr.size();i++){
        assert(sr[i] <= 1);
        rank[i] += sr[i];
      }
      for (int i=root;i<(int)rank.size();i++){
        assert(rank[i] <= 3);
        if (rank[i] >= 2) {
          assert(i+1>root);
          assert(i+1<(int)rank.size());
          root=i+1;
          rank[i] = 0;
          rank[i+1]++;
        }
        assert(rank[i] <= 1);
      }
      for (int i=0;i<root;i++) {
        rank[i] = 0;
      }
      while (rank.back() == 0){
        rank.pop_back();
      }
    }
  }
  assert(rank.size() > 0 && rank.back() == 1);
  tree_lab_[v] = root;
  assert(root < (int)rank.size());
  assert(rank[root] == 1);
  return rank;
}

void Preprocessor::GetRvs(int v, int ml, const SparseGraph& graph, const vector<int>& parent, vector<vector<int>>& rvs) {
  assert(tree_lab_[v] < (int)rvs.size());
  if (tree_lab_[v] > ml) {
    ml = tree_lab_[v];
    rvs[ml].push_back(v);
  }
  for (int nv : graph.Neighbors(v)) {
    if (parent[nv] == v) {
      GetRvs(nv, ml, graph, parent, rvs);
    }
  }
}

pair<int, int> Preprocessor::CTPDFS(int v, int p, int de, const SparseGraph& graph, const vector<int>& parent) {
  assert(tree_par_[v] == -1);
  assert(tree_lab_[v] >= 0);
  assert(tree_lab_[v] <= de);
  pair<int, int> ml = {v, tree_lab_[v]};
  for (int nv : graph.Neighbors(v)) {
    if (parent[nv] == v || parent[v] == nv) {
      if (tree_lab_[nv] <= de && nv != p) {
        auto cml = CTPDFS(nv, v, de, graph, parent);
        if (cml.S > ml.S) {
          ml = cml;
        } else if (cml.S == ml.S) {
          ml.F = -1;
        }
      }
    }
  }
  return ml;
}

int Preprocessor::ConstructTreePar(int v, int p, int de, const SparseGraph& graph, const vector<int>& parent) {
  assert(tree_par_[v] == -1);
  assert(tree_lab_[v] >= 0);
  assert(tree_lab_[v] <= de);
  pair<int, int> r = CTPDFS(v, -1, de, graph, parent);
  assert(r.F >= -1 && r.S >= 0 && r.S <= de);
  v = r.F;
  de = r.S;
  assert(tree_par_[v] == -1);
  assert(tree_lab_[v] == de);
  if (p != -1) {
    tree_par_[v] = p;
  }
  for (int nv : graph.Neighbors(v)) {
    if (parent[nv] == v || parent[v] == nv) {
      if (tree_lab_[nv] < de) {
        int cr = ConstructTreePar(nv, v, de-1, graph, parent);
        assert(tree_par_[cr] == v);
        assert(tree_lab_[cr] < tree_lab_[v]);
      }
    }
  }
  return v;
}

void Preprocessor::ParseTrees(SparseGraph& graph) {
  int n = graph.n();
  vector<int> dgs(n);
  vector<int> parent(n);
  vector<int> isols;
  queue<int> proc;
  tree_par_.resize(n);
  tree_lab_.resize(n);
  for (int i=0;i<n;i++) {
    tree_par_[i] = -1;
    tree_lab_[i] = -1;
    parent[i] = -1;
    dgs[i] = graph.Degree(i);
    if (dgs[i] == 0) {
      isols.push_back(i);
    } else if (dgs[i] == 1) {
      proc.push(i);
    }
  }
  while (!proc.empty()) {
    int v = proc.front();
    proc.pop();
    int cnt = 0;
    for (int nv : graph.Neighbors(v)) {
      if (parent[nv] == -1) {
        cnt++;
      }
    }
    assert(cnt == dgs[v]);
    assert(dgs[v] <= 1);
    assert(parent[v] == -1);
    for (int nv : graph.Neighbors(v)) {
      if (parent[nv] == -1) {
        parent[v] = nv;
        dgs[nv]--;
        if (dgs[nv] == 1) {
          proc.push(nv);
        }
      } else {
        assert(parent[nv] == v);
      }
    }
  }
  for (int i=0;i<n;i++) {
    if (parent[i] == -1) {
      int sz = SubtreeSize(i, graph, parent);
      if (sz >= 3) {
        auto sol = SolveTree(i, graph, parent);
        vector<vector<int>> rvs(sol.size());
        GetRvs(i, -1, graph, parent, rvs);
        for (int j=0;j<(int)sol.size();j++){
          if (sol[j] == 0) {
            assert(rvs[j].size() == 0);
          } else {
            assert(rvs[j].size() > 0);
          }
        }
        int iss = isols.size();
        DelTree(i, graph, parent, isols);
        int v = i;
        vector<int> dummys;
        vector<int> tvs;
        vector<vector<int>> clqs;
        for (int j = iss; j < (int)isols.size(); j++) {
          tvs.push_back(isols[j]);
        }
        for (int j = 0; j < (int)sol.size(); j++) {
          assert(sol[j] >= 0 && sol[j] <= 1);
          if (sol[j] == 0) continue;
          vector<int> clq = {v};
          for (int jj = 0; jj < j; jj++) {
            assert(isols.size() > 0);
            clq.push_back(isols.back());
            dummys.push_back(isols.back());
            isols.pop_back();
          }
          assert((int)clq.size() == j+1);
          clqs.push_back(clq);
          for (int a : clq) {
            for (int b : clq) {
              if (a < b) {
                assert(!graph.HasEdge(a, b));
                graph.AddEdge(a, b);
              }
            }
          }
          if (j+1 < (int)sol.size()) {
            assert(isols.size() > 0);
            graph.AddEdge(v, isols.back());
            v = isols.back();
            dummys.push_back(isols.back());
            isols.pop_back();
          }
        }
        assert(tvs.size() >= dummys.size());
        trees_.push_back({dummys, tvs, rvs, clqs});
      }
    }
  }
}

SparseGraph Preprocessor::Preprocess(SparseGraph graph) {
  org_graph = graph;
  n_ = graph.n();
  ParseTrees(graph);
  SparseGraph ppg(graph.Edges());
  vertex_map_.resize(ppg.n());
  for (int i = 0; i < ppg.n(); i++) {
    vertex_map_[i] = ppg.MapBack(i);
  }
  return ppg;
}

SparseGraph Preprocessor::TamakiRules(SparseGraph graph, int k) {
  org_graph2 = graph;
  bool fo = true;
  while (fo) {
    fo = false;
    for (int x = 0; x < graph.n(); x++) {
      if (graph.Degree(x) < k) continue;
      for (int y = x+1; y < graph.n(); y++) {
        if (graph.Degree(y) < k) continue;
        if (graph.HasEdge(x, y)) continue;
        if (graph.Mincut(x, y) >= k) {
          graph.AddEdge(x, y);
          Log::Write(3, "cutrule ", x, " ", y);
          fo = true;
        }
      }
    }
  }
  vector<pair<int, int>> dgo;
  for (int i=0;i<graph.n();i++){
    dgo.push_back({graph.Degree(i), i});
  }
  sort(dgo.begin(), dgo.end());
  fo = true;
  while (fo) {
    fo = false;
    for (int i=0;i<graph.n();i++){
      int x = dgo[i].second;
      if (graph.Degree(x) == 0) continue;
      if (!graph.IsClique(graph.Neighbors(x))) continue;
      bool ok = true;
      for (int nx : graph.Neighbors(x)) {
        if (graph.Degree(nx) <= k) {
          ok = false;
          break;
        }
      }
      if (ok) {
        fo = true;
        Log::Write(3, "tamakirule ", x, " ", graph.Degree(x));
        auto nbs = graph.Neighbors(x);
        tamaki_elim_.push_back({x, nbs});
        for (int y : nbs) {
          graph.RemoveEdge(x, y);
        }
        assert(graph.Degree(x) == 0);
      }
    }
  }
  SparseGraph ppg(graph.Edges());
  vertex_map2_.resize(ppg.n());
  for (int i = 0; i < ppg.n(); i++) {
    vertex_map2_[i] = ppg.MapBack(i);
  }
  return ppg;
}

vector<int> Preprocessor::Reconstruct(vector<int> colors) const {
  if (vertex_map2_.size() > 0) {
    assert(colors.size() == vertex_map2_.size());
    int nn = org_graph2.n();
    vector<int> col0(nn);
    for (int i=0;i<nn;i++){
      col0[i] = -1;
    }
    for (int i = 0; i < (int)colors.size(); i++) {
      assert(vertex_map2_[i] >= 0 && vertex_map2_[i] < nn);
      col0[vertex_map2_[i]] = colors[i];
    }
    for (int i = (int)tamaki_elim_.size() - 1; i >= 0; i--) {
      int x = tamaki_elim_[i].F;
      auto nbs = tamaki_elim_[i].S;
      assert(x >= 0 && x < nn && col0[x] == -1);
      int mic = nn;
      for (int nx : nbs) {
        mic = std::min(mic, col0[nx]);
      }
      assert(mic > 0);
      col0[x] = mic-1;
    }
    for (int i=0;i<nn;i++){
      assert(col0[i] >= 0);
    }
    colors = col0;
  }

  vector<int> col(n_);
  for (int i = 0; i < n_; i++) {
    col[i] = -1;
  }
  assert(colors.size() == vertex_map_.size());
  int td = 0;
  for (int i = 0; i < (int)colors.size(); i++) {
    col[vertex_map_[i]] = colors[i];
    td = max(td, colors[i]+1);
  }
  for (const auto& tree : trees_) {
    int depth = tree.clqs.back().size();
    int attach = tree.clqs[0][0];
    int maxc = col[attach];
    for (int v : tree.dummys) {
      maxc = max(maxc, col[v]);
    }
    assert(maxc >= depth-1);
    if (maxc >= depth) {
      for (int v : tree.vertices) {
        assert(tree_lab_[v] < maxc);
        col[v] = tree_lab_[v];
      }
      col[attach] = maxc;
    } else {
      vector<int> thr(depth);
      int mr = -1;
      for (const auto& clq : tree.clqs) {
        if (col[clq[0]] > mr) {
          mr = col[clq[0]];
          thr[mr] = 1;
        }
        for (int v : clq) {
          if (col[v] > mr) {
            thr[col[v]] = 1;
          }
        }
      }
      int fo = -1;
      for (int i = depth-1;i>=0;i--){
        if (tree.rvs[i].size() > 0) {
          assert(thr[i]>0);
        } else if (thr[i]>0) {
          fo = i;
          break;
        }
      }
      for (int v : tree.vertices) {
        col[v] = tree_lab_[v];
      }
      if (fo < tree_lab_[attach]) {
        col[attach] = tree_lab_[attach];
      } else {
        col[attach] = fo;
      }
    }
  }
  for (int i=0;i<n_;i++){
    assert(col[i] >= 0 && col[i] < td);
  }
  return col;
}

std::vector<int> ColToPar(const SparseGraph& graph, const std::vector<int>& col) {
  assert((int)col.size() == graph.n());
  int td = 0;
  for (int i = 0; i < graph.n(); i++) {
    assert(col[i] >= 0);
    td = max(td, col[i]+1);
  }
  vector<vector<int>> vs(td);
  vector<int> un(graph.n());
  vector<int> par(graph.n());
  for (int i = 0; i < graph.n(); i++) {
    vs[col[i]].push_back(i);
    un[i] = i;
    par[i] = -1;
  }
  for (int i = 0; i < td; i++) {
    for (int v : vs[i]) {
      assert(col[v] == i);
      for (int nv : graph.Neighbors(v)) {
        if (col[nv] > i) continue;
        assert(col[nv] < i);
        int uu = utils::GetU(nv, un);
        if (uu == v) continue;
        assert(col[uu] < i);
        assert(par[uu] == -1);
        par[uu] = v;
        un[uu] = v;
      }
    }
  }
  int cnt = 0;
  for (int i = 0; i < graph.n(); i++) {
    if (par[i] == -1) {
      cnt++;
    } else {
      assert(par[i] >= 0 && par[i] < graph.n());
    }
  }
  assert(cnt == 1);
  return par;
}
} // namespace sms