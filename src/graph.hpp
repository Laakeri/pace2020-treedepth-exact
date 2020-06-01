#pragma once

#include <vector>
#include <ostream>
#include <set>
#include <queue>

#include "utils.hpp"
#include "staticset.hpp"
#include "bitset.hpp"

namespace triangulator {
namespace {
using std::vector;
using std::queue;
using std::min;

template<size_t chunks>
struct MaxFlow2 {
  vector<FBitset<chunks>> g;
  FBitset<chunks> u;
  MaxFlow2(int n) : g(n) {
    u.FillTrue();
  }
  void addEdge(int a, int b) {
    g[a].SetTrue(b);
  }
  int flow(int x, int sink) {
    if (x==sink) return 1;
    assert(u.Get(x));
    u.SetFalse(x);
    for (size_t ch = 0; ch < chunks; ch++) {
      while (u.data_[ch] & g[x].data_[ch]) {
        int y = __builtin_ctzll(u.data_[ch] & g[x].data_[ch]) + ch*BITS;
        if (flow(y, sink)) {
          g[x].SetFalse(y);
          g[y].SetTrue(x);
          return 1;
        }
      }
    }
    return 0;
  }
  int getMaxFlow(int source, int sink, int ma) {
    int r=0;
    while (int t=flow(source, sink)) {
      r += t;
      if (r > ma) return r;
      u.FillTrue();
    }
    return r;
  }
};

template<size_t chunks, size_t chunks2>
FBitset<chunks2> getMinCut(MaxFlow2<chunks>& mf, int source, int sink, int n) {
  mf.u.FillTrue();
  queue<int> q;
  q.push(source);
  mf.u.SetFalse(source);
  while (!q.empty()) {
    int x = q.front();
    q.pop();
    assert(x != sink);
    for (size_t ch = 0; ch < chunks; ch++) {
      while (mf.u.data_[ch] & mf.g[x].data_[ch]) {
        int y = __builtin_ctzll(mf.u.data_[ch] & mf.g[x].data_[ch]) + ch*BITS;
        q.push(y);
        mf.u.SetFalse(y);
      }
    }
  }
  FBitset<chunks2> ret;
  for (int i = 0; i < n; i ++) {
    if (!mf.u.Get(2+i*2) && mf.u.Get(3+i*2)) {
      ret.SetTrue(i);
    }
  }
  return ret;
}
} // namespace

typedef std::pair<int, int> Edge;

class SparseGraph;

template<size_t chunks>
class FGraph;

// SPARSEGRAPH
class SparseGraph {
 public:
  explicit SparseGraph(int n);
  explicit SparseGraph(std::vector<Edge> edges);
  template<size_t chunks>
  explicit SparseGraph(const FGraph<chunks>& graph);

  int n() const;
  int m() const;
  bool HasEdge(int v, int u) const;
  bool HasEdge(Edge e) const;

  void AddEdge(int v, int u);
  void AddEdge(Edge e);
  void AddEdges(const std::vector<Edge>& edges);

  void RemoveEdge(int v, int u);

  const std::vector<int>& Neighbors(int v) const;

  int Degree(int v) const;
  std::vector<int> Distances(const std::vector<int>& start) const;

  void Print(std::ostream& out) const;

  std::vector<Edge> Edges() const;

  int MapBack(int v) const;

  bool IsConnected() const;
  std::vector<std::vector<int> > Components(const std::vector<int>& separator) const;
  std::vector<int> FindComponentAndMark(int v, std::vector<char>& block) const;

  StaticSet<int> VertexMap() const;

  bool IsClique(const std::vector<int>& vs) const;
  void ShuffleAdjList(std::mt19937& gen);

  int Mincut(int a, int b) const;

 private:
  int n_,m_;
  StaticSet<int> vertex_map_;
  std::vector<std::vector<int> > adj_list_;
  void Dfs(int v, std::vector<char>& blocked, std::vector<int>& component) const;
};



// FIXEDGRAPH
template <size_t chunks>
class FGraph {
 public:
  FGraph(int n);
  FGraph(std::vector<Edge> edges);
  FGraph(const SparseGraph& graph);
  int n() const;
  int m() const;
  bool HasEdge(int v, int u) const;
  bool HasEdge(Edge e) const;
  void AddEdge(int v, int u);
  void AddEdge(Edge e);
  const std::vector<int>& Neighbors(int v) const;
  std::vector<FBitset<chunks>> CompNeighsBit(const FBitset<chunks>& block) const;
  void Dfs2Bit(FBitset<chunks>& vis, FBitset<chunks>& ne) const;
  std::vector<FBitset<chunks>> SmallMinsepsHeuristic(int sz) const;
  void Dfs(int v, std::vector<char>& block, std::vector<int>& component) const;
  std::vector<int> FindComponentAndMark(int v, std::vector<char>& block) const;
  bool IsConnectedOrIsolated() const;
  std::vector<std::vector<int> > Components(const std::vector<int>& separator) const;
  uint64_t Hash() const;
  std::vector<Edge> Edges() const;
  StaticSet<int> VertexMap() const;
  std::vector<Edge> FillEdges(FBitset<chunks> bs) const;
  int FillSize(FBitset<chunks> bs) const;
  int Degree(int v) const;
  FBitset<chunks> Neighbors(const FBitset<chunks>& vs) const;
  uint64_t Hash2(const FBitset<chunks>& vert) const;
  std::vector<uint64_t> Labels(const FBitset<chunks>& vert) const;
  std::vector<uint64_t> RefinedLabels(const FBitset<chunks>& vert) const;
  std::vector<FBitset<chunks>> BitComps(FBitset<chunks> vis) const;
  std::vector<FBitset<chunks>> FullComponentsWithSep(const FBitset<chunks>& minsep) const;
  void ShuffleAdjList(std::mt19937& gen);

  int MaxCompSize(const FBitset<chunks>& minsep, const FBitset<chunks>& vert) const;

  bool IsStar(const FBitset<chunks>& vs) const;
  std::vector<FBitset<chunks>> StarMinsep(int sz) const;

  void Print(std::ostream& out) const;

  std::vector<FBitset<chunks>> adj_mat2_;
 private:
  int n_, m_;
  StaticSet<int> vertex_map_;
  std::vector<std::vector<int> > adj_list_;
};






// IMPLEMENTATION
template<size_t chunks>
SparseGraph::SparseGraph(const FGraph<chunks>& graph) {
  n_ = graph.n();
  m_ = graph.m();
  adj_list_.resize(n_);
  for (int i = 0; i < n_; i++) {
    adj_list_[i] = graph.Neighbors(i);
  }
  vertex_map_ = graph.VertexMap();
}


template <size_t chunks>
FGraph<chunks>::FGraph(int n) : n_(n), m_(0), adj_list_(n) {
  assert(BITS * chunks >= n_);
  adj_mat2_.resize(n_);
  std::vector<int> identity(n);
  for (int i = 0; i < n; i++) {
    identity[i] = i;
    adj_mat2_[i].SetTrue(i);
  }
  vertex_map_.Init(identity);
}

template <size_t chunks>
FGraph<chunks>::FGraph(std::vector<Edge> edges) : vertex_map_(edges) {
  n_ = vertex_map_.Size();
  assert(BITS * chunks >= n_);
  m_ = 0;
  adj_list_.resize(n_);
  adj_mat2_.resize(n_);
  for (int i = 0; i < n_; i++) {
    adj_mat2_[i].SetTrue(i);
  }
  for (auto edge : edges) {
    AddEdge(vertex_map_.Rank(edge.first), vertex_map_.Rank(edge.second));
  }
}

template <size_t chunks>
FGraph<chunks>::FGraph(const SparseGraph& graph) {
  n_ = graph.n();
  assert(BITS * chunks >= n_);
  m_ = 0;
  adj_list_.resize(n_);
  adj_mat2_.resize(n_);
  for (int i = 0; i < n_; i++) {
    adj_mat2_[i].SetTrue(i);
  }
  for (auto edge : graph.Edges()) {
    assert(0 <= edge.first && edge.first < edge.second && edge.second < n_);
    AddEdge(edge.first, edge.second);
  }
  vertex_map_ = graph.VertexMap();
}

template <size_t chunks>
int FGraph<chunks>::n() const {
  return n_;
}

template <size_t chunks>
int FGraph<chunks>::m() const {
  return m_;
}

template <size_t chunks>
bool FGraph<chunks>::HasEdge(int v, int u) const {
  return adj_mat2_[v].Get(u);
}

template <size_t chunks>
bool FGraph<chunks>::HasEdge(Edge e) const {
  return HasEdge(e.first, e.second);
}

template <size_t chunks>
void FGraph<chunks>::AddEdge(int v, int u) {
  if (HasEdge(v, u)) return;
  assert(v != u);
  m_++;
  adj_list_[v].push_back(u);
  adj_list_[u].push_back(v);
  adj_mat2_[v].SetTrue(u);
  adj_mat2_[u].SetTrue(v);
}

template <size_t chunks>
void FGraph<chunks>::AddEdge(Edge e) {
  AddEdge(e.first, e.second);
}

template <size_t chunks>
const std::vector<int>& FGraph<chunks>::Neighbors(int v) const {
  return adj_list_[v];
}

template <size_t chunks>
std::vector<FBitset<chunks>> FGraph<chunks>::CompNeighsBit(const FBitset<chunks>& block) const {
  FBitset<chunks> vis = ~block;
  std::vector<FBitset<chunks>> ret;
  FBitset<chunks> ne;
  FBitset<chunks> sep;
  for (int i=0;i<n_;i++){
    if (vis.Get(i)) {
      ne = adj_mat2_[i];
      Dfs2Bit(vis, ne);
      sep.SetAnd(block, ne);
      if (sep.Popcount() > 0) {
        ret.push_back(sep);
      }
    }
  }
  return ret;
}

template <size_t chunks>
void FGraph<chunks>::Dfs2Bit(FBitset<chunks>& vis, FBitset<chunks>& ne) const {
  bool fo = true;
  while (fo) {
    fo = false;
    for (size_t j = 0; j < chunks; j++) {
      uint64_t gv = vis.data_[j] & ne.data_[j];
      vis.data_[j] &= ~gv;
      if (gv) fo = true;
      while (gv) {
        int x = __builtin_ctzll(gv) + j*BITS;
        gv &= ~-gv;
        ne |= adj_mat2_[x];
      }
    }
  }
}

template <size_t chunks>
std::vector<FBitset<chunks>> FGraph<chunks>::SmallMinsepsHeuristic(int sz) const {
  assert(IsConnectedOrIsolated());
  std::vector<FBitset<chunks>> minseps;
  FBitsetSet<chunks> ff(n_, 2);
  for (int i = 0; i < n_; i++) {
    if (Neighbors(i).empty()) continue;
    for (const FBitset<chunks>& nbs : CompNeighsBit(adj_mat2_[i])) {
      if (nbs.Popcount() <= sz && ff.Insert(nbs)) {
        minseps.push_back(nbs);
      }
    }
  }
  FBitset<chunks> vis;
  FBitset<chunks> sep;
  FBitset<chunks> ne;
  FBitset<chunks> mask;
  for (int i=0;i<n_;i++){
    if (!Neighbors(i).empty()) mask.SetTrue(i);
  }
  for (int i = 0;i<(int)minseps.size(); i++) {
    if (i > 0 && i%1000000 == 0) Log::Write(5, "F enum minseps ", i, " ", minseps.size());
    const FBitset<chunks> tsep = minseps[i];
    for (int j : tsep) {
      FBitset<chunks> block = minseps[i];
      block |= adj_mat2_[j];
      vis.SetNegAnd(block, mask);
      for (int ch = 0; ch < chunks; ch++) {
        while (vis.data_[ch] > 0) {
          int k = __builtin_ctzll(vis.data_[ch]) + ch*BITS;
          sep = block;
          ne = adj_mat2_[k];
          Dfs2Bit(vis, ne);
          sep.SetAnd(ne, block);
          if (sep.Popcount() <= sz && ff.Insert(sep)) {
            minseps.push_back(sep);
          }
        }
      }
    }
  }
  for (int i=0;i<(int)minseps.size();i++) {
    assert(minseps[i].Popcount() <= sz);
  }
  return minseps;
}

template <size_t chunks>
void FGraph<chunks>::Dfs(int v, std::vector<char>& block, std::vector<int>& component) const {
  block[v] = true;
  component.push_back(v);
  for (int nv : adj_list_[v]) {
    if (!block[nv]) {
      Dfs(nv, block, component);
    }
  }
}

template <size_t chunks>
std::vector<int> FGraph<chunks>::FindComponentAndMark(int v, std::vector<char>& block) const {
  std::vector<int> component;
  Dfs(v, block, component);
  return component;
}

template <size_t chunks>
bool FGraph<chunks>::IsConnectedOrIsolated() const {
  auto cs = Components({});
  int f = 0;
  for (const auto& c : cs) {
    if ((int)c.size() > 1) f++;
  }
  return f <= 1;
}

template <size_t chunks>
std::vector<std::vector<int> > FGraph<chunks>::Components(const std::vector<int>& separator) const {
  std::vector<char> blocked(n_);
  for (int v : separator) {
    blocked[v] = true;
  }
  std::vector<std::vector<int> > components;
  for (int i = 0; i < n_; i++) {
    if (!blocked[i]) {
      components.push_back(FindComponentAndMark(i, blocked));
    }
  }
  return components;
}

template <size_t chunks>
uint64_t FGraph<chunks>::Hash() const {
  PolyHash p;
  for (int i = 0; i < n_; i++) {
    for (int a : adj_mat2_[i]) {
      if (a > i) {
        p.Add(i);
        p.Add(a);
      }
    }
  }
  return p.Value();
}

template <size_t chunks>
std::vector<Edge> FGraph<chunks>::Edges() const {
  std::vector<Edge> ret;
  for (int i = 0; i < n_; i++) {
    for (int a : adj_list_[i]) {
      if (a > i) ret.push_back({i, a});
    }
  }
  return ret;
}

template <size_t chunks>
StaticSet<int> FGraph<chunks>::VertexMap() const {
  return vertex_map_;
}

template <size_t chunks>
std::vector<Edge> FGraph<chunks>::FillEdges(FBitset<chunks> bs) const {
  std::vector<Edge> ret;
  for (int i=0;i<chunks;i++){
    while (bs.data_[i]) {
      int v = i*BITS + __builtin_ctzll(bs.data_[i]);
      bs.data_[i] &= ~-bs.data_[i];
      for (int j=i;j<chunks;j++){
        uint64_t td = bs.data_[j] & (~adj_mat2_[v].data_[j]);
        while (td) {
          int u = j*BITS + __builtin_ctzll(td);
          td &= ~-td;
          ret.push_back({v, u});
        }
      }
    }
  }
  return ret;
}

template <size_t chunks>
int FGraph<chunks>::FillSize(FBitset<chunks> bs) const {
  int ans = 0;
  for (int i=0;i<chunks;i++){
    while (bs.data_[i]) {
      int v = i*BITS + __builtin_ctzll(bs.data_[i]);
      bs.data_[i] &= ~-bs.data_[i];
      for (int j=i;j<chunks;j++){
        ans += __builtin_popcountll(bs.data_[j] & (~adj_mat2_[v].data_[j]));
      }
    }
  }
  return ans;
}

template<size_t chunks>
int FGraph<chunks>::Degree(int v) const {
  return adj_list_[v].size();
}

template<size_t chunks>
FBitset<chunks> FGraph<chunks>::Neighbors(const FBitset<chunks>& vs) const {
  FBitset<chunks> nbs;
  for (int v : vs) {
    nbs |= adj_mat2_[v];
  }
  nbs.TurnOff(vs);
  return nbs;
}

template<size_t chunks>
std::vector<uint64_t> FGraph<chunks>::RefinedLabels(const FBitset<chunks>& vert) const {
  std::vector<uint64_t> vh(n_);
  std::vector<FBitset<chunks>> reach(n_);
  int tn = vert.Popcount();
  for (int v : vert) {
    vh[v] = 1;
  }
  int iters = 0;
  while (1) {
    for (int v : vert) {
      reach[v].Clear();
      reach[v].SetTrue(v);
    }
    for (int r = 0; r < tn; r++) {
      std::vector<uint64_t> vh_new(n_);
      bool fo = false;
      for (int v : vert) {
        if (reach[v] != vert) {
          fo = true;
          reach[v] |= Neighbors(reach[v]);
          reach[v] &= vert;
          uint64_t nhv = 1;
          for (int u : adj_list_[v]) {
            if (vert.Get(u)) {
              nhv *= (vh[u] + 1ull);
              nhv %= 1000000007ull;
            }
          }
          PolyHash p;
          p.Add(nhv);
          p.Add(vh[v]);
          vh_new[v] = p.Value();
        } else {
          vh_new[v] = vh[v];
        }
      }
      vh = vh_new;
      if (!fo) break;
    }
    bool bad = false;
    uint64_t min_bad = 0;
    std::set<uint64_t> hl;
    for (int v : vert) {
      if (hl.count(vh[v])) {
        if (!bad) {
          bad = true;
          min_bad = vh[v];
        } else {
          min_bad = std::min(min_bad, vh[v]);
        }
      } else {
        hl.insert(vh[v]);
      }
    }
    if (!bad) {
      return vh;
    } else {
      iters++;
      assert(iters <= 2*tn + 10);
      for (int v : vert) {
        if (vh[v] == min_bad) {
          PolyHash p;
          p.Add(vh[v]);
          p.Add(vh[v]);
          vh[v] = p.Value();
          break;
        }
      }
    }
  }
}

template<size_t chunks>
std::vector<uint64_t> FGraph<chunks>::Labels(const FBitset<chunks>& vert) const {
  std::vector<uint64_t> vh(n_);
  std::vector<FBitset<chunks>> reach(n_);
  int tn = vert.Popcount();
  for (int v : vert) {
    reach[v].Clear();
    reach[v].SetTrue(v);
    vh[v] = 1;
  }
  for (int r = 0; r < tn; r++) {
    std::vector<uint64_t> vh_new(n_);
    bool fo = false;
    for (int v : vert) {
      if (reach[v] != vert) {
        fo = true;
        reach[v] |= Neighbors(reach[v]);
        reach[v] &= vert;
        uint64_t nhv = 1;
        for (int u : adj_list_[v]) {
          if (vert.Get(u)) {
            nhv *= (vh[u] + 1ull);
            nhv %= 1000000007ull;
          }
        }
        nhv += 1ll;
        nhv *= (vh[v] + 2ull);
        nhv %= 1000000007ull;
        vh_new[v] = nhv;
      } else {
        vh_new[v] = vh[v];
      }
    }
    vh = vh_new;
    if (!fo) break;
  }
  return vh;
}

template<size_t chunks>
uint64_t FGraph<chunks>::Hash2(const FBitset<chunks>& vert) const {
  auto vh = Labels(vert);
  std::sort(vh.begin(), vh.end());
  PolyHash p;
  for (uint64_t v : vh) {
    if (v > 0) p.Add(v);
  }
  return p.Value();
}

template<size_t chunks>
std::vector<FBitset<chunks>> FGraph<chunks>::BitComps(FBitset<chunks> vis) const {
  FBitset<chunks> ne;
  std::vector<FBitset<chunks>> ret;
  bool fo = false;
  while (1) {
    if (!fo) {
      for (int j = 0; j < chunks; j++) {
        if (vis.data_[j]) {
          int x = __builtin_ctzll(vis.data_[j]) + j*BITS;
          ne.SetTrue(x);
          ret.push_back(FBitset<chunks>());
          fo = true;
          break;
        }
      }
      if (!fo) return ret;
    }
    fo = false;
    for (int j = 0; j < chunks; j++) {
      uint64_t gv = vis.data_[j] & ne.data_[j];
      while (gv) {
        fo = true;
        vis.data_[j] &= (~(gv&-gv));
        int x = __builtin_ctzll(gv) + j*BITS;
        ne |= adj_mat2_[x];
        ret.back().SetTrue(x);
        gv &= ~-gv;
      }
    }
  }
}

template<size_t chunks>
std::vector<FBitset<chunks>> FGraph<chunks>::FullComponentsWithSep(const FBitset<chunks>& minsep) const {
  FBitset<chunks> vis;
  vis.FillUpTo(n_);
  vis.TurnOff(minsep);
  auto comps = BitComps(vis);
  for (int i = 0; i < (int)comps.size(); i++) {
    if (Neighbors(comps[i]) != minsep) {
      std::swap(comps[i], comps.back());
      comps.pop_back();
    }
  }
  assert(comps.size() >= 2);
  for (auto& c : comps) {
    c |= minsep;
  }
  return comps;
}

template<size_t chunks>
void FGraph<chunks>::ShuffleAdjList(std::mt19937& gen) {
  for (int i = 0; i < n_; i++) {
    std::shuffle(adj_list_[i].begin(), adj_list_[i].end(), gen);
  }
}

template<size_t chunks>
bool FGraph<chunks>::IsStar(const FBitset<chunks>& vs) const {
  bool fb = false;
  for (int v : vs) {
    if (vs.IntersectionPopcount(adj_mat2_[v]) > 2) {
      if (fb) {
        return false;
      } else {
        fb = true;
      }
    }
  }
  return true;
}

template<size_t chunks>
std::vector<FBitset<chunks>> FGraph<chunks>::StarMinsep(int sz) const {
  assert(IsConnectedOrIsolated());
  FBitset<chunks> mask;
  for (int i=0;i<n_;i++){
    if (!Neighbors(i).empty()) mask.SetTrue(i);
  }
  for (int i = 0; i < n_; i++) {
    if (Degree(i) == 0) continue;
    FBitset<chunks> rch = mask;
    for (int v : adj_mat2_[i]) {
      rch.TurnOff(adj_mat2_[v]);
    }
    bool ok1 = true;
    for (const auto& comp : BitComps(rch)) {
      if (!IsStar(comp)) {
        ok1 = false;
        break;
      }
    }
    if (!ok1) continue;
    std::vector<FBitset<chunks>> minseps;
    FBitsetSet<chunks> ff(2, 2);
    for (const FBitset<chunks>& nbs : CompNeighsBit(adj_mat2_[i])) {
      if (ff.Insert(nbs)) {
        minseps.push_back(nbs);
      }
    }
    Timer tt;
    tt.start();
    for (int it = 0; it < (int)minseps.size(); it++) {
      if (it > 0 && it%1000 == 0) Log::Write(5, "star minseps ", it, " ", minseps.size(), " ", tt.get());
      FBitset<chunks> tsep = minseps[it];
      FBitset<chunks> vv = mask;
      vv.TurnOff(tsep);
      bool ok = false;
      for (const auto& comp : BitComps(vv)) {
        if (comp.Get(i) && IsStar(comp) && Neighbors(comp) == tsep) {
          ok = true;
          break;
        }
      }
      if (!ok) continue;
      if (tsep.Popcount() <= sz) {
        assert(ok);
        for (const auto& comp : BitComps(vv)) {
          if (!IsStar(comp)) {
            ok = false;
            break;
          }
        }
        if (ok) {
          return {tsep};
        }
      }
      for (int j : tsep) {
        if (!adj_mat2_[i].Get(j)) continue;
        FBitset<chunks> block = tsep;
        block |= adj_mat2_[j];
        FBitset<chunks> vis;
        vis.SetNegAnd(block, mask);
        for (size_t ch = 0; ch < chunks; ch++) {
          while (vis.data_[ch] > 0) {
            int k = __builtin_ctzll(vis.data_[ch]) + ch*BITS;
            FBitset<chunks> sep = block;
            FBitset<chunks> ne = adj_mat2_[k];
            Dfs2Bit(vis, ne);
            sep.SetAnd(ne, block);
            if (ff.Insert(sep)) {
              minseps.push_back(sep);
            }
          }
        }
      }
    }
  }
  return {};
}

template<size_t chunks>
void FGraph<chunks>::Print(std::ostream& out) const {
  out<<"v e: "<<n_<<" "<<m_<<std::endl;
  for (int i = 0; i < n_; i++) {
    for (int a : adj_list_[i]) {
      out<<i<<" "<<a<<std::endl;
    }
  }
}

template<size_t chunks>
void SepRec(const FGraph<chunks>& graph, int a, int b, FBitset<chunks> neA, FBitset<chunks> neB, FBitset<chunks> F, std::vector<FBitset<chunks>>& minseps, int sz, int n) {
  // Log::Write(3, "seprec ", a, " ", b, " ", A.Elements(), " ", F.Elements(), " ", sz);
  assert(F.Popcount() <= sz);
  FBitset<chunks> inter = neA;
  inter &= neB;
  assert(inter.Subsumes(F));
  if (inter == F) {
    minseps.push_back(F);
    return;
  }
  if (F.Popcount() == sz) return;

  int szthr = std::min(neB.Popcount() - inter.Popcount(), (n-F.Popcount())/2);
  int a_size = neA.Popcount() - inter.Popcount();
  if (a_size > szthr) return;
  if (neA.Popcount() - sz > szthr) return;
  inter.TurnOff(F);
  for (int v : inter) {
    if (graph.HasEdge(b, v)) {
      F.SetTrue(v);
      SepRec(graph, a, b, neA, neB, F, minseps, sz, n);
      return;
    }
  }

  // pruning
  if (a_size + 3*(inter.Popcount() - (sz - F.Popcount())) > szthr) {
    std::vector<std::vector<int>> paths;
    for (int v : inter) {
      paths.push_back({v});
    }
    FBitset<chunks> space = neB;
    space.TurnOff(F);
    space.TurnOff(inter);
    int d = 0;
    int cantake = (int)paths.size() - (sz - F.Popcount());
    for (int len = 1; !paths.empty() && cantake > 0; len++) {
      for (int i = 0; i < (int)paths.size(); i++) {
        assert((int)paths[i].size() == len);
        if (graph.adj_mat2_[paths[i].back()].Intersects(space)) {
          FBitset<chunks> lol = graph.adj_mat2_[paths[i].back()] & space;
          int x = lol.First();
          assert(space.Get(x) && x != paths[i].back() && graph.HasEdge(x, paths[i].back()));
          paths[i].push_back(x);
          space.SetFalse(x);
        } else {
          if (cantake) {
            cantake--;
            d += len;
          }
          std::swap(paths[i], paths.back());
          paths.pop_back();
          i--;
        }
      }
    }
    if (a_size + d > szthr) {
      return;
    }
  }
  int x = inter.First();
  F.SetTrue(x);
  SepRec(graph, a, b, neA, neB, F, minseps, sz, n);
  F.SetFalse(x);

  FBitset<chunks> vis = neB;
  vis.TurnOff(neA);
  vis.TurnOff(graph.adj_mat2_[x]);
  neB = graph.adj_mat2_[b];
  graph.Dfs2Bit(vis, neB);
  if (!neB.Subsumes(F)) return;

  vis.FillTrue();
  vis.TurnOff(neB);
  graph.Dfs2Bit(vis, neA);
  SepRec(graph, a, b, neA, neB, F, minseps, sz, n);
}

template<size_t chunks>
inline std::vector<FBitset<chunks>> NibbleSmallMinseps(FGraph<chunks> graph, int sz) {
  assert(graph.IsConnectedOrIsolated());
  int mfi = graph.n() * graph.n();
  int mfv = graph.n();
  FBitset<chunks> vert;
  for (int i = 0; i < graph.n(); i++) {
    if (graph.Degree(i) > 0) {
      vert.SetTrue(i);
      int fi = graph.FillSize(graph.adj_mat2_[i]);
      if (fi < mfi) {
        mfi = fi;
        mfv = i;
      }
    }
  }
  if (vert.Popcount() <= 2) return {};
  Timer tmr;
  bool wl = false;
  if (vert.Popcount() == graph.n()) {
    Log::Write(5, "nibbleseps ", sz);
    tmr.start();
    wl = true;
  }
  assert(mfv < graph.n() && mfi < graph.n() * graph.n());
  std::vector<FBitset<chunks>> minseps;
  for (int a : graph.Neighbors(mfv)) {
    for (int b : graph.Neighbors(mfv)) {
      if (a == b || graph.HasEdge(a, b)) continue;
      FBitset<chunks> F = graph.adj_mat2_[a];
      F &= graph.adj_mat2_[b];
      if (F.Popcount() <= sz) {
        FBitset<chunks> vis = vert;
        vis.TurnOff(graph.adj_mat2_[a]);
        FBitset<chunks> neB = graph.adj_mat2_[b];
        graph.Dfs2Bit(vis, neB);
        SepRec(graph, a, b, graph.adj_mat2_[a], neB, F, minseps, sz, vert.Popcount());
      }
    }
  }
  vert.TurnOff(graph.adj_mat2_[mfv]);
  for (auto comp : graph.BitComps(vert)) {
    FBitset<chunks> nbs = graph.Neighbors(comp);
    if (nbs.Popcount() <= sz) {
      minseps.push_back(nbs);
    }
    FGraph<chunks> ngraph(graph.n());
    comp |= nbs;
    for (auto e : graph.Edges()) {
      if (comp.Get(e.first) && comp.Get(e.second)) {
        ngraph.AddEdge(e);
      }
    }
    for (auto e : graph.FillEdges(nbs)) {
      ngraph.AddEdge(e);
    }
    auto rms = NibbleSmallMinseps(ngraph, sz);
    for (const auto& sep : rms) {
      minseps.push_back(sep);
    }
  }

  utils::SortAndDedup(minseps);
  if (wl) {
    Log::Write(5, "nibbleseps ret ", minseps.size(), " ", tmr.get());
  }
  return minseps;
}

template<size_t chunks>
int FGraph<chunks>::MaxCompSize(const FBitset<chunks>& minsep, const FBitset<chunks>& vert) const {
  FBitset<chunks> vis = vert;
  vis.TurnOff(minsep);
  int ret = 0;
  for (size_t ch = 0; ch < chunks; ch++) {
    while (vis.data_[ch] > 0) {
      int k = __builtin_ctzll(vis.data_[ch]) + ch*BITS;
      FBitset<chunks> ne = adj_mat2_[k];
      Dfs2Bit(vis, ne);
      for (size_t ch2 = 0; ch2 < chunks; ch2++) {
        ne.data_[ch2] &= ~vis.data_[ch2];
      }
      ret = std::max(ret, ne.Popcount());
      if (ret > vis.Popcount()) {
        return ret;
      }
    }
  }
  return ret;
}
} // namespace triangulator
