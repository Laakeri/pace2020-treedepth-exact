#include "graph.hpp"

#include <vector>
#include <algorithm>
#include <set>
#include <cassert>
#include <ostream>
#include <iostream>
#include <queue>

#include "utils.hpp"
#include "bitset.hpp"

namespace triangulator {


SparseGraph::SparseGraph(int n)
  : n_(n), m_(0), adj_list_(n) {
  std::vector<int> identity(n);
  for (int i = 0; i < n; i++) {
    identity[i] = i;
  }
  vertex_map_.Init(identity);
}

SparseGraph::SparseGraph(std::vector<Edge> edges) : vertex_map_(edges) {
  n_ = vertex_map_.Size();
  m_ = 0;
  adj_list_.resize(n_);
  for (auto edge : edges) {
    AddEdge(vertex_map_.Rank(edge.first), vertex_map_.Rank(edge.second));
  }
}

int SparseGraph::n() const {
  return n_;
}

int SparseGraph::m() const {
  return m_;
}

std::vector<Edge> SparseGraph::Edges() const {
  std::vector<Edge> ret;
  for (int i = 0; i < n_; i++) {
    for (int a : adj_list_[i]) {
      if (a > i) ret.push_back({i, a});
    }
  }
  return ret;
}

bool SparseGraph::HasEdge(int v, int u) const {
  if (adj_list_[v].size() > adj_list_[u].size()) {
    return HasEdge(u, v);
  }
  for (int nx : adj_list_[v]) {
    if (nx == u) return true;
  }
  return false;
}

bool SparseGraph::HasEdge(Edge e) const {
  return HasEdge(e.first, e.second);
}

void SparseGraph::AddEdge(int v, int u) {
  if (HasEdge(v, u)) return;
  assert(v != u);
  m_++;
  adj_list_[v].push_back(u);
  adj_list_[u].push_back(v);
}

void SparseGraph::AddEdge(Edge e) {
  AddEdge(e.first, e.second);
}

void SparseGraph::AddEdges(const std::vector<Edge>& edges) {
  for (auto& edge : edges) AddEdge(edge);
}

StaticSet<int> SparseGraph::VertexMap() const {
  return vertex_map_;
}

bool SparseGraph::IsConnected() const {
  auto cs = Components({});
  return (cs.size() == 1) && ((int)cs[0].size() == n_);
}

void SparseGraph::Dfs(int v, std::vector<char>& block, std::vector<int>& component) const {
  block[v] = true;
  component.push_back(v);
  for (int nv : adj_list_[v]) {
    if (!block[nv]) {
      Dfs(nv, block, component);
    }
  }
}

std::vector<int> SparseGraph::FindComponentAndMark(int v, std::vector<char>& block) const {
  std::vector<int> component;
  Dfs(v, block, component);
  return component;
}

std::vector<std::vector<int> > SparseGraph::Components(const std::vector<int>& separator) const {
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

const std::vector<int>& SparseGraph::Neighbors(int v) const {
  return adj_list_[v];
}

int SparseGraph::MapBack(int v) const {
  return vertex_map_.Kth(v);
}

int SparseGraph::Degree(int v) const {
  return adj_list_[v].size();
}

void SparseGraph::RemoveEdge(int v, int u) {
  m_--;
  int fo = 0;
  for (int i = 0; i < adj_list_[v].size(); i++) {
    if (adj_list_[v][i] == u) {
      std::swap(adj_list_[v][i], adj_list_[v].back());
      adj_list_[v].pop_back();
      fo++;
      break;
    }
  }
  for (int i = 0; i < adj_list_[u].size(); i++) {
    if (adj_list_[u][i] == v) {
      std::swap(adj_list_[u][i], adj_list_[u].back());
      adj_list_[u].pop_back();
      fo++;
      break;
    }
  }
  assert(fo == 2);
}

void SparseGraph::Print(std::ostream& out) const {
  out<<"v e: "<<n_<<" "<<m_<<std::endl;
  for (int i = 0; i < n_; i++) {
    for (int a : adj_list_[i]) {
      out<<i<<" "<<a<<std::endl;
    }
  }
}

bool SparseGraph::IsClique(const std::vector<int>& vs) const {
  for (int i = 0; i < (int)vs.size(); i++) {
    for (int ii = i+1; ii < (int)vs.size(); ii++) {
      if (!HasEdge(vs[i], vs[ii])) {
        return false;
      }
    }
  }
  return true;
}

void SparseGraph::ShuffleAdjList(std::mt19937& gen) {
  for (int i = 0; i < n_; i++) {
    std::shuffle(adj_list_[i].begin(), adj_list_[i].end(), gen);
  }
}

namespace {
struct MaxFlow3 {
  std::vector<std::vector<char>> g;
  std::vector<char> u;
  MaxFlow3(int n) : g(n), u(n) {
    for (int i = 0; i < n; i++) {
      u[i] = true;
      g[i].resize(n);
    }
  }
  void addEdge(int a, int b) {
    g[a][b] = 1;
  }
  int flow(int x, int sink) {
    if (x==sink) return 1;
    assert(u[x]);
    u[x] = false;
    for (int y = 0; y < (int)g[x].size(); y++) {
      while (u[y] && g[x][y]) {
        if (flow(y, sink)) {
          g[x][y] = false;
          g[y][x] = true;
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
      for (int i = 0; i < (int)u.size(); i++) {
        u[i] = true;
      }
    }
    return r;
  }
};
} // namespace

int SparseGraph::Mincut(int a, int b) const {
  assert(a >= 0 && b >= 0 && a < n_ && b < n_ && a != b && !HasEdge(a, b));
  MaxFlow3 mf(n_ * 2 + 2);
  for (int i = 0; i < n_; i++) {
    mf.addEdge(i*2+2, i*2+3);
    for (int nx : adj_list_[i]) {
      mf.addEdge(i*2+3, nx*2+2);
    }
  }
  for (int na : adj_list_[a]) {
    mf.addEdge(0, na*2+2);
  }
  for (int nb : adj_list_[b]) {
    mf.addEdge(nb*2+3, 1);
  }
  int fl = mf.getMaxFlow(0, 1, n_);
  return fl;
}
} // namespace triangulator                                           
