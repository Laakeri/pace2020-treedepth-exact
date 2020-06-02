#pragma once

#include <vector>
#include <random>

#include "graph.hpp"

namespace sms {
namespace mcs {
template<size_t chunks>
std::vector<int> Mcs(const FGraph<chunks>& graph);

template<size_t chunks>
void LbTriang(FGraph<chunks>& graph, std::mt19937& gen, int vari, int upd_cnt);

template<size_t chunks>
int Treewidth(const FGraph<chunks>& graph);

template<size_t chunks>
std::vector<FBitset<chunks>> ChordalMinseps(const FGraph<chunks>& graph);

template<size_t chunks>
std::vector<int> Mcs(const FGraph<chunks>& graph) {
  std::vector<int> order(graph.n());
  static std::vector<int> label;
  static std::vector<char> rm;
  static std::vector<std::vector<int> > labels;
  utils::InitZero(label, graph.n());
  utils::InitZero(rm, graph.n());
  if (labels.size() < graph.n()) labels.resize(graph.n());
  for (int i = 0; i < graph.n(); i++) labels[i].clear();
  for (int i = 0; i < graph.n(); i++) labels[0].push_back(i);
  int max_label = 0;
  for (int it = graph.n() - 1; it >= 0; it--) {
    if (labels[max_label].size() == 0) {
      max_label--;
      it++;
      continue;
    }
    int x = labels[max_label].back();
    labels[max_label].pop_back();
    if (rm[x]) {
      it++;
      continue;
    }
    order[it] = x;
    for (int nx : graph.Neighbors(x)) {
      if (!rm[nx]) {
        label[nx]++;
        labels[label[nx]].push_back(nx);
        max_label = std::max(max_label, label[nx]);
      }
    }
    rm[x] = true;
  }
  return order;
}

template<size_t chunks>
int HeurLbt(const FGraph<chunks>& graph, int v) {
  int fec = 0;
  for (const FBitset<chunks>& cn : graph.CompNeighsBit(graph.adj_mat2_[v])) {
    fec += graph.FillSize(cn);
  }
  return fec;
}

template<size_t chunks>
void LbTriang(FGraph<chunks>& graph, std::mt19937& gen, int vari, int upd_cnt) {
  Timer lbt;
  lbt.start();
  std::vector<int> hs(graph.n());
  std::vector<int> upd(graph.n());
  std::vector<int> score(graph.n());
  for (int i=0;i<graph.n();i++) {
    hs[i] = HeurLbt(graph, i)*3 + utils::GetRand(0, vari, gen)*3 + utils::GetRand(0, 2, gen);
  }
  Log::Write(10, "lbt init ", lbt.get());
  lbt.stop();
  for (int it=0;it<graph.n();it++){
    int v = 0;
    int be = 1e9;
    for (int i=0;i<graph.n();i++){
      if (hs[i] >= 0 && hs[i] < be) {
        v = i;
        be = hs[i];
      }
    }
    assert(hs[v] >= 0 && hs[v] == be);
    hs[v] = -1;
    for (const FBitset<chunks>& cn : graph.CompNeighsBit(graph.adj_mat2_[v])) {
      for (auto fe : graph.FillEdges(cn)) {
        graph.AddEdge(fe);
        upd[fe.first] = 1;
        upd[fe.second] = 1;
      }
    }
    if (upd_cnt == 0) continue;
    for (int i=0;i<graph.n();i++){
      if (upd[i]) {
        for (int nb : graph.Neighbors(i)) {
          score[nb]++;
        }
        upd[i] = 0;
        score[i] += 2;
      }
    }
    if (upd_cnt >= graph.n()) {
      for (int i=0;i<graph.n();i++){
        if (hs[i] >= 0 && score[i] > 0) {
          hs[i] = HeurLbt(graph, i)*3 + utils::GetRand(0, vari, gen)*3 + utils::GetRand(0, 2, gen);
          score[i] = 0;
        }
      }
      continue;
    }
    std::vector<std::pair<int, int>> to_upd;
    to_upd.reserve(graph.n());
    for (int i=0;i<graph.n();i++){
      if (hs[i] >= 0 && score[i] > 0) {
        to_upd.push_back({score[i], i});
      }
    }
    std::sort(to_upd.rbegin(), to_upd.rend());
    for (int i=0;i<std::min((int)to_upd.size(), upd_cnt);i++) {
      int u = to_upd[i].second;
      assert(hs[u] >= 0 && score[u] > 0);
      hs[u] = HeurLbt(graph, u)*3 + utils::GetRand(0, vari, gen)*3 + utils::GetRand(0, 2, gen);
      score[u] = 0;
    }
  }
}

template<size_t chunks>
int Treewidth(const FGraph<chunks>& graph) {
  if (graph.m() == 0) return 0;
  std::vector<int> order = Mcs(graph);
  std::vector<int> inv_order = utils::PermInverse(order);
  int treewidth = 0;
  for (int i = 0; i < graph.n(); i++) {
    int x = order[i];
    int nb = 0;
    for (int nx : graph.Neighbors(x)) {
      if (inv_order[nx] > i) {
        nb++;
      }
    }
    treewidth = std::max(treewidth, nb);
  }
  return treewidth;
}

template<size_t chunks>
std::vector<FBitset<chunks>> ChordalMinseps(const FGraph<chunks>& graph) {
  if (graph.m() == 0) return {};
  std::vector<int> order = Mcs(graph);
  std::vector<int> inv_order = utils::PermInverse(order);
  std::vector<FBitset<chunks>> ret;
  int prev = 0;
  for (int i=graph.n()-2;i>=0;i--){
    int cnt = 0;
    for (int nx : graph.Neighbors(order[i])) {
      if (inv_order[nx] > i) {
        cnt++;
      }
    }
    if (cnt <= prev && cnt > 0) {
      FBitset<chunks> ms;
      for (int nx : graph.Neighbors(order[i])) {
        if (inv_order[nx] > i) {
          ms.SetTrue(nx);
        }
      }
      ret.push_back(ms);
    }
    prev = cnt;
  }
  return ret;
}
} // namespace mcs
} // namespace sms
