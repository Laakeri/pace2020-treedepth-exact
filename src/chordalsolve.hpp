#pragma once

#include <map>

#include "graph.hpp"
#include "bitset.hpp"

#define F first
#define S second

namespace sms {
using std::max;
using std::vector;

template<size_t chunks>
struct MinsepInfo {
  FBitset<chunks> minsep;
  int comp_size, sep_size;
  bool operator==(const MinsepInfo& other) const {
    return minsep == other.minsep;
  }
  int ExpSize() const {
    return comp_size;
  }
};

template<size_t chunks>
struct PieceC {
  PieceC () {
    lb = -1;
    ub = (int)1e9;
  }
  int lb, ub;
  std::vector<MinsepInfo<chunks>> minseps;
};

template<size_t chunks>
class ChordalSolve {
 public:
  ChordalSolve(const FGraph<chunks>& graph);
  int Solve(int goal, int tune, double timelimit);
  std::vector<int> Get(int goal);
 private:
  int PieceId(const FBitset<chunks>& piece, bool insert, bool expect);
  int glo_goal_;
  int tune_;
  double time_limit_;
  Timer timer_;
  FGraph<chunks> graph_;
  std::vector<int> resu_;
  std::vector<PieceC<chunks>> pcs_;
  FBitsetMap<chunks> bs_cac_;
  FLBSieve<chunks> lb_sieve_;
  bool Go(FBitset<chunks> vert, int k, const std::vector<Edge>& parent_edges, const std::vector<MinsepInfo<chunks>>& parent_minseps, int parent_n);
  bool Reco(FBitset<chunks> vert, int k, const std::vector<Edge>& parent_edges);
};

template<size_t chunks>
ChordalSolve<chunks>::ChordalSolve(const FGraph<chunks>& graph) : graph_(graph) {
  assert(!bs_cac_.Inited());
  bs_cac_ = FBitsetMap<chunks>(2, 2);
  assert(bs_cac_.Inited());
  pcs_.push_back({});
  lb_sieve_ = FLBSieve<chunks>(graph.n(), graph.n()+1);
}

template<size_t chunks>
int ChordalSolve<chunks>::PieceId(const FBitset<chunks>& piece, bool insert, bool expect) {
  assert(pcs_.size() >= 1);
  if (insert) {
    int id = bs_cac_.Insert(piece, (int)pcs_.size()).first;
    assert(id > 0);
    if (id == (int)pcs_.size()) {
      pcs_.push_back({});
      assert(!expect);
      return id;
    } else {
      assert(id < (int)pcs_.size());
      return id;
    }
  } else {
    int id = bs_cac_.Get(piece, expect);
    assert(id >= 0);
    if (expect) {
      assert(id > 0);
    }
    assert(id < (int)pcs_.size());
    return id;
  }
}

template<size_t chunks>
bool ChordalSolve<chunks>::Reco(FBitset<chunks> vert, int k, const vector<Edge>& parent_edges) {
  int pc_id = bs_cac_.Get(vert, true);
  assert(pcs_[pc_id].ub <= k);
  int n = vert.Popcount();
  if (n <= k) {
    for (int v : vert) {
      assert(k >= 1);
      resu_[v] = k--;
    }
    return true;
  }
  FGraph<chunks> t_graph(graph_.n());
  std::vector<Edge> tg_edges;
  if (parent_edges.empty()) {
    for (auto e : graph_.Edges()) {
      if (vert.Get(e.F) && vert.Get(e.S)) {
        t_graph.AddEdge(e);
        tg_edges.push_back(e);
      }
    }
  } else {
    for (auto e : parent_edges) {
      if (vert.Get(e.F) && vert.Get(e.S)) {
        t_graph.AddEdge(e);
        tg_edges.push_back(e);
      }
    }
  }

  assert(t_graph.m() < n*(n-1)/2);
  if (k == n-1) {
    for (int u : vert) {
      for (int v : vert) {
        if (u != v && !t_graph.HasEdge(u, v)) {
          resu_[u] = 1;
          resu_[v] = 1;
          vert.SetFalse(u);
          vert.SetFalse(v);
          for (int w : vert) {
            assert(k >= 2);
            resu_[w] = k--;
          }
          return true;
        }
      }
    }
  }
  assert(k <= n-2);
  assert(n >= 4 && t_graph.m() >= 3);
  for (const MinsepInfo<chunks>& ms : pcs_[pc_id].minseps) {
    if (ms.sep_size >= k) continue;
    vert.TurnOff(ms.minsep);
    auto bcs = t_graph.BitComps(vert);
    vert |= ms.minsep;
    assert(bcs.size() >= 2);
    auto cmp = [](const FBitset<chunks>& a, const FBitset<chunks>& b) {
      return a.Popcount() > b.Popcount();
    };
    std::sort(bcs.begin(), bcs.end(), cmp);
    bool ok = true;
    for (const FBitset<chunks>& bc : bcs) {
      size_t bcid = bs_cac_.Get(bc, false);
      if (bcid == 0) {
        ok = false;
        break;
      }
      if (pcs_[bcid].ub > k - ms.sep_size) {
        ok = false;
        break;
      }
    }
    if (ok) {
      for (const FBitset<chunks>& bc : bcs) {
        assert(Reco(bc, k - ms.sep_size, tg_edges));
      }
      for (int v : ms.minsep) {
        assert(k >= 1);
        resu_[v] = k--;
      }
      return true;
    }
  }
  assert(false);
}

Timer cmst1,cmst2,cmst3,csubt,ctott,ctg;

template<size_t chunks>
bool ChordalSolve<chunks>::Go(FBitset<chunks> vert, int k, const vector<Edge>& parent_edges, const vector<MinsepInfo<chunks>>& parent_minseps, int parent_n) {
  ctott.start();
  int pc_id = PieceId(vert, true, false);
  int n = vert.Popcount();
  if (n == 1) {
    pcs_[pc_id].lb = 1;
  } else {
    pcs_[pc_id].lb = max(pcs_[pc_id].lb, 2);
  }
  if (pcs_[pc_id].lb > k) {
    return false;
  }
  pcs_[pc_id].ub = min(pcs_[pc_id].ub, n);
  if (pcs_[pc_id].ub <= k) {
    return true;
  }
  csubt.start();
  int ss_lb = lb_sieve_.Get(vert, k+1);
  if (ss_lb > 0) {
    assert(ss_lb > k);
    pcs_[pc_id].lb = ss_lb;
    csubt.stop();
    return false;
  }
  csubt.stop();
  FGraph<chunks> t_graph(graph_.n());
  std::vector<Edge> tg_edges;
  for (auto e : parent_edges) {
    if (vert.Get(e.F) && vert.Get(e.S)) {
      tg_edges.push_back(e);
      t_graph.AddEdge(e);
    }
  }
  if ((int64_t)t_graph.m() == (int64_t)n*(n-1)/2) {
    pcs_[pc_id].ub = n;
    pcs_[pc_id].lb = n;
    csubt.start();
    lb_sieve_.Insert(vert, pcs_[pc_id].lb);
    csubt.stop();
    return false;
  } else {
    assert((int64_t)t_graph.m() < (int64_t)n*(n-1)/2);
    pcs_[pc_id].ub = min(pcs_[pc_id].ub, n - 1);
    if (pcs_[pc_id].ub <= k) {
      return true;
    }
  }
  int tw = mcs::Treewidth(t_graph);
  pcs_[pc_id].lb = max(pcs_[pc_id].lb, tw + 1);
  if (pcs_[pc_id].lb > k) {
    csubt.start();
    lb_sieve_.Insert(vert, pcs_[pc_id].lb);
    csubt.stop();
    return false;
  }
  assert(n >= 4 && t_graph.m() >= 3);
  int glo_k = glo_goal_ - graph_.Neighbors(vert).Popcount();
  assert(glo_k >= k);
  if (pcs_[pc_id].minseps.empty()) {
    if (glo_k == k && n > parent_n/2 && !parent_minseps.empty()) {
      cmst1.start();
      for (const auto& pms : parent_minseps) {
        int sep_size = vert.IntersectionPopcount(pms.minsep);
        if (sep_size == 0) continue;
        assert(sep_size < glo_k);
        FBitset<chunks> vert2 = vert;
        vert2.TurnOff(pms.minsep);
        FBitset<chunks> vert4 = pms.minsep;
        vert4 &= vert;
        int fcs = 0;
        int comp_size = 0;
        int sz = n - sep_size;
        for (int i=0;i<t_graph.n();i++) {
          if (vert2.Get(i)) {
            FBitset<chunks> vert3 = t_graph.adj_mat2_[i];
            t_graph.Dfs2Bit(vert2, vert3);
            int nsz = vert2.Popcount();
            comp_size = max(comp_size, sz-nsz);
            sz = nsz;
            if (vert3.Subsumes(vert4)) {
              fcs++;
            }
          }
        }
        if (fcs < 2) continue;
        pcs_[pc_id].minseps.push_back({vert4, comp_size, sep_size});
      }
      cmst1.stop();
    } else {
      cmst2.start();
      for (const FBitset<chunks>& ms : mcs::ChordalMinseps(t_graph)) {
        int sep_size = ms.Popcount();
        assert(sep_size < k && sep_size >= 1);
        FBitset<chunks> vert2 = vert;
        vert2.TurnOff(ms);
        int fcs = 0;
        int comp_size = 0;
        int sz = n - sep_size;
        for (int i=0;i<t_graph.n();i++) {
          if (vert2.Get(i)) {
            FBitset<chunks> vert3 = t_graph.adj_mat2_[i];
            t_graph.Dfs2Bit(vert2, vert3);
            int nsz = vert2.Popcount();
            comp_size = max(comp_size, sz-nsz);
            sz = nsz;
            if (vert3.Subsumes(ms)) {
              fcs++;
            }
          }
        }
        assert(fcs >= 2);
        pcs_[pc_id].minseps.push_back({ms, comp_size, sep_size});
      }
      cmst2.stop();
    }
    cmst3.start();
    auto cmp = [&](const MinsepInfo<chunks>& a, const MinsepInfo<chunks>& b) {
      if (a.comp_size == b.comp_size) {
        if (a.sep_size == b.sep_size) {
          return a.minsep < b.minsep;
        }
        return a.sep_size < b.sep_size;
      }
      return a.comp_size < b.comp_size;
    };
    std::sort(pcs_[pc_id].minseps.begin(), pcs_[pc_id].minseps.end(), cmp);
    pcs_[pc_id].minseps.erase(std::unique(pcs_[pc_id].minseps.begin(), pcs_[pc_id].minseps.end()), pcs_[pc_id].minseps.end());
    cmst3.stop();
  }
  int branches = 0;
  for (const MinsepInfo<chunks>& ms : pcs_[pc_id].minseps) {
    if (branches >= 1 && timer_.get() > time_limit_) {
      break;
    }
    if (ms.comp_size >= n/2 + k - ms.sep_size + tune_) {
      continue;
    }
    assert(ms.sep_size < k);
    branches++;
    assert(vert.Subsumes(ms.minsep));
    vert.TurnOff(ms.minsep);
    auto bcs = t_graph.BitComps(vert);
    assert(bcs.size() >= 2);
    auto cmp = [](const FBitset<chunks>& a, const FBitset<chunks>& b) {
      return a.Popcount() > b.Popcount();
    };
    std::sort(bcs.begin(), bcs.end(), cmp);
    bool ok = true;
    for (const FBitset<chunks>& bc : bcs) {
      int bcid = PieceId(bc, true, false);
      if (pcs_[bcid].lb > k - ms.sep_size) {
        ok = false;
        break;
      }
    }
    int tans = 0;
    if (ok) {
      for (const FBitset<chunks>& bc : bcs) {
        if (!Go(bc, k - ms.sep_size, tg_edges, pcs_[pc_id].minseps, n)) {
          ok = false;
          break;
        } else {
          tans = max(tans, ms.sep_size + pcs_[PieceId(bc, false, true)].ub);
        }
      }
    }
    vert |= ms.minsep;
    if (ok) {
      assert(tans <= pcs_[pc_id].ub);
      pcs_[pc_id].ub = tans;
      return true;
    }
  }
  pcs_[pc_id].lb = k+1;
  csubt.start();
  lb_sieve_.Insert(vert, pcs_[pc_id].lb);
  csubt.stop();
  return false;
}

template<size_t chunks>
int ChordalSolve<chunks>::Solve(int goal, int tune, double timelimit) {
  ctott.stop();
  csubt.clear();
  cmst1.clear();
  cmst2.clear();
  cmst3.clear();
  ctott.clear();
  ctg.clear();
  time_limit_ = timelimit;
  tune_ = tune;
  timer_.start();
  int tw = mcs::Treewidth(graph_);
  if (tw+1 > goal) return goal+1;
  FBitset<chunks> vert;
  vert.FillUpTo(graph_.n());
  glo_goal_ = goal;
  const auto graph_edges = graph_.Edges();
  for (int k = goal; k >= 0; k--) {
    glo_goal_ = k;
    if (!Go(vert, k, graph_edges, {}, graph_.n())) {
      if (timer_.get() > time_limit_ || tune_ >= graph_.n()/2) {
        Log::Write(10, "fin ", k+1, " ", tune_);
        Log::Write(10, "times ", csubt.get(), " ", cmst1.get(), " ", cmst2.get(), " ", cmst3.get(), " ", ctott.get(), " ", ctg.get());
        return k+1;
      }
      tune_++;
      lb_sieve_ = FLBSieve<chunks>(graph_.n(), graph_.n()+1);
      for (int i = 1; i < (int)pcs_.size(); i++) {
        pcs_[i].lb = 0;
      }
      k++;
    }
  }
  return 0;
}

template<size_t chunks>
vector<int> ChordalSolve<chunks>::Get(int goal) {
  resu_.resize(graph_.n());
  FBitset<chunks> vert;
  vert.FillUpTo(graph_.n());
  Timer recot;
  recot.start();
  assert(Reco(vert, goal, {}));
  Log::Write(5, "recot ", recot.get());
  for (int i=0;i<graph_.n();i++) {
    assert(resu_[i] >= 1 && resu_[i] <= goal);
    resu_[i]--;
  }
  return resu_;
}
} // namespace sms