#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
#include <set>
#include <cassert>
#include <random>
#include <sys/resource.h>

#include "graph.hpp"
#include "io.hpp"
#include "utils.hpp"
#include "mcs.hpp"
#include "staticset.hpp"
#include "bitset.hpp"
#include "chordalsolve.hpp"
#include "best.hpp"
#include "preprocessor.hpp"
#include "ms_solve.hpp"

using namespace triangulator;

#define F first
#define S second

using std::vector;

std::mt19937 gen(1337);

void SetStackSize(int64_t sz) {
  struct rlimit rl;
  assert(getrlimit(RLIMIT_STACK, &rl) == 0);
  Log::Write(3, "Cur stack size ", rl.rlim_cur);
  if (rl.rlim_cur < sz) {
    rl.rlim_cur = sz;
    Log::Write(3, "Setting stack size ", sz);
    assert(setrlimit(RLIMIT_STACK, &rl) == 0);
  }
}

template<size_t chunks>
int HeurComp(const FGraph<chunks>& graph, int best, double time, const Preprocessor& pp) {
  Timer timer;
  timer.start();
  int it=0;
  std::set<uint64_t> gs;
  int vari = 0;
  int upd_cnt = 0;
  int last_add = 0;
  while (timer.get() < time) {
    double dupls = 0;
    if (it > 0) {
      dupls = (double)(it - (int)gs.size()) / (double)it;
    }
    if (dupls > 0.5 && upd_cnt == graph.n() && it - last_add > 10) {
      vari++;
      last_add = it;
    }
    it++;
    Timer triang_tmr;
    triang_tmr.start();
    FGraph<chunks> lol_g = graph;
    mcs::LbTriang(lol_g, gen, vari, upd_cnt);
    triang_tmr.stop();
    double est = (double)gs.size() * (double)it / ((double)it - (double)gs.size());
    Log::Write(10, "min tri ", triang_tmr.get(), " ", best, " ", est, " ", lol_g.m(), " ", dupls, " ", vari, " ", upd_cnt);
    upd_cnt = upd_cnt * 2 + 1;
    upd_cnt = std::min(upd_cnt, graph.n());
    if (gs.count(lol_g.Hash())) {
      Log::Write(10, "Same triang ", gs.size(), " ", it, " ", est);
      continue;
    }
    gs.insert(lol_g.Hash());
    {
      Timer td_tmr;
      td_tmr.start();
      ChordalSolve<chunks> cs(lol_g);
      int td = cs.Solve(best-1, vari, std::min(time - timer.get(), triang_tmr.get() + 0.01));
      if (td < best) {
        best = td;
        Log::Write(3, "Treedepth: ", best);
        auto resu = cs.Get(best);
        resu = pp.Reconstruct(resu);
        resu = ColToPar(pp.org_graph, resu);
        int got = best::SetBest(resu, true);
        assert(got <= td);
        best = got;
        Log::Write(3, "Got ", got);
      }
    }
  }
  return best;
}

template<size_t chunks> int DoSolve2(const SparseGraph& graph, int best, const Preprocessor& pp) {
  assert(graph.n() <= chunks * BITS && graph.n() > (chunks-1) * BITS);
  FGraph<chunks> ppg(graph);
  Log::Write(3, "Solve2 n:", ppg.n(), " m:", ppg.m());
  {
    MSSolve<chunks> mss(ppg);
    mss.incorrect_msenum_ = true;
    int ans = mss.Solve(best-1, true);
    if (ans < best) {
      best = ans;
      Log::Write(3, "Heur ans ", ans);
      auto sol = mss.Get(ans);
      sol = pp.Reconstruct(sol);
      sol = ColToPar(pp.org_graph, sol);
      int got = best::SetBest(sol, true);
      assert(got <= ans);
      Log::Write(3, "Ans valid ", got, " ", ans);
      best = got;
      Log::Write(3, "Re preprocess");
      return best;
    }
  }
  MSSolve<chunks> mss2(ppg);
  int ans2 = mss2.Solve(best-1, false);
  if (ans2 < best) {
    best = ans2;
    Log::Write(3, "Exact ans ", ans2);
    auto sol = mss2.Get(ans2);
    sol = pp.Reconstruct(sol);
    sol = ColToPar(pp.org_graph, sol);
    int got = best::SetBest(sol, true);
    assert(got == ans2);
    Log::Write(3, "Ans valid ", ans2);
  }
  return -1;
}

template<size_t chunks> void DoSolve1(const SparseGraph& graph, int best, const Preprocessor& pp) {
  assert(graph.n() <= chunks * BITS && graph.n() > (chunks-1) * BITS);
  const FGraph<chunks> ppg(graph);
  Log::Write(3, "Dosolve1 n:", ppg.n(), " m:", ppg.m());

  // This is variable mostly to reduce the total time
  double pp_time = 40;
  if (ppg.n() <= 50) {
    pp_time = 1;
  } else if (ppg.n() <= 75) {
    pp_time = 5;
  } else if (ppg.n() <= 100) {
    pp_time = 20;
  } else if (ppg.n() <= 150) {
    pp_time = 30;
  } else if (ppg.n() <= 200) {
    pp_time = 40;
  } else if (ppg.n() <= 250) {
    pp_time = 50;
  } else {
    pp_time = 60;
  }

  best = HeurComp<chunks>(ppg, best, pp_time, pp);

  while (true) {
    Preprocessor pp2 = pp;
    SparseGraph pp_graph = pp2.TamakiRules(SparseGraph(ppg), best-1);
    int nbest = best;
    if (pp_graph.n() <= BITS) {
      nbest = DoSolve2<1>(pp_graph, best, pp2);
    } else if (pp_graph.n() <= 2*BITS) {
      nbest = DoSolve2<2>(pp_graph, best, pp2);
    } else if (pp_graph.n() <= 3*BITS) {
      nbest = DoSolve2<3>(pp_graph, best, pp2);
    } else if (pp_graph.n() <= 4*BITS) {
      nbest = DoSolve2<4>(pp_graph, best, pp2);
    } else if (pp_graph.n() <= 5*BITS) {
      nbest = DoSolve2<5>(pp_graph, best, pp2);
    } else if (pp_graph.n() <= 6*BITS) {
      nbest = DoSolve2<6>(pp_graph, best, pp2);
    } else if (pp_graph.n() <= 7*BITS) {
      nbest = DoSolve2<7>(pp_graph, best, pp2);
    } else if (pp_graph.n() <= 8*BITS) {
      nbest = DoSolve2<8>(pp_graph, best, pp2);
    } else if (pp_graph.n() <= 9*BITS) {
      nbest = DoSolve2<9>(pp_graph, best, pp2);
    } else if (pp_graph.n() <= 10*BITS) {
      nbest = DoSolve2<10>(pp_graph, best, pp2);
    } else { // Assume that the graph has n <= 640
      assert(0);
    }
    if (nbest == -1) return;
    assert(nbest >= 0 && nbest < best);
    best = nbest;
    Log::Write(3, "Re solve ", best);
  }
}

int main() {
  // Set stack size equal to the memory limit (8GB).
  SetStackSize(8ll * 1024 * 1024);
  // Should be 3 in the final submission. 5 is reasonable, 10 prints a lot.
  Log::SetLogLevel(3);
  Io io;
  SparseGraph graph = io.ReadGraph(std::cin);

  Log::Write(3, "Input n:", graph.n(), " m:", graph.m());
  best::InitBest(graph);
  assert(graph.IsConnected());
  int best = graph.n();

  Preprocessor pp;
  SparseGraph pp_graph = pp.Preprocess(graph);

  if (pp_graph.n() <= BITS) {
    DoSolve1<1>(pp_graph, best, pp);
  } else if (pp_graph.n() <= 2*BITS) {
    DoSolve1<2>(pp_graph, best, pp);
  } else if (pp_graph.n() <= 3*BITS) {
    DoSolve1<3>(pp_graph, best, pp);
  } else if (pp_graph.n() <= 4*BITS) {
    DoSolve1<4>(pp_graph, best, pp);
  } else if (pp_graph.n() <= 5*BITS) {
    DoSolve1<5>(pp_graph, best, pp);
  } else if (pp_graph.n() <= 6*BITS) {
    DoSolve1<6>(pp_graph, best, pp);
  } else if (pp_graph.n() <= 7*BITS) {
    DoSolve1<7>(pp_graph, best, pp);
  } else if (pp_graph.n() <= 8*BITS) {
    DoSolve1<8>(pp_graph, best, pp);
  } else if (pp_graph.n() <= 9*BITS) {
    DoSolve1<9>(pp_graph, best, pp);
  } else if (pp_graph.n() <= 10*BITS) {
    DoSolve1<10>(pp_graph, best, pp);
  } else { // Assume that the graph has n <= 640
    assert(0);
  }
  best::PrintBest();
}