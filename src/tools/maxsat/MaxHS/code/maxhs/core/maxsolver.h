/***********[maxsolver.h]
Copyright (c) 2012-2013 Jessica Davies, Fahiem Bacchus

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

***********/

#ifndef MaxSolver_h
#define MaxSolver_h

#include <iostream>
#include <string>
#include <vector>

#ifdef GLUCOSE
#include "glucose/core/SolverTypes.h"
#else
#include "minisat/core/SolverTypes.h"
#endif

#include "maxhs/core/bvars.h"
#include "maxhs/core/maxsolvertypes.h"
#include "maxhs/ds/packed.h"

class Wcnf;
class Assumps;
class SumManager;

namespace MaxHS_Iface {
class Cplex;
class SatSolver;
class Muser;
class GreedySolver;
}  // namespace MaxHS_Iface

namespace MaxHS {

// note that each Seed_type implies the previous, so Seed_type::all implies
// cores, non_cores mixed, and ordinary
  enum class Seed_type { none = 0, cores, non_cores, mixed, all, all_no_limit };

class MaxSolver {
 public:
  MaxSolver(Wcnf* f);
  ~MaxSolver();
  void solve();  // Solve the initialized Wcnf
  Weight UB();
  Weight LB() { return lower_bnd; }
  bool isSolved() { return solved; }
  bool isUnsat() { return unsat; }
  const std::vector<lbool>& getBestModel() { return UBmodel; }
  const Wcnf* getWcnf() { return theWcnf; }

  // Not private but not for general users.
  // Used by interrupt trap
  void printStatsAndExit(int signum, int exitType);

  // Called by the SAT solver to tell whether or not a learnt clause
  // can be deleted
  bool deleteLearntTest(const std::vector<Lit>& c) const;
  Weight updateUB();
  Weight updateUB_from_true_lits(vector<Lit>&);
  void updateLB(Weight wt) {
    if (wt > lower_bnd) lower_bnd = wt;
  }
  bool check_termination();

 protected:
  Wcnf* theWcnf;
  Bvars bvars;
  MaxHS_Iface::SatSolver* satsolver{};
  MaxHS_Iface::GreedySolver* greedysolver{};
  MaxHS_Iface::Muser* muser{};
  MaxHS_Iface::Cplex* cplex{};
  SumManager* summations{};
  int top_freq{};

  // solvers
  bool solve_hitting_set();
  bool solve_lsu();
  bool solve_unwt_lsu();
  bool solve_wt_lsu();

  // BOUNDS
  Weight sat_wt{};  // wt of soft clauses known to be satisfiable.
  Weight forced_wt{};

  Weight lower_bnd{};  // lower bound on wt false soft clauses.
  Weight absGap;  // Stop when (UB() - LB()) <= absGap (absGap will be zero for
                  // integer weights).
  bool check_termination(const std::string& location);
  Weight getForcedWt();  // wt of softs falsified by forced units
  Weight getSatClsWt();  // wt of satisfied softs in current model
  Weight checkSolWt(const std::vector<Lit>& sol);
  Weight getWtOfBvars(const std::vector<Lit>& blits);  // sum a vector of wts.
  std::vector<lbool> UBmodel;  // Holds the current best (upper bound) model
  std::vector<lbool>
      UBmodelSofts;  // Holds true/false for satisfied softs in the UB model.
  std::vector<lbool> tmpModelSofts;  // Temp vector for holding satisfied status
                                     // of softs in latest model.

  bool haveUBModel{false};
  bool have_new_UB_model{false};
  void setUBModel();

  // SAT solver interaction
  //  set up bvars
  void configBvar(Var, MaxHS_Iface::SatSolver*);

  // add clauses to sat solver
  void addHards(MaxHS_Iface::SatSolver*);
  void addSofts(MaxHS_Iface::SatSolver*);

  // b-var equivalent clauses for Fbeq (potentially removable).
  // removable b-vars implemented by adding extra control variable to all eq
  // clauses. If the eq clauses are specified as being removable one must assume
  // the literal returned by "activateSoftEqLit" during sat solving. rmSoftEqs
  // removes "removable" eqs.
  void addSoftEqs(MaxHS_Iface::SatSolver*, bool removable);
  void addSoftEqs(MaxHS_Iface::SatSolver*, bool removable,
                  const std::vector<int>& indicies);
  Var eqCvar{};  // control variable for b-var equivalences.
  Lit eqCvarPos{Minisat::lit_Undef};  // Postive literal of eqCvar
  void rmSoftEqs(MaxHS_Iface::SatSolver* slv) {
    // assert eqCvar and thus remove all eq clauses after
    // a call to simplify. Can't add eq clauses again
    if (eqCvarPos != Minisat::lit_Undef) slv->addClause(eqCvarPos);
  }
  Lit activateSoftEqLit() {
    // return assumption that activates the soft equivalent clauses
    return mkLit(eqCvar, true);
  }

  // If not using Fbeq allow forcing negated b-vars after each sat call.
  void satSolverAddBvarsFromSofts();  // with Fb force negated b-vars when soft
                                      // is satisfied.
  void initSftSatisfied();  // for var --> satisfied softs map for forcing neg
                            // b-vars in Fb

  // Manage Cplex and Greedy solver Clauses.
  std::vector<vector<Lit>> dsjnt_cores;
  std::vector<vector<Lit>> find_seeding_clauses(Seed_type seed_type);
  bool allClausesSeeded{false};

  std::vector<int> bLitOccur;
  Packed_vecs<Lit> cplexClauses{};
  Packed_vecs<Lit> greedyClauses{};
  Packed_vecs<int>
      sftSatisfied{};  // map from ordinary var --> satisfied sft clauses

  int cplexAddNewClauses();
  int greedyAddNewClauses();
  void cplexAddCls(std::vector<Lit>&& cls);
  void store_cplex_greedy_cls(const std::vector<Lit>& cls);

  // Transfer units between sub-solvers.
  int greedyAddNewForcedBvars();
  int cplexAddNewForcedBvars();

  struct GetNewUnits {
    int not_seen{0};
    std::vector<Lit> forced;
    void update(MaxHS_Iface::SatSolver* slv) {
      forced = slv->getForced(not_seen);
      not_seen += forced.size();
    }
    std::vector<Lit>& new_units(MaxHS_Iface::SatSolver* slv) {
      update(slv);
      return forced;
    }
  };

  GetNewUnits greedyNU, satBvarNU, cplexNU;

  // Statistics
  int amountConflictMin{};
  double globalStartTime{};

  // output routines
  void printErrorAndExit(const char* msg);
  void printSolution(const std::vector<lbool>& model);
  void printCurClause(const std::vector<Lit>& cls);
  void reportCplex(Weight cplexLB, Weight solnWt);
  void reportSAT_min(lbool result, double iTime, size_t orig_size, int nMins,
                     double mTime, size_t final_size);
  void reportForced(const std::vector<Lit>& forced, Weight wt);
  void outputConflict(const std::vector<Lit>& conf);
  void optFound(const std::string& reason);
  void unsatFound();
  void checkModel(const std::string& location);
  void printNclausesInSatSolver(const std::string& msg);

  bool printStatsExecuted{false};

  // status flags
  bool solved{false};  // True when finished solving.
  bool unsat{false};   // Hards are UNSAT or no solution of cost < dimacs top

  // internal Subroutines
  void disjointPhase(double time_limit = params.noLimit);
  void allClausesSeeded_maxsat();
  void seqOfSAT_maxsat();
  std::vector<Lit> extract_costBoundinglits(const std::vector<Lit>& soln) const;
  std::vector<Lit> extract_NonCoreBvars(const std::vector<Lit>& soln) const;
  std::vector<Lit> extract_UnValued(const std::vector<Lit>& soln) const;
  bool any_forced_false(const std::vector<Lit>& soln) const;
  bool all_bvars(const std::vector<Lit>& soln) const;

  std::vector<Lit> abstractSolution(const std::vector<Lit>&);
  bool tryHarden_with_lp_soln(const Weight lp_objval,
                              std::vector<double>& lp_redvals,
                              std::vector<Var>& cplex_vars);
  int n_softs_forced_hard_not_in_cplex{};
  int n_softs_forced_hard{};
  int n_softs_forced_false{};

  int n_ovars_forced_true{};
  int n_ovars_forced_false{};
  int n_touts_forced_true{};
  int n_touts_forced_false{};

  bool get_cplex_conflicts(const std::vector<Lit>&, char type, double timeout,
                           int max_cores);
  bool get_greedy_conflicts(double timeout, int max_cores);
  bool get_populate_conflicts(double timeout, int max_cores, double gap);
  bool get_ub_conflicts(double timeout);
  int get_seq_of_conflicts(const std::vector<Lit>&, double first_core_cpu_lim,
                           double other_core_cpu_lim, double timeout,
                           int max_cores);
  bool assumps_are_blocked(const std::vector<Lit>&);

  // mutexes
  void processMutexes();
  const uint8_t inMxdvar = 2;
  const uint8_t inMxbvar = 1;
  std::vector<uint8_t> inMx;
  void setInMutex(Lit b, int type) {
    ensureInMx(b);
    inMx[bvars.toIndex(b)] = (uint8_t)type;
  }
  void ensureInMx(Lit b) {
    if (static_cast<size_t>(bvars.toIndex(b)) >= inMx.size())
      inMx.resize(bvars.toIndex(b) + 1, (uint8_t)0);
  }

  // Noncore stuff
  bool vec_isCore(const std::vector<Lit>& core);

  // Accumulate Cores
  std::vector<Lit> getAssumpUpdates(int sinceCplex, std::vector<Lit>& core);
  std::vector<Lit> greedySoln(bool need_soln);
  std::vector<Lit> fracOfCore(int nCplex, std::vector<Lit>& core);
  Lit maxOccurring(const std::vector<Lit>& core);
  void incrBLitOccurrences(const std::vector<Lit>& core);
  lbool satsolve_min(const Assumps& inAssumps, std::vector<Lit>& outConflict,
                     double sat_cpu_lim, double mus_cpu_lim);
  void minimize_muser(std::vector<Lit>& con, double mus_cpu_lim);
  double m_sum_reduced_frac{};
  double mtime{};
  int mcalls{};
  bool doMin{true};
  bool blit_true_warning{true};

  void check_mus(std::vector<Lit>& con);  // debugging

  // Functor to sort blits for assumptions, minimization, and
  // agumenting the hitting set during Karp.
  struct BLitOrderLt {
    // on function call (b1, b2) return true if b1 satisfies fewer
    // cores per unit weight, or if it satifies the same number of
    // cores per unit weight and its corresponding soft clause is
    // longer.  If b1 < b2 then relaxing b1's soft clause is not as
    // good as relaxing b2's soft clause for making the formula
    // satisfiable at lower cost.

    const std::vector<int>* const occurCount;  // How many cores
                                               // will blit relax.
    Bvars& bvars;  // reference to solver's bvar structure
    bool gt;       // true: > ordering; false < ordering

    BLitOrderLt(const std::vector<int>* const o, Bvars& b, bool g)
        : occurCount{o}, bvars(b), gt{g} {}

    bool operator()(const Lit l1, const Lit l2) const;
  };

  BLitOrderLt blit_lt;
  BLitOrderLt blit_gt;

  // data for tracking triggers for abstraction.
  struct AbsTracking {
    std::vector<Weight> deltas{};
    std::vector<Weight> gaps{};
    std::vector<int> n_constraints{};
    std::vector<bool> did_abstraction{};
  };

  AbsTracking track_abs_triggers{};
  bool do_abstraction() const;
  bool iteration_is_bad(Weight desired_delta, Weight delta, Weight gap) const;

};  // namespace MaxHS

}  // namespace MaxHS

#endif
