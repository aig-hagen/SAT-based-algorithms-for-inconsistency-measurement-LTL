/***********[cplex.h]
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

#ifndef CPLEX_H
#define CPLEX_H

#include <ilcplex/cplexx.h>
#include <map>
#include <utility>
#include <vector>

#ifdef GLUCOSE
#include "glucose/core/SolverTypes.h"
#include "glucose/utils/System.h"
#else
#include "minisat/core/SolverTypes.h"
#include "minisat/utils/System.h"
#endif

#include "maxhs/core/bvars.h"
#include "maxhs/core/maxsolvertypes.h"
#include "maxhs/core/summanager.h"
#include "maxhs/core/wcnf.h"

#ifdef GLUCOSE
namespace Minisat = Glucose;
#endif

using Minisat::var_Undef;
using std::vector;

namespace MaxHS_Iface {

enum class Cplex_solver_type { with_lp, no_lp };
enum class Cplex_callback_type {
  no_callback,
  stop_less_than_ub,
  stop_less_than_equal_ub
};
enum class Cplex_soln_status { error, no_solution, feasible, optimal };

class Cplex {
 public:
  Cplex(Bvars& b, SumManager* s, vector<lbool>& ubModelSofts,
        vector<lbool>& ubModel, bool integerWts,
        Cplex_solver_type solver_type = Cplex_solver_type::with_lp);
  ~Cplex();

  std::pair<Weight, Cplex_soln_status> solveBudget(
      vector<Lit>& solution, double UB, double timeLimit,
      Cplex_callback_type callback_type =
          Cplex_callback_type::stop_less_than_ub) {
    // Find best solution within time and call back limits.  return
    // <wt, Cplex_soln_status> a pair, where wt is the objective cost
    // of the best solution found (if one was found) and
    // Cplex_soln_status is set appropriately. If a solution was
    // returned (Cplex_soln_status = optimal or feasible) in solution
    // return a setting of all bvars in best solution. If the bvar is
    // not mentioned in the cplex problem, return the literal that has
    // zero weight in the solution. Also return the setting of all
    // ordinary variables that have been added to cplex.
    auto stime = cpuTime();
    prevTotalTime = totalTime;

    auto val = solve_(solution, UB, timeLimit, callback_type);
    totalTime += cpuTime() - stime;
    numSolves++;
    return val;
  }

  std::pair<Weight, Cplex_soln_status> solve(
      vector<Lit>& solution, double UB,
      Cplex_callback_type callback_type =
          Cplex_callback_type::stop_less_than_ub) {
    double timeLimit{-1.0};
    auto ret = solveBudget(solution, UB, timeLimit, callback_type);
    return ret;
  }

  // try to populate cplex solutions using a time limit...don't return
  // them return the number of solutions.
  int populate(double timeLimit, double gap);
  int populate(double gap) { return populate(-1, gap); }
  void getPopulatedSolution(int, vector<Lit>&);

  // solve the lp relaxation...return objective value of LP relaxation
  // vector of solution weights and reduced costs indexed by soft clause index
  // e.g., solution[i] = weight of blit of i-th soft clause.
  Weight solve_lp_relaxation(vector<double>& solution,
                             vector<double>& reduced_costs,
                             vector<Var>& cplex_vars) {
    auto startTime = cpuTime();
    auto val = solve_lp_relaxation_(solution, reduced_costs, cplex_vars);
    totalLPTime += cpuTime() - startTime;
    numLPSolves++;
    return val;
  }

  bool is_valid() { return solver_valid; }
  bool filter_by_units(vector<Lit>& theCon);
  void add_clausal_constraint(vector<Lit>& con);
  void add_processed_clause(const vector<Lit>& theCon);

  bool add_mutex_constraint(const SC_mx& mx);
  bool add_sum_constraint(vector<Lit>&, int64_t);
  bool var_in_cplex(Var v) { return ex2in(v) != var_Undef; }
  bool lit_in_cplex(Lit l) { return var_in_cplex(var(l)); }

  // stats
  int nCnstr() { return numConstraints; }
  uint64_t totalCnstrSize() { return totalConstraintSize; }
  int nNonCores() { return numNonCoreConstraints; }
  uint64_t totalNonCore() { return totalNonCoreSize; }

  double solveTime() { return totalTime - prevTotalTime; }
  double total_time() { return totalTime; }

  double total_lp_time() { return totalLPTime; }
  int nSolves() { return numSolves; }
  int nLPSolves() { return numLPSolves; }

  // public for call back access --- friend should work?
  void processError(int status, bool terminal, const char* msg);

 protected:
  Bvars& bvars;
  SumManager* summations;

  vector<lbool>& ubModelSofts;
  vector<lbool>& ubModel;

  static CPXENVptr env;
  bool init_env();

  CPXLPptr mip{nullptr};
  // use this lp to do lp-relaxation and trial hardening
  CPXLPptr linp{nullptr};
  bool solver_valid{true};
  bool intWts;
  double LB{};
  double absGap;

  void add_sum_unit(Lit tout);
  void add_sum_output_defn(Lit lt);
  void add_sum_output_constraint(Lit lt);
  // bool add_sum_atmost(Lit lt);

  // forced units (in external ordering)
  vector<lbool> exUnits;
  void setExUnits(Lit l);
  lbool getExUnits(Lit l);

  // main processing code
  void addNewVar(Var ex);
  std::pair<Weight, Cplex_soln_status> getSolution(vector<Lit>& solution,
                                                   bool optimal);

  std::pair<Weight, Cplex_soln_status> solve_(
      vector<Lit>& solution, double UB, double timeLimit,
      Cplex_callback_type callback_type);

  Weight solve_lp_relaxation_(vector<double>& solution,
                              vector<double>& reduced_costs,
                              vector<Var>& cplex_vars);

  // internal cplex routines
  void writeCplexModel();
  void useBestModelStart(CPXLPptr);

  // Stats
  int numSolves{};
  double totalTime{}, prevTotalTime{};
  int numLPSolves{};
  double totalLPTime{};

  int numConstraints{}, numNonCoreConstraints{};
  uint64_t totalConstraintSize{0}, totalNonCoreSize{0};

  // External to Internal Mapping
  vector<Var> in2ex_map;
  vector<int> ex2in_map;

  void ensure_mapping(const Var ex);
  void ensure_mapping(const Lit lt) { ensure_mapping(var(lt)); }

  Var in2ex(int v) const {
    if (v >= (int)in2ex_map.size())
      return var_Undef;
    else
      return in2ex_map[v];
  }

  // In most applications every internal variable of the Cplex Solver
  // is associated with an external literal on creation.
  // So this array function is safe...i.e., won't add var_Undef to output
  // vector. An array version of ex2in is typically not safe in this way.
  // so is not provided.
  void in2ex(const vector<int>& from, vector<Var>& to) const {
    to.clear();
    for (size_t i = 0; i < from.size(); i++) to.push_back(in2ex(from[i]));
  }

  int ex2in(Var v) const {
    if (v >= static_cast<int>(ex2in_map.size()))
      return var_Undef;
    else
      return ex2in_map[v];
  }

  int ex2in(Lit lt) const { return ex2in(var(lt)); }
};
}  // namespace MaxHS_Iface
#endif
