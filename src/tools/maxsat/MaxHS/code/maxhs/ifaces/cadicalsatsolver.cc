/***********[cadicalSatSolver.cc]
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

/* Interface to minisat.
 */

#include <algorithm>
#include <ostream>
#include <vector>

#include "minisat/core/SolverTypes.h"

#include "cadical/src/cadical.hpp"
#include "maxhs/ifaces/cadicalsatsolver.h"
#include "maxhs/utils/clause_utils.h"
#include "maxhs/utils/io.h"
#include "maxhs/utils/params.h"

using namespace MaxHS_Iface;
using sat_solver = CaDiCaL::Solver;
using Minisat::mkLit;
using Minisat::var;

/*************************************************/
// Private
int cadicalSolver::var2cadi(const Var v) { return v + 1; }
int cadicalSolver::lit2cadi(const Lit l) {
  int v = var(l) + 1;
  return sign(l) ? -v : v;
}

Var cadicalSolver::cadi2var(const int c) { return abs(c) - 1; }
Lit cadicalSolver::cadi2lit(const int c) {
  int v = cadi2var(c);
  return c < 0 ? mkLit(v, true) : mkLit(v);
}

void cadicalSolver::mark_var_in_solver(const Var v) {
  if (static_cast<size_t>(v) >= v_in_s.size()) {
    v_in_s.resize(v + 1, 0);
    var_fixed_vals.resize(v + 1, l_Undef);
  }
  v_in_s[v] = 1;
}

/*************************************************/
// Solver interface routines

lbool cadicalSolver::solve_(const vector<Lit>& assumps, vector<Lit>& conflict,
                            int64_t confBudget, int64_t propagationBudget) {
  // cout << "cadicalSolver confBudget = " << confBudget << " propagationBudget
  // "
  //      << propagationBudget << "\n";
  // reset_assumptions();
  conflict.clear();
  for (auto l : assumps) {
    // if lit not in solver then it can't impose any restriction when assumed
    mark_lit_in_solver(l);
    sat_solver::assume(lit2cadi(l));
  }
  limit("conflicts", confBudget);
  limit("propagations", propagationBudget);
  int res = sat_solver::solve();
  if (res == 0) return l_Undef;
  if (res == 10) return l_True;

  for (auto l : assumps)
    // find all lits in conflict (they are all negated assumptions)
    if (lit_in_solver(l) && in_conflict(~l)) conflict.push_back(~l);
  return l_False;
}

class Time_terminator : public CaDiCaL::Terminator {
 public:
  Time_terminator(double st) : stop_time{st} {}
  ~Time_terminator() = default;
  bool terminate() {
    if (cpuTime() >= stop_time) {
      return true;
    } else {
      return false;
    }
  }
  double stop_time;
};

lbool cadicalSolver::solve_(const vector<Lit>& assumps, vector<Lit>& conflict,
                            double timeLimit) {
  // cout << "cadicalSolver confBudget = " << confBudget << " propagationBudget
  // "
  //      << propagationBudget << "\n";
  // reset_assumptions();
  log(2, "called cadicalsolve time budget = ", time_fmt(timeLimit), '\n');

  conflict.clear();
  for (auto l : assumps) {
    // if lit not in solver then it can't impose any restriction when assumed
    mark_lit_in_solver(l);
    sat_solver::assume(lit2cadi(l));
  }
  limit("conflicts", -1);
  limit("propagations", -1);
  auto time_term{Time_terminator(cpuTime() + timeLimit)};
  connect_terminator(&time_term);
  int res = sat_solver::solve();
  disconnect_terminator();

  log(2, "finished cadicalsolve res = ", res, "\n");

  if (res == 0) return l_Undef;
  if (res == 10) return l_True;

  for (auto l : assumps)
    // find all lits in conflict (they are all negated assumptions)
    if (lit_in_solver(l) && in_conflict(~l)) conflict.push_back(~l);
  return l_False;
}

bool cadicalSolver::in_conflict(Lit l) {
  // l is in the conflict if ~l is a failed assumption
  return (sat_solver::failed(lit2cadi(~l)));
}

lbool cadicalSolver::simplify(int rounds) {
  int res = sat_solver::simplify(rounds);
  if (res == 10) return l_True;
  if (res == 20)
    return l_False;
  else
    return l_Undef;
}

lbool cadicalSolver::no_elim_simplify() {
  int res = sat_solver::no_elim_simplify();
  if (res == 10) return l_True;
  if (res == 20)
    return l_False;
  else
    return l_Undef;
}

bool cadicalSolver::findImplications(const vector<Lit>& assumps,
                                     vector<Lit>& imps) {
  vector<int> iassumps, iimps;
  for (auto lit : assumps) {
    if (lit_in_solver(lit)) iassumps.push_back(lit2cadi(lit));
  }
  auto res = sat_solver::find_up_implicants(iassumps, iimps);
  imps.clear();
  if (res)
    for (auto ilit : iimps) imps.push_back(cadi2lit(ilit));
  // else
  // cout << "findImplications bad result " << assumps << " res " << res <<
  // '\n';
  return res;
}

bool cadicalSolver::unit_propagate() {
  auto res = sat_solver::root_level_propagate();
  return res;
}

void cadicalSolver::addClause(const vector<Lit>& lts) {
  for (auto lt : lts) {
    mark_lit_in_solver(lt);
    add(lit2cadi(lt));
  }
  add(0);
}

void cadicalSolver::setPolarity(lbool pol, Var v) {
  // cout << "setPolarity(" << v << "," << pol << ")\n";
  if (!var_in_solver(v)) return;
  int cvar = var2cadi(v);
  if (pol == l_Undef)
    sat_solver::unphase(cvar);
  else
    sat_solver::phase(pol == l_False ? -cvar : cvar);
  // cout << "cadi = " << var2cadi(v)
  //      << " cadilit = " << lit2cadi(mkLit(v, pol == l_True ? true : false))
  //      << "\n";
}

lbool cadicalSolver::modelValue(Lit p) const {
  if (!lit_in_solver(p) || status() != 10) return l_Undef;
  int res = val(lit2cadi(p));
  return res > 0 ? l_True : l_False;
}

lbool cadicalSolver::modelValue(Var x) const {
  if (!var_in_solver(x) || status() != 10) return l_Undef;
  int res = val(var2cadi(x));
  return res > 0 ? l_True : l_False;
}

lbool cadicalSolver::fixedValue(const Var x) const {
  if (!var_in_solver(x)) return l_Undef;
  int res = sat_solver::fixed(var2cadi(x));
  return !res ? l_Undef : (res > 0 ? l_True : l_False);
}

lbool cadicalSolver::fixedValue(const Lit p) const {
  if (!lit_in_solver(p)) return l_Undef;
  int res = sat_solver::fixed(lit2cadi(p));
  return !res ? l_Undef : (res > 0 ? l_True : l_False);
}

Lit cadicalSolver::eqLit(Lit p) {
  if (!lit_in_solver(p)) return p;
  return cadi2lit(sat_solver::eqLit(lit2cadi(p)));
}

bool cadicalSolver::get_ith_clause(vector<Lit>& cls, size_t i) const {
  vector<int> icls;
  cls.clear();
  auto res = sat_solver::get_ith_clause(icls, i);
  for (auto ilt : icls) cls.push_back(cadi2lit(ilt));
  return res;
}

// wrap interface level FnOnClause class in cadical specific clause
// iterators

class cadical_ClauseIterator : public CaDiCaL::ClauseIterator {
 public:
  cadical_ClauseIterator(FnOnClause& f) : fn{f} {}
  bool clause(const std::vector<int>& icls) {
    cls.clear();
    for (auto ilt : icls) cls.push_back(cadicalSolver::cadi2lit(ilt));
    return fn.on_clause(cls);
  }

 private:
  FnOnClause& fn;
  vector<Lit> cls;
};

class cadical_EX_ClauseIterator : public CaDiCaL::WitnessIterator {
  // the extension clauses are not always reduced by units and
  // equalities, if, e.g., these clauses were added to the extension
  // stack before a unit was discovered. So the returned clauses might
  // need to be post-processed

 public:
  cadical_EX_ClauseIterator(FnOnClause& f) : fn{f} {}

  bool witness(const std::vector<int>& icls,
               [[maybe_unused]] const std::vector<int>& witness) {
    cls.clear();
    for (auto ilt : icls) cls.push_back(cadicalSolver::cadi2lit(ilt));
    return fn.on_clause(cls);
  }

 private:
  FnOnClause& fn;
  vector<Lit> cls;
};

bool cadicalSolver::reduce_clause_by_eqs_and_fixed_values(vector<Lit>& cls) {
  // clauses on extension stack have not necessarily been reduced by
  // units nor by equality substitutions. Do those steps on the clause
  // return false if clause is satisfied (no longer have a clause)
  // else return true (clause returned might be empty)
  size_t j{0}, i{0};
  for (; i < cls.size(); ++i) {
    auto eq = eqLit(cls[i]);
    auto tval = fixedValue(eq);
    if (tval == l_True) {
      cls.clear();
      return false;
    } else if (tval == l_Undef) {
      cls[j++] = eq;
    }
  }
  cls.resize(j);
  return true;
}

std::pair<MaxHS::Packed_vecs<Lit>, vector<Lit>>
cadicalSolver::get_all_clauses() {
  MaxHS::Packed_vecs<Lit> all_cls;
  vector<Lit> units;
  vector<Lit> new_cls;
  auto ex_clauses_fn = [&](vector<Lit>& cls) {
    if (reduce_clause_by_eqs_and_fixed_values(cls) &&
        MaxHS::rm_dups_and_check_tautology(cls)) {
      if (cls.size() == 1)
        units.push_back(cls[0]);
      else
        all_cls.addVec(cls);
    }
    return true;
  };
  FnTraverse collect_ex{ex_clauses_fn};
  log(1, "c getting all clauses from sat solver\n");
  traverse_extension_clauses(collect_ex);
  struct {
    size_t ncls, nunits;
  } before, after;
  before = {all_cls.size(), units.size()};
  log(1, "c obtained ", before.ncls, " clauses and ", before.nunits,
      " units from extension stack\n");
  auto in_clauses_fn = [&](vector<Lit>& cls) {
    // internal clauses are already equality and fixedvalue reduced by cadical
    if (cls.size() == 1)
      units.push_back(cls[0]);
    else
      all_cls.addVec(cls);
    return true;
  };
  FnTraverse collect_in{in_clauses_fn};
  traverse_solver_clauses(collect_in);
  after = {all_cls.size(), units.size()};
  log(1, "c obtained ", after.ncls - before.ncls, " clauses and ",
      after.nunits - before.nunits, " units from internal clauses\n");
  log(1, "c obtains ", after.ncls, " clauses and ", after.nunits,
      " units from sat solver\n");

  return {all_cls, units};
}

bool cadicalSolver::traverse_solver_clauses(FnOnClause& fn) const {
  cadical_ClauseIterator ci{fn};
  return traverse_clauses(ci);
}

bool cadicalSolver::traverse_extension_clauses(FnOnClause& fn) const {
  cadical_EX_ClauseIterator ci{fn};
  return traverse_witnesses_forward(ci);
}

// This can be called whenever the sat solver could fix more
// variables, e.g., after solve.

void cadicalSolver::updateFixed() {
  if (n_internal_fixed() == n_ifixed) {
    // nothing new was fixed in solver since last call.
    return;
  }
  auto prev_fixed{fixed_lits.size()};
  const vector<int>& ex_units = solver_units();
  for (; n_ifixed < ex_units.size(); n_ifixed++) {
    Lit fixedLit = cadi2lit(ex_units[n_ifixed]);
    if (var_fixed_vals[var(fixedLit)] == l_Undef) {
      fixed_lits.push_back(fixedLit);
      var_fixed_vals[var(fixedLit)] = sign(fixedLit) ? l_False : l_True;
    }
  }

  if (params.verbosity > 2 && prev_fixed < fixed_lits.size())
    cout << "c cadical found " << fixed_lits.size() - prev_fixed
         << " forced lits\n";
  // auto prev1_fixed{fixed_lits.size()};

  // Units_in_witnesses uw{this};
  // sat_solver::traverse_witnesses_backward(uw);

  // if (params.verbosity > 2 && prev1_fixed < fixed_lits.size())
  //   cout << "c cadical extension stack yielded " << fixed_lits.size() -
  //   prev1_fixed
  //        << " additional forced lits\n"
  //        << "fixed_lits.size() = " << fixed_lits.size() << '\n';
}

vector<Lit> cadicalSolver::getForced(int index) {
  updateFixed();
  return vector<Lit>(fixed_lits.begin() + index, fixed_lits.end());
}
