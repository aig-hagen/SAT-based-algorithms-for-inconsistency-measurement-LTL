/***********[wcnf.h]
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

#ifndef WCNF_H
#define WCNF_H

#include <iostream>
#include <limits>
#include <string>
#include <vector>

#ifdef GLUCOSE
#include "glucose/core/SolverTypes.h"
#include "glucose/mtl/Heap.h"
#else
#include "minisat/core/SolverTypes.h"
#include "minisat/mtl/Heap.h"
#endif

#include "maxhs/core/maxsolvertypes.h"
#include "maxhs/ds/packed.h"
#include "maxhs/ifaces/satsolver.h"
#include "maxpre/src/preprocessorinterface.hpp"

using std::cout;

#ifdef GLUCOSE
namespace Minisat = Glucose;
#endif

using MaxHS_Iface::SatSolver;
using Minisat::Heap;
using Minisat::l_Undef;
using Minisat::lbool;
using Minisat::Lit;
using Minisat::toInt;
using Minisat::var;

class Bvars;  // Bvars.h loads this header.

class SC_mx {
  /* The blits are such that if they are made
     true the correesponding soft clause is
     relaxed. Thus we incur the cost of
     the soft clause.
  */
 public:
  SC_mx(const vector<Lit> blits, bool is_core, Lit dlit)
      : _blits{blits}, _dlit{dlit}, _is_core{is_core} {}
  // if is_core then
  //  at most one of the blits can be true (at most one of the
  //  corresponding soft clauses can be falsified)
  //  and if dlit is true one of the blits is true.

  // if !is_core then
  //  at most one of the blits can be false (at most one of the
  //  corresponding soft clauses can be satisfied)
  //  and if dlit is false then one of the blits is false.

  const vector<Lit>& soft_clause_lits() const { return _blits; }
  bool is_core() const { return _is_core; }
  Lit encoding_lit() const { return _dlit; }
  // modifying versions
  vector<Lit>& soft_clause_lits_mod() { return _blits; }
  Lit& encoding_lit_mod() { return _dlit; }

 private:
  vector<Lit> _blits;
  Lit _dlit;
  bool _is_core;
};

inline std::ostream& operator<<(std::ostream& os, const SC_mx& mx) {
  os << (mx.is_core() ? "Core Mx: " : "Non-Core-Mx: ")
     << "Defining Lit = " << mx.encoding_lit()
     << " blits = " << mx.soft_clause_lits();
  return os;
}

class Wcnf {
 public:
  bool inputDimacs(const std::string input_name, const std::string file_name,
                   bool preprocess);
  // input clauses
  void set_dimacs_params(int nvars, int nclauses,
                         Weight top = std::numeric_limits<Weight>::max());
  // time to parse the instance
  double parseTime() const { return parsing_time; }
  // total time to parse and preprocess the instance.
  double inputTime() const { return input_time; }

  // api-support for adding hard or soft clauses
  void addHardClause(vector<Lit>& lits);
  void addSoftClause(vector<Lit>& lits, Weight w);
  void addHardClause(Lit p) {
    vector<Lit> tmp{p};
    addHardClause(tmp);
  }
  void addHardClause(Lit p, Lit q) {
    vector<Lit> tmp{p, q};
    addHardClause(tmp);
  }
  void addSoftClause(Lit p, Weight w) {
    vector<Lit> tmp{p};
    addSoftClause(tmp, w);
  }
  void addSoftClause(Lit p, Lit q, Weight w) {
    vector<Lit> tmp{p, q};
    addSoftClause(tmp, w);
  }

  void addSumConst();
  // modify Wcnf. This might add new variables. But all variables
  // in the range 0--dimacs_nvars-1 are preserved. New variables
  // only exist above that range.
  // TODO: for incremental solving we need a map from external to internal
  //      variables so that new clauses with new variables can be
  //      added on top of any added new variables arising from transformations

  void simplify();
  bool orig_all_lits_are_softs() const { return orig_all_lits_soft; }
  bool is_hitting_set_problem() { return is_hitting_set; }

  // print info
  void printFormulaStats();
  void printSimpStats();
  void printFormula(std::ostream& out = std::cout) const;
  void printFormula(Bvars&, std::ostream& out = std::cout) const;
  void printDimacs(std::ostream& out = std::cout) const;

  // data setters and getters mainly for solver
  vector<lbool> rewriteModelToInput(
      const vector<lbool>& ubmodel);  // convert model to model of input formula

  //  check model against input formula. Only available if option
  //  chkSoln is used (on by default). Since we can read instance from
  //  stdin we must store a copy of the instance---which can use a
  //  non-trivial amount of memory. This copy is only available if
  //  chkSoln is on.
  Weight checkModel(const vector<lbool>& ubmodel, int& nfalseSofts);

  const MaxHS::Packed_vecs<Lit>& hards() const { return hard_cls; }
  const MaxHS::Packed_vecs<Lit>& softs() const { return soft_cls; }
  const vector<Weight>& softWts() const { return soft_clswts; }
  vector<Lit> getSoft(int i) const { return soft_cls.getVec(i); }
  vector<Lit> getHard(int i) const { return hard_cls.getVec(i); }

  Weight getWt(int i) const { return soft_clswts[i]; }
  size_t softSize(int i) const { return soft_cls.ithSize(i); }
  size_t hardSize(int i) const { return hard_cls.ithSize(i); }

  Weight totalWt() const { return baseCost() + totalClsWt(); }
  Weight totalClsWt() const { return total_cls_wt; }
  Weight baseCost() const { return base_cost; }

  // info about wcnf
  size_t nHards() const { return hard_cls.size(); }
  size_t nSofts() const { return soft_cls.size(); }
  size_t n_non_unit_softs() const;
  size_t nCls() const { return nHards() + nSofts(); }
  //sat solver does not get feed unit softs
  size_t n_clauses_for_solvers() const { return nHards() + n_non_unit_softs(); }
  // including extra variables added via transformations
  size_t nVars() const { return maxvar + 1; }

  Var maxVar() const {
    return maxvar;
  }  // Users should regard this as being the number of vars

  Weight minSftWt() const { return wt_min; }
  Weight maxSftWt() const { return wt_max; }
  int nDiffWts() const { return diffWts.size(); }
  bool is_unweighted() const { return diffWts.size() == 1;}
  bool is_weighted() const {return diffWts.size() > 1; }
  const vector<Weight>& getDiffWts() { return diffWts; }
  const vector<int>& getDiffWtCounts() { return diffWtCounts; }
  const vector<Weight>& getTransitionWts() { return transitionWts; }
  Weight aveSftWt() const { return wt_mean; }
  Weight varSftWt() const { return wt_var; }

  bool isUnsat() { return unsat; }
  bool integerWts() { return intWts; }
  const std::string& fileName() const { return instance_file_name; }

  // mutexes
  const vector<SC_mx>& get_SCMxs() const { return mutexes; }
  int n_mxes() const { return mutexes.size(); }
  const SC_mx& get_ith_mx(int i) const { return mutexes[i]; }
  int ith_mx_size(int i) const { return mutexes[i].soft_clause_lits().size(); }

  // get input file literal
  Lit input_lit(Lit l) const {
    if (static_cast<size_t>(var(l)) >= in2ex.size() ||
        in2ex[var(l)] == Minisat::var_Undef)
      return Minisat::lit_Undef;
    return Minisat::mkLit(in2ex[var(l)], sign(l));
  }

  vector<Lit> vec_to_file_lits(const vector<Lit>& v) {
    vector<Lit> fv;
    for (auto l : v) fv.push_back(input_lit(l));
    return fv;
  }

 private:
  bool parse_DIMACS(std::istream& inStream);
  void parse_clause(const std::string& ln);
  void parse_p_line(const std::string& ln);
  // Weight get_wt(std::istringstream& ss);
  // void get_lits(vector<Lit>& lits, std::istringstream& ss);

  // input clauses from file (as +/- ints)

  Weight get_wt(const char*& cur_pos, const std::string& ln);
  void get_file_lits(vector<int>& lits, const char*& cur_pos,
                     const std::string& ln);
  bool seen_clause, seen_p_line, isWcnf;
  vector<int> input_cls;
  void add_file_hard_clause(const vector<int>&);
  void add_file_soft_clause(const vector<int>&, Weight);
  void add_file_dimacs_clause(const vector<int>&, Weight);

  void update_maxorigvar(const vector<Lit>& lits);
  void update_maxorigvar(Lit lt);
  void _addHardClause(vector<Lit>& lits);
  void _addHardClause(Lit p) {
    vector<Lit> tmp{p};
    _addHardClause(tmp);
  }
  void _addHardClause(Lit p, Lit q) {
    vector<Lit> tmp{p, q};
    _addHardClause(tmp);
  }
  void _addSoftClause(vector<Lit>& lits, Weight w);
  void _addSoftClause(Lit p, Weight w) {
    vector<Lit> tmp{p};
    _addSoftClause(tmp, w);
  }
  void _addSoftClause(Lit p, Lit q, Weight w) {
    vector<Lit> tmp{p, q};
    _addSoftClause(tmp, w);
  }
  void computeWtInfo();

  // preprocessing routines
  bool simplify_has_been_run{false};
  void subEqsAndUnits();
  MaxHS::Packed_vecs<Lit> reduce_by_eqs_and_units(
      SatSolver* sat_solver, const MaxHS::Packed_vecs<Lit>& cls, bool softs);
  void rm_subsumed_softs();
  bool test_all_lits_are_softs();
  bool test_if_input_is_hitting_set();
  void feed_sat_solver_hards(SatSolver* const);
  void simpleHarden();
  void mxBvars();
  vector<vector<Lit>> mxFinder(Bvars&);
  void processMxs(const vector<vector<Lit>>, Bvars&);
  void remapVars();
  Var maxOrigVar() const { return maxorigvar; }        // input variables
  size_t nOrigVars() const { return maxorigvar + 1; }  // are for private use.
  maxPreprocessor::PreprocessorInterface* pif{};
  void try_maxpre();
  void use_orig_instance();
  bool using_maxpre{false};

  // Packed_vecs<Lit> reduce_by_units(Packed_vecs<Lit>& cls, SatSolver_uniqp&
  // slv,
  //                                  bool softs);
  int maxorigvar{}, maxvar{};
  int dimacs_nvars{};
  int dimacs_nclauses{};
  int nfalse_softs{};
  double parsing_time{};
  double input_time{};
  Weight total_cls_wt{};  // Weight of soft clauses after simplifications.
  Weight base_cost{};
  // weight of a hardclause...typically sum of soft clause weights + 1
  Weight dimacs_top{std::numeric_limits<Weight>::max()};
  Weight wt_var{}, wt_mean{}, wt_min{}, wt_max{};
  std::string instance_file_name;
  bool unsat{false};
  bool intWts{true};
  bool orig_all_lits_soft{false};
  bool is_hitting_set{false};
  vector<Weight> diffWts;
  vector<int> diffWtCounts;

  vector<char> orig_unit_soft;
  vector<Weight> transitionWts;  // weights w s.t. sum of soft clauses with
                                 // weight less that w is less than w
  MaxHS::Packed_vecs<Lit> hard_cls;
  MaxHS::Packed_vecs<Lit> soft_cls;
  vector<Weight> soft_clswts;
  // store unchanged copy of input formula if verifying (stored in
  // original file format--as + and negative numbers)
  MaxHS::Packed_vecs<int> file_hards;
  MaxHS::Packed_vecs<int> file_softs;
  vector<Weight> file_wts;

  // store preprocessing computation for remaping.
  int nOrigUnits{};

  vector<lbool> fixed_tvals;
  vector<Lit> eqLit;
  void set_up_eq_and_fixed_maps();

  lbool fixed_value(Lit l) {
    return sign(l) ? fixed_value(var(l)).neg() : fixed_value(var(l));
  }
  lbool fixed_value(Var v) {
    if (static_cast<size_t>(v) >= fixed_tvals.size())
      return l_Undef;
    else
      return fixed_tvals[v];
  }
  void set_fixed_var(Var v, lbool tv) {
    if (static_cast<size_t>(v) < fixed_tvals.size()) fixed_tvals[v] = tv;
  }
  void set_fixed(Lit l) {
    lbool tv = sign(l) ? l_False : l_True;
    set_fixed_var(var(l), tv);
  }
  void set_fixed_equal(Lit l, Lit e) {
    if (sign(l) == sign(e))
      set_fixed_var(var(l), fixed_value(var(e)));
    else
      set_fixed_var(var(l), fixed_value(var(e)).neg());
  }
  Lit get_eqLit(Lit l) {
    if (static_cast<size_t>(toInt(l)) >= eqLit.size())
      return l;
    else
      return eqLit[toInt(l)];
  }
  void set_eqLit(Lit l, Lit e) {
    if (static_cast<size_t>(toInt(l)) < eqLit.size()) {
      eqLit[toInt(l)] = e;
      eqLit[toInt(~l)] = ~e;
    }
  }

  // for remapping the variables
  vector<char> flipped_vars;  // convert unit softs to contain positive lit.
                              // must remove dups first
  vector<Var> ex2in, in2ex;
  Lit map_in2ex(Lit l) const {
    assert(static_cast<size_t>(var(l)) < in2ex.size() &&
           in2ex[var(l)] != Minisat::var_Undef);
    return Minisat::mkLit(in2ex[var(l)], sign(l));
  }

  struct SftClsData {
    uint32_t index;
    uint32_t hash;
    Weight w;  // weight == 0 ==> redundant
    SftClsData(uint32_t i, uint32_t h, Weight wt) : index{i}, hash{h}, w{wt} {}
  };

  void initSftClsData(vector<SftClsData>& cdata);
  bool eqSoftCls(const SftClsData& a, const SftClsData& b);

  // mutexes
  vector<SC_mx> mutexes;
};

#endif
