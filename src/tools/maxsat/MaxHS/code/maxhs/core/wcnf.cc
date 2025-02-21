/***********[wcnf.cc]
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

#include <charconv>
#include <cmath>
#include <set>

#ifdef GLUCOSE
#include "glucose/core/SolverTypes.h"
#include "glucose/utils/System.h"
#else
#include "minisat/core/SolverTypes.h"
#include "minisat/utils/System.h"
#endif

#include "maxhs/core/bvars.h"
#include "maxhs/core/wcnf.h"
#include "maxhs/ds/onewatch.h"
#include "maxhs/ifaces/cadicalsatsolver.h"
#include "maxhs/utils/clause_utils.h"
#include "maxhs/utils/hash.h"
#include "maxhs/utils/io.h"
#include "maxhs/utils/params.h"
#include "maxhs/utils/parse.h"
#include "maxhs/utils/zfstream.h"

#ifdef GLUCOSE
namespace Minisat = Glucose;
#endif

using MaxHS::rm_dups_and_check_tautology;
using MaxHS_Iface::cadicalSolver;
using MaxHS_Iface::FnTraverse;
using MaxHS_Iface::SatSolver_uniqp;
using Minisat::lbool;
using Minisat::lit_Undef;
using Minisat::mkLit;
using Minisat::sign;
using Minisat::toInt;
using Minisat::toLit;
using Minisat::var_Undef;
using std::cin;
using std::cout;
using std::runtime_error;
using std::string;

constexpr bool parse_DEBUG{false};
constexpr unsigned prepro_DEBUG{false};
constexpr unsigned maxpre_DEBUG{true};

void Wcnf::set_dimacs_params(int nvars, int nclauses, Weight top) {
  dimacs_nvars = nvars;
  dimacs_nclauses = nclauses;
  dimacs_top = top;
}

bool Wcnf::inputDimacs(const std::string input_name,
                       const std::string file_name, bool preprocess) {
  instance_file_name = file_name;
  double start_time = Minisat::cpuTime();

  gzifstream ifs;
  if (!input_name.empty()) {
    ifs.open(input_name, std::ios::in | std::ios::binary);
    if (!ifs.is_open()) {
      cout << "c ERROR: problem opening input file: " << input_name << "\n";
      return false;
    }
  }
  DEBUG_PR<parse_DEBUG>("Parsing from \"", input_name, "\" file_name = \"",
                        file_name, "\"\n");
  if (!parse_DIMACS(input_name.empty() ? cin : ifs)) {
    cout << "c ERROR: Parsing error on input "
         << (input_name.empty() ? "<stdin>" : input_name) << "\n";
    return false;
  }

  is_hitting_set = test_if_input_is_hitting_set();
  if (is_hitting_set) log(1, "c WCNF instance is hitting set\n");

  if (!is_hitting_set && params.use_maxpre &&
      (file_hards.total_size() + file_softs.total_size() <=
       256 * 1024 * 1024) &&
      file_softs.size() >= 5000)
    try_maxpre();
  if (!using_maxpre) use_orig_instance();

  orig_all_lits_soft = test_all_lits_are_softs();
  if (orig_all_lits_soft)
    log(1, "c WCNF all variables are softs in instance\n");

  if (!params.chkSoln) {
    file_hards.clear();
    file_softs.clear();
    vector<Weight> tmp;
    std::swap(file_wts, tmp);
  }

  computeWtInfo();
  parsing_time = Minisat::cpuTime() - start_time;
  printFormulaStats();

  set_up_eq_and_fixed_maps();

  if (preprocess) {
    simplify();
    if (params.verbosity > 0) printSimpStats();
  }
  input_time = Minisat::cpuTime() - start_time;
  return true;
}

bool Wcnf::parse_DIMACS(std::istream& inStream) {
  // throw exception on inStream badbit; auto clean up
  Istream_guard guard(inStream, std::ios_base::badbit);
  string ln;
  seen_clause = seen_p_line = false;
  isWcnf = true;
  try {
    while (std::getline(inStream, ln)) {
      if (ln.empty() || ln.front() == 'c' ||
          ln.find_first_not_of(" \t") == string::npos)
        continue;
      else if (ln.front() == 'p') {
        parse_p_line(ln);
      } else {
        parse_clause(ln);
      }
    }
  } catch (std::runtime_error& ex) {
    cout << "c ERROR: in input " << ex.what();
    return false;
  }
  return true;
}

void Wcnf::parse_p_line(const string& ln) {
  if (seen_p_line) throw runtime_error("Double p-line: \"" + ln + "\"\n");
  if (seen_clause) throw runtime_error("Clauses appeared before p-line\n");
  size_t pos{1};
  while (std::isspace(ln[pos])) ++pos;
  if (ln.rfind("cnf", pos) == pos) {
    isWcnf = false;
    pos += 3;
  } else if (ln.rfind("wcnf", pos) == pos) {
    isWcnf = true;
    pos += 4;
  } else
    throw runtime_error("Invalid pline: \"" + ln + "\"\n");
  auto ptr = ln.c_str() + pos;
  parse_num(ln, ptr, dimacs_nvars);
  parse_num(ln, ptr, dimacs_nclauses);
  if (!isWcnf)
    dimacs_top = std::numeric_limits<Weight>::max();
  else
    try {
      parse_num(ln, ptr, dimacs_top);
    } catch (...) { dimacs_top = std::numeric_limits<Weight>::max(); }
  seen_p_line = true;
  DEBUG_PR<parse_DEBUG>("parsed pline: \"p ", (isWcnf ? "wcnf" : "cnf"),
                        " nvars = ", dimacs_nvars, " ncls = ", dimacs_nclauses,
                        " top = ",
                        (dimacs_top < std::numeric_limits<Weight>::max()
                             ? std::to_string(dimacs_top)
                             : "no top (all clauses soft)\n"));
}

Weight Wcnf::get_wt(const char*& cur_pos, const string& ln) {
  if (isWcnf) {
    double wt;
    parse_num(ln, cur_pos, wt);
    return wt;
  } else
    return 1;
}

void Wcnf::get_file_lits(vector<int>& lits, const char*& cur_pos,
                         const string& ln) {
  int64_t lit;
  lits.clear();
  while (true) {
    parse_num(ln, cur_pos, lit);
    if (lit == 0)
      break;
    else
      lits.push_back(lit);
  }
}

void Wcnf::parse_clause(const string& ln) {
  auto cur_pos = ln.c_str();
  if (!seen_p_line) {
    // parse new format with hard clauses marked with initial h
    if (ln.front() == 'h' || ln.front() == 'H') {
      get_file_lits(input_cls, ++cur_pos, ln);
      add_file_hard_clause(input_cls);
    } else {
      Weight wt;
      wt = get_wt(cur_pos, ln);
      get_file_lits(input_cls, cur_pos, ln);
      add_file_soft_clause(input_cls, wt);
    }
  } else {
    // old format with weights >= top marking hard clauses
    Weight wt;
    wt = get_wt(cur_pos, ln);
    get_file_lits(input_cls, cur_pos, ln);
    add_file_dimacs_clause(input_cls, wt);
  }
  seen_clause = true;
}

void Wcnf::add_file_dimacs_clause(const vector<int>& lits, Weight w) {
  // This routine needs to know dimacs_top (so set_dimacs_params should
  // have been called first) to determine if the clause is soft or hard
  if (w >= dimacs_top)
    add_file_hard_clause(lits);
  else
    add_file_soft_clause(lits, w);
}

void Wcnf::add_file_hard_clause(const vector<int>& lits) {
  file_hards.addVec(lits);
}

void Wcnf::add_file_soft_clause(const vector<int>& lits, Weight w) {
  file_softs.addVec(lits);
  file_wts.push_back(w);
}

bool Wcnf::test_if_input_is_hitting_set() {
  // tests input file clauses (i.e., before conversion to Lits)
  std::set<int> soft_lits;
  for (const auto sc : file_softs) {
    if (sc.size() != 1) return false;
    soft_lits.insert(sc[0]);
  }

  for (const auto hc : file_hards)
    for (auto l : hc)
      if (soft_lits.find(-l) == soft_lits.end()) return false;
  return true;
}

void Wcnf::use_orig_instance() {
  log(1, "c Using original input instance\n");
  vector<Lit> cls;
  for (auto hc : file_hards) {
    cls.clear();
    for (auto l : hc) {
      auto var = abs(l) - 1;
      cls.push_back((l > 0) ? mkLit(var) : mkLit(var, true));
    }
    addHardClause(cls);
  }
  int i{0};
  for (auto sc : file_softs) {
    cls.clear();
    for (auto l : sc) {
      auto var = abs(l) - 1;
      cls.push_back((l > 0) ? mkLit(var) : mkLit(var, true));
    }
    addSoftClause(cls, file_wts[i++]);
  }
}

void Wcnf::update_maxorigvar(const vector<Lit>& lits) {
  for (auto l : lits) update_maxorigvar(l);
}

void Wcnf::update_maxorigvar(Lit l) {
  if (var(l) > maxorigvar) maxorigvar = var(l);
  if (maxorigvar > maxvar) maxvar = maxorigvar;
}

void Wcnf::addHardClause(vector<Lit>& lits) {
  update_maxorigvar(lits);
  if (lits.size() == 1) ++nOrigUnits;
  _addHardClause(lits);
}

void Wcnf::_addHardClause(vector<Lit>& lits) {
  // Use this routine when adding a clause not contained in the
  // original formula, e.g., adding preprocessing clause
  if (lits.empty()) unsat = true;
  if (unsat) return;
  if (!rm_dups_and_check_tautology(lits)) return;  // skip tautologies
  for (auto l : lits)
    if (maxvar < var(l)) maxvar = var(l);
  hard_cls.addVec(lits);
}

void Wcnf::addSoftClause(vector<Lit>& lits, Weight w) {
  // zero weight clauses discarded by this interface function.
  if (w < 0)
    cout << "c ERROR: soft clause cannot have negative weight: " << w << "\n";
  else if (w > 0) {
    update_maxorigvar(lits);
    _addSoftClause(lits, w);
  }
}

void Wcnf::_addSoftClause(vector<Lit>& lits, Weight w) {
  // Use this routine when adding a clause not contained in the
  // original formula, e.g., adding preprocessing clause
  if (unsat) return;
  if (!rm_dups_and_check_tautology(lits)) return;  // skip tautologies
  Weight intPart;
  if (lits.size() > 0) {
    if (modfl(w, &intPart) > 0) intWts = false;
    soft_cls.addVec(lits);
    soft_clswts.push_back(w);
    total_cls_wt += w;
    for (auto l : lits)
      if (maxvar < var(l)) maxvar = var(l);
  } else {
    base_cost += w;
    ++nfalse_softs;
  }
}

size_t Wcnf::n_non_unit_softs() const {
  size_t n{};
  for(const auto cls: soft_cls)
    if(cls.size() > 1) ++n;
  return n;
}

void Wcnf::try_maxpre() {
  vector<vector<int>> clauses;
  vector<uint64_t> weights;
  uint64_t top{1};
  for (auto w : file_wts) top += w;

  for (size_t i = 0; i < file_hards.size(); ++i) {
    clauses.push_back(file_hards.getVec(i));
    weights.push_back(top);
  }
  for (size_t i = 0; i < file_softs.size(); ++i) {
    clauses.push_back(file_softs.getVec(i));
    weights.push_back(file_wts[i]);
  }

  pif = new maxPreprocessor::PreprocessorInterface(clauses, weights, top);
  string techniques = "[bu]#[buvsrgc]";
  auto startTime = cpuTime();
  log(1, "c Running MAXPRE\n");
  pif->preprocess(techniques, 0, params.maxpre_time);
  auto duration = cpuTime() - startTime;
  log(1, "c MAXPRE took ", time_fmt(duration), "s.\n");

  clauses.clear();
  weights.clear();
  vector<int> ignore_labels;
  pif->getInstance(clauses, weights, ignore_labels);
  size_t ns{0}, nh{0};  // count number different types of clauses in
                        // preprocessed instance
  for (auto& w : weights)
    if (w == top)
      ++nh;
    else
      ++ns;
  log(1, "c MAXPRE returns ", ns,
      " soft clauses. Original #soft = ", file_softs.size(), " and ", nh,
      " hard clauses. Original #hards = ", file_hards.size(), '\n');
  if (ns < file_softs.size() * .95) {
    log(1, "c using MAXPRE transformed formula\n");
    using_maxpre = true;
    vector<Lit> cls;
    for (size_t i = 0; i < clauses.size(); ++i) {
      cls.clear();
      for (auto file_lit : clauses[i]) {
        auto var = abs(file_lit) - 1;
        cls.push_back((file_lit > 0) ? mkLit(var) : mkLit(var, true));
      }
      if (weights[i] == top)
        addHardClause(cls);
      else
        addSoftClause(cls, weights[i]);
    }
  } else {
    using_maxpre = false;
    log(1, "c not using MAXPRE formula\n");
    if (pif) delete pif;
  }
}

void Wcnf::set_up_eq_and_fixed_maps() {
  fixed_tvals.clear();
  fixed_tvals.resize(nVars(), l_Undef);
  eqLit.clear();
  eqLit.resize(nVars() * 2);
  for (size_t i = 0; i < nVars(); i++) {
    auto lt = mkLit(i);
    set_eqLit(lt, lt);
  }
}

void Wcnf::simplify() {
  if (simplify_has_been_run)
    log_warn("trying to run WCNF simplify more than once...not supported\n");

  // We allow a Wcnf to transform itself in various model equivalent
  // ways.  Only the remaining hard and soft clauses after
  // simplification will be passed to the solver. After the solver has
  // found a model for the transformed wcnf, it will need to invoke the
  // method 'rewriteModelToInput' to convert that model to a model of
  // the original input formula. (E.g., if eliminated variables must be
  // added back).

  // wcnf_harden --- test some softs can be hardened because of their high
  // weight

  if (params.wcnf_harden) simpleHarden();

  // Look for units and/or equalities implied by the hard clauses.
  // Simplify by hard units and replace y by x if x==y.
  if (params.wcnf_eqs || params.wcnf_units) subEqsAndUnits();

  // New b-variables are not added to soft units, e.g., (x). Instead
  // we reuse the literal in the soft unit as its own b-variable. In this case
  // we need to ensure that we have no duplicate softs for correctness.
  else
    // subEqsAndUnits() already calls rm_subsumed_softs (hence the else
    // condition)
    rm_subsumed_softs();

  // Find groups of mutually exclusive b-variables. Both positive
  // \sum bi <= 1 and negative \sum -bi <= 1.
  //  These are given a special treatment.

  if (params.mx_find_mxes) {
    if (params.mx_find_mxes == 3) {
      params.mx_find_mxes = 2;
      mxBvars();
      params.mx_find_mxes = 1;
      mxBvars();
    } else
      mxBvars();
  }

  // Some of these transformations might increase the forced (base) cost of
  // the WCNF and change the weights
  computeWtInfo();

  // remap the variables
  remapVars();

  if (params.simplify_and_exit) {
    // write simplified wcnf to cout then exit.
    printDimacs(cout);
  }
}

bool Wcnf::test_all_lits_are_softs() {
  // tests clauses after conversion to Lits
  vector<char> is_soft(nVars(), false);
  for (const auto cls : soft_cls) {
    if (cls.size() == 1) is_soft[var(cls[0])] = true;
  }
  for (const auto cls : hard_cls)
    for (auto l : cls)
      if (!is_soft[var(l)]) return false;
  for (const auto cls : soft_cls)
    for (auto l : cls)
      if (!is_soft[var(l)]) return false;
  return true;
}

void Wcnf::simpleHarden() {
  // Check if we can make some soft clauses hard before we do other
  // preprocessing that might be aided by additional hard clauses.
  // This hardening uses the transition weights if there are any.
  // In particular W is a transition weight if the sum of the weights
  // of the clauses with weight less than W is less than W:
  // \sum_{soft clauses c s.t. wt(c) < W} wt(c) < W
  //
  // This means that we would prefer to falsify all clauses of weight
  // less than W rather than a single clause of weight >= W.  We check
  // to see if the set of soft clauses with weight >= W are satisfiable
  //(in conjunction with the hard clauses), i.e., the set
  // HIGHWT = {c| wt(c) >= W or c is hard}.

  // If these clauses HIGHWT are satisfiable then there is no need to
  // ever falsify any clause in HIGHWT, and we can harden these
  // clauses.

  if (unsat) return;
  computeWtInfo();
  log(1, "c transitionWts = ", transitionWts, '\n');
  if (transitionWts.empty()) return;

  SatSolver_uniqp sat_solver{new cadicalSolver};
  feed_sat_solver_hards(sat_solver.get());
  if (unsat) return;

  if constexpr (prepro_DEBUG > 1) {
    DEBUG_PR<(prepro_DEBUG > 1)>("Diffwts = [");
    for (size_t i = 0; i < diffWts.size(); i++)
      DEBUG_PR<(prepro_DEBUG > 1)>("(", diffWts[i], ", ", diffWtCounts[i],
                                   "), ");
    DEBUG_PR<(prepro_DEBUG > 1)>("]\n");
  }

  Weight wt_added{wt_max + 1}, wt_hardened{wt_max + 1};
  for (int i = transitionWts.size() - 1; i >= 0; i--) {
    DEBUG_PR<prepro_DEBUG>("Processing wt transitionWts[", i,
                           "] = ", transitionWts[i], '\n');
    int n{0};
    for (size_t c = 0; c < nSofts(); c++) {
      // this is potentially quadratic. But in practice transitionWts is
      // quite small.
      if (transitionWts[i] <= soft_clswts[c] && soft_clswts[c] < wt_added) {
        n++;
        sat_solver->addClause(getSoft(c));
        if (sat_solver->theory_is_unsat()) break;
      }
    }
    wt_added = transitionWts[i];
    DEBUG_PR<prepro_DEBUG>("Added ", n, " soft clauses to sat solver\n");
    if (sat_solver->theory_is_unsat()) break;
    auto canHarden = sat_solver->solvePropBudget(1024 * 1024);
    DEBUG_PR<prepro_DEBUG>("sat solver returns ", canHarden, '\n');
    if (canHarden == l_True)
      wt_hardened = wt_added;
    else
      break;
  }
  DEBUG_PR<prepro_DEBUG>("wt_hardened = ", wt_hardened, '\n');

  if (wt_hardened > wt_max) {
    log(1, "c WCNF hardened 0 soft clauses\n");
  } else {
    MaxHS::Packed_vecs<Lit> tmp;
    vector<Weight> tmpWts;
    int nHardened{0};
    for (size_t i = 0; i < nSofts(); i++)
      if (soft_clswts[i] >= wt_hardened) {
        nHardened++;
        auto sftcls = getSoft(i);
        _addHardClause(sftcls);
      } else {
        tmp.addVec(getSoft(i));
        tmpWts.push_back(soft_clswts[i]);
      }
    soft_cls = std::move(tmp);
    soft_clswts = std::move(tmpWts);
    total_cls_wt = 0;
    for (auto wt : soft_clswts) total_cls_wt += wt;
    log(1, "c WCNF hardened ", nHardened,
        " soft clauses. New total_cls wt = ", wt_fmt(total_cls_wt), '\n');
  }
}

void Wcnf::subEqsAndUnits() {
  if (unsat || !nHards()) return;
  // use cadical sat solver to find units and equalities. Use these to
  // reduce the formula.

  SatSolver_uniqp sat_solver{new cadicalSolver};
  feed_sat_solver_hards(sat_solver.get());
  if (sat_solver->no_elim_simplify() == l_False) {
    unsat = true;
    return;
  }

  // 1) replace hards with hards as reduced by sat solver
  // auto hard_units = sat_solver->getForced(0);
  // log(1, "c WCNF found ", hard_units.size(),
  //     " hard units after simplification\n");
  sat_solver->updateFixed();
  log(1, "c WCNF sat_solver found ", sat_solver->getNFixed(),
      " hard units after simplification\n");

  struct {
    size_t ncls, total_size;
  } before, after;
  before = {hard_cls.size(), hard_cls.total_size()};
  // hard_cls = reduce_by_eqs_and_units(sat_solver.get(), hard_cls, false);
  auto [new_hard_cls, new_units] = sat_solver->get_all_clauses();
  for (auto e : new_units) sat_solver->addClause(e);
  hard_cls = std::move(new_hard_cls);
  after = {hard_cls.size(), hard_cls.total_size()};

  log(1, "c WCNF: Units and equality reduced #hard clauses by ",
      before.ncls - after.ncls, ", total lits in hard clauses by ",
      before.total_size - after.total_size, " and increased units by ",
      new_units.size(), "\n");

  // 2. reduce softs by found equalities and units.
  struct {
    size_t ncls, total_size;
    Weight base_cost;
  } before_s, after_s;
  nfalse_softs = 0;
  before_s = {soft_cls.size(), soft_cls.total_size(), base_cost};
  soft_cls = reduce_by_eqs_and_units(sat_solver.get(), soft_cls, true);
  after_s = {soft_cls.size(), soft_cls.total_size(), base_cost};
  log(1, "c WCNF: Units and equality reduced #soft clauses by ",
      before_s.ncls - after_s.ncls, ", total lits in soft clauses by ",
      before_s.total_size - after_s.total_size, ", and increased base cost by ",
      after_s.base_cost - before_s.base_cost, " (base_cost = ", base_cost,
      ") falsifying ", nfalse_softs, " softs\n");

  // 3. Remove subsumed softs (hards should have already been reduced by
  //    cadical.
  before_s = {soft_cls.size(), soft_cls.total_size(), base_cost};
  rm_subsumed_softs();
  after_s = {soft_cls.size(), soft_cls.total_size(), base_cost};
  log(1, "c WCNF: rm_subsumed_softs reduced #soft clauses by",
      before_s.ncls - after_s.ncls, ", total lits in soft clauses by ",
      before_s.total_size - after_s.total_size, ", and increased base cost by ",
      after_s.base_cost - before_s.base_cost, " (base_cost = ", base_cost,
      ")\n");

  // 4. Reduction by pures.
  int pures{0};
  vector<uint8_t> appears(nVars(), 0);
  for (const auto cls : hard_cls)
    for (auto l : cls) appears[var(l)] |= sign(l) ? 0b01 : 0b10;
  for (auto cls : soft_cls)
    for (auto l : cls) appears[var(l)] |= sign(l) ? 0b01 : 0b10;
  for (size_t v = 0; v != nVars(); v++) {
    if (appears[v] == 0b01) {  // pure -v
      ++pures;
      sat_solver->addClause(mkLit(v, true));
    } else if (appears[v] == 0b10) {  // pure v
      ++pures;
      sat_solver->addClause(mkLit(v, false));
    }
  }

  cout << "c WCNF: found " << pures << " pure literals\n";
  if (pures) {
    sat_solver->unit_propagate();
    hard_cls = reduce_by_eqs_and_units(sat_solver.get(), hard_cls, false);
    soft_cls = reduce_by_eqs_and_units(sat_solver.get(), soft_cls, true);
  }

  // 5. Set up eqLit and fixed_tvals maps---so solutions to original
  //    formula can be reconstructed (sat_solver will go out of scope)
  for (size_t i = 0; i < nVars(); i++) {
    auto lt = mkLit(i);
    set_eqLit(lt, sat_solver->eqLit(lt));
    if constexpr (prepro_DEBUG) {
      if (sat_solver->eqLit(lt) != ~sat_solver->eqLit(~lt))
        log(0,
            "c ERROR: sat_solver eq map not setting opposite literals "
            "properly\n");
    }
  }
  for (auto lt : sat_solver->getForced(0)) {
    if constexpr (prepro_DEBUG) {
      if (fixed_value(get_eqLit(lt)) == l_False)
        log(0, "c ERROR: eqLit root of unit already set to false\n");
    }
    set_fixed(get_eqLit(lt));
    if constexpr (prepro_DEBUG) {
      if (fixed_value(get_eqLit(~lt)) != l_False)
        log(0, "c ERROR: fixed value not properly set\n");
    }
  }
  for (size_t i = 0; i < nVars(); i++) {
    auto lt = mkLit(i);
    if constexpr (prepro_DEBUG) {
      if (get_eqLit(lt) != lt && fixed_value(lt) != l_Undef)
        log(0, "c error: non-eq root set in fixed value before it should be\n");
    }
    set_fixed_equal(lt, get_eqLit(lt));
    if constexpr (prepro_DEBUG) {
      if (fixed_value(lt) != fixed_value(get_eqLit(lt)) ||
          fixed_value(~lt) != fixed_value(~get_eqLit(lt)) ||
          fixed_value(~lt) != fixed_value(get_eqLit(~lt)))
        log(0, "c ERROR: set_fixed_equal has an error\n");
    }
  }

  if constexpr (prepro_DEBUG) {
    // check eqLit and fixed_tvals map
    for (size_t i = 0; i < nVars(); i++) {
      auto lt = mkLit(i);
      auto eq = sat_solver->eqLit(lt);
      if (eq != ~sat_solver->eqLit(~lt))
        log(0, "c ERROR: sat_solver->eqLit(", lt, ") = ", sat_solver->eqLit(lt),
            "but sat_solver->eqLit(", ~lt, ") = ", sat_solver->eqLit(~lt),
            "\n");
      if (sat_solver->eqLit(eq) != eq)
        log(0, "c ERROR: sat_solver->eqLit(", lt, ") = ", sat_solver->eqLit(lt),
            "but sat_solver->eqLit(sat_solver->eqLIT(", lt,
            ") = ", sat_solver->eqLit(sat_solver->eqLit(lt)),
            " (not collapsed) \n");
    }
    for (size_t i = 0; i < nVars(); i++) {
      auto lt = mkLit(i);
      auto eq = get_eqLit(lt);
      if (~eq != get_eqLit(~lt))
        log(0, "c ERROR: get_eqLit(", lt, ") = ", get_eqLit(lt),
            "but get_eqLit(", ~lt, ") = ", get_eqLit(~lt), "\n");
      if (eq != get_eqLit(eq))
        log(0, "c ERROR: get_eqLit(", lt, ") = ", get_eqLit(lt),
            "but get_eqLit(get_eqLit(", lt, ") = ", get_eqLit(get_eqLit(lt)),
            " (not collapsed) \n");
      if (fixed_value(lt) != fixed_value(eq) ||
          fixed_value(~lt) != fixed_value(~eq))
        log(0, "c  ERROR: fixed_valuet(", lt, ") = ", fixed_value(lt),
            " fixed_value(get_eqLit(", lt, ")) = ", fixed_value(get_eqLit(lt)),
            " fixed_value(", ~lt, ") = ", fixed_value(~lt),
            " fixed_value(get_eqLit(", ~lt,
            ")) = ", fixed_value(get_eqLit(~lt)), "\n");
    }
  }
}

MaxHS::Packed_vecs<Lit> Wcnf::reduce_by_eqs_and_units(
    SatSolver* sat_solver, const MaxHS::Packed_vecs<Lit>& clauses, bool softs) {
  if (unsat) return MaxHS::Packed_vecs<Lit>{};
  MaxHS::Packed_vecs<Lit> tmp;
  size_t j = 0;
  for (size_t i = 0; i < clauses.size(); i++) {
    vector<Lit> cls{clauses.getVec(i)};
    if (sat_solver->reduce_clause_by_eqs_and_fixed_values(cls) &&
        rm_dups_and_check_tautology(cls)) {
      if (cls.empty()) {
        if (softs) {
          base_cost += soft_clswts[i];
          ++nfalse_softs;
        } else {
          unsat = true;
          return MaxHS::Packed_vecs<Lit>{};
        }
      } else {
        tmp.addVec(cls);
        if (softs) soft_clswts[j++] = soft_clswts[i];
      }
    }
  }

  if (softs) {
    soft_clswts.resize(j);
    soft_clswts.shrink_to_fit();
    total_cls_wt = 0;
    for (auto wt : soft_clswts) total_cls_wt += wt;
  }
  return tmp;
}

void Wcnf::rm_subsumed_softs() {
  if (unsat) return;
  // hards have been subsumption reduced by sat solver But softs need
  // special treatment as they can only be subsumed by a hard

  // A soft can also be shortened by self-subsuming resolution: let
  // soft s = (A, x) and hard h = (A', ~x) where A' \subseq A we can
  // replace s with the soft (A) (having the same weight) If ~x,
  // then s will yield a cost iff (A) If x, then A' is a hard and
  // and (A) is always satisfied so there is no cost (as in the
  // original).

  // Identical softs can have their weights merged. Unit softs that
  // conflict can be resolved adding the minimal wt to the basecost and
  // the residue to the polarity that had higher weight.

  // 1. remove duplicate softs (and do unit resolutions among the
  // softs)
  vector<SftClsData> sft_cdata;
  bool softs_are_all_units = true;
  int n_res{}, n_rm{}, nsubsumed{}, nself_subsumed{};

  initSftClsData(sft_cdata);
  auto byHash = [](const SftClsData& a, const SftClsData& b) {
    return a.hash < b.hash;
  };
  sort(sft_cdata.begin(), sft_cdata.end(), byHash);
  for (size_t i = 0; i < sft_cdata.size(); i++) {
    if (sft_cdata[i].w == 0) continue;  // w==0 indicates clause is deleted;
    auto vi = soft_cls[sft_cdata[i].index];
    if (vi.size() > 1) softs_are_all_units = false;
    for (size_t j = i + 1;
         j < sft_cdata.size() && sft_cdata[j].hash == sft_cdata[i].hash; j++) {
      if (sft_cdata[j].w == 0) continue;
      auto vj = soft_cls[sft_cdata[j].index];
      if (vi.size() == 1 && vj.size() == 1 && vi[0] == ~vj[0]) {
        ++n_res, ++n_rm;
        Weight cost{}, residue{};
        if (sft_cdata[i].w < sft_cdata[j].w) {
          vi[0] = vj[0];  // higher cost unit is preserved
          cost = sft_cdata[i].w;
          residue = sft_cdata[j].w - cost;
        } else {
          cost = sft_cdata[j].w;
          residue = sft_cdata[i].w - cost;
        }
        // note that if costs are equal residue becomes 0 and
        // both clauses vanish
        base_cost += cost;
        ++nfalse_softs;
        sft_cdata[i].w = residue;
        sft_cdata[j].w = 0;
        if (residue == 0) ++n_rm;
      } else if (eqSoftCls(sft_cdata[i], sft_cdata[j])) {
        // equal clauses are merged.
        sft_cdata[i].w += sft_cdata[j].w;  // join softs.
        sft_cdata[j].w = 0;
        ++n_rm;
      }
    }
  }
  auto byIndex = [](const SftClsData& a, const SftClsData& b) {
    return a.index < b.index;
  };
  sort(sft_cdata.begin(), sft_cdata.end(), byIndex);
  log(1, "c WCNF: removed ", n_rm, " softs by resolution and found ", n_res,
      " unit soft resolutions\n");

  if (!hard_cls.empty() && !softs_are_all_units) {
    OneWatch occur{hard_cls};
    vector<uint8_t> marks(nVars() * 2, false);
    // subroutines
    auto chk_watch = [&](const Watch& w) -> std::pair<bool, Lit> {
      // check hard clause in watch to see if it subsumes or self-subsumes
      // soft clause s = soft_cls[sc]. Note that w is on watch list of s
      // so the first literal of the hard clause it points is contained in s.
      Lit sub_idx = lit_Undef;
      if (!marks[toInt(w.blocker)] && !marks[toInt(~w.blocker)])
        return {false, lit_Undef};
      for (size_t i{}; i < hard_cls[w.cls_idx].size(); ++i) {
        Lit x = hard_cls[w.cls_idx][i];
        if (marks[toInt(x)])
          continue;
        else if (marks[toInt(~x)]) {
          if (sub_idx != lit_Undef)
            // contains more than one negated lit of hard
            return {false, lit_Undef};
          else
            sub_idx = x;
        } else
          return {false, lit_Undef};
      }
      return {true, sub_idx};
    };

    auto chk_lit = [&](size_t sz, Lit x) -> std::pair<bool, Lit> {
      // return {true, NOIDX} if s = soft_cls[sc] is subsumed by some
      // clause on x's watch list. return {true, i>=0} if we can
      // remove the i'th literal of s via self-subsumption from some
      // lit on x's watch list.
      Lit sub_idx = lit_Undef;
      for (auto& w : occur[x]) {
        // watch list sorted by clause size, no need to check hard clauses
        // that are too large to be subsets
        if (w.cls_size > sz) return {false, lit_Undef};
        auto [subsumed, self_subsumed] = chk_watch(w);
        if (subsumed && self_subsumed == lit_Undef)
          return {subsumed, self_subsumed};
        if (sub_idx == lit_Undef && subsumed && self_subsumed != lit_Undef)
          sub_idx = self_subsumed;
      }
      return {sub_idx != lit_Undef, sub_idx};
    };

    auto chk_subsumed = [&](size_t sc) -> std::pair<bool, Lit> {
      // Return {true, NOIDX} if clause s = soft_cls[sc] is
      // subsumed. Return {true, i>=0} if we can remove the i'th
      // literal of s via self-subsumption.  To determine this we look
      // at the watch lists of each of s's literals. If a subset of s
      // is a hard clause h, then h will appear on one of these watch
      // lists. Checking for self-subsumption is more work---we also
      // have to check the watch lists of -l for l in s. Do complete
      // check for subsumption, leaving self-subsumption as backup.
      auto cls = soft_cls[sc];
      Lit sub_idx = lit_Undef;
      for (auto x : cls) {
        auto [subsumed, self_subsumed] = chk_lit(cls.size(), x);
        if (subsumed && self_subsumed == lit_Undef) {
          return {subsumed, self_subsumed};
        }
        if (sub_idx == lit_Undef && subsumed && self_subsumed != lit_Undef)
          sub_idx = self_subsumed;
      }
      if (sub_idx != lit_Undef) return {true, sub_idx};
      for (auto x : cls) {
        auto [subsumed, self_subsumed] = chk_lit(cls.size(), ~x);
        if (subsumed && self_subsumed != lit_Undef)
          return {subsumed, self_subsumed};
      }
      return {false, lit_Undef};
    };

    // main loop
    for (size_t i{0}; i < sft_cdata.size(); ++i) {
      if (sft_cdata[i].w == 0) continue;
      auto sc = sft_cdata[i].index;
      for (auto x : soft_cls[sc]) marks[toInt(x)] = true;
      auto [subsumed, self_subsumed] = chk_subsumed(sc);
      if (subsumed && self_subsumed == lit_Undef) {
        ++nsubsumed;
        sft_cdata[i].w = 0;
      } else if (subsumed && self_subsumed != lit_Undef) {
        ++nself_subsumed;
        if (std::remove(soft_cls[sc].begin(), soft_cls[sc].end(),
                        ~self_subsumed) != soft_cls[sc].end() - 1) {
          log(1, "c WCNF: ERROR self subsuming did not remove one literal\n");
          cout << " self_subsumed = " << self_subsumed << "\nsoft_cls[" << sc
               << "] = [";
          for (auto z : soft_cls[sc]) cout << z << " ";
          cout << "]\n";
        }
        soft_cls[sc].shrink(1);
      }
      for (auto x : soft_cls[sc]) marks[toInt(x)] = false;
    }

    log(1, "c WCNF: subsumed ", nsubsumed, " soft clauses with the hards\n");
    log(1, "c WCNF: self_subsumed ", nself_subsumed,
        " soft clauses with the hards\n");
  }

  // rewrite the softs
  if (nsubsumed + nself_subsumed + n_rm > 0) {
    MaxHS::Packed_vecs<Lit> tmpS;
    vector<Weight> tmp_wts;
    vector<Lit> new_soft;

    for (size_t i{0}; i < sft_cdata.size(); ++i) {
      if (sft_cdata[i].w > 0) {
        new_soft = soft_cls.getVec(sft_cdata[i].index);
        if (new_soft.empty()) {
          base_cost += sft_cdata[i].w;
          ++nfalse_softs;
        } else {
          tmpS.addVec(new_soft);
          tmp_wts.push_back(sft_cdata[i].w);
        }
      }
    }
    soft_cls = std::move(tmpS);
    soft_clswts = std::move(tmp_wts);
    total_cls_wt = 0;
    for (auto wt : soft_clswts) total_cls_wt += wt;
  }
}

bool Wcnf::eqSoftCls(const SftClsData& a, const SftClsData& b) {
  // relies on all clauses being sorted on input.
  if (a.index == b.index) return true;
  if (soft_cls[a.index].size() != soft_cls[b.index].size()) return false;
  for (size_t x = 0; x < soft_cls[a.index].size(); x++)
    if (soft_cls[a.index][x] != soft_cls[b.index][x]) return false;
  return true;
}

void Wcnf::initSftClsData(vector<SftClsData>& cdata) {
  // set up hash codes for all soft clauses to check for duplicate
  // softs. For unit softs use the same hash code for -x and x (to
  // support resolutions.
  for (size_t i = 0; i < nSofts(); i++)
    if (softSize(i) == 1) {
      vector<Var> unit{var(soft_cls[i][0])};
      cdata.emplace_back(static_cast<uint32_t>(i),
                         hashCode(unit.begin(), unit.end()), getWt(i));
    } else
      cdata.emplace_back(static_cast<uint32_t>(i),
                         hashCode(soft_cls[i].begin(), soft_cls[i].end()),
                         getWt(i));

  if (params.verbosity > 2) {
    vector<uint32_t> hcodes;
    for (auto& item : cdata) hcodes.push_back(item.hash);
    sort(hcodes.begin(), hcodes.end());
    auto it = std::unique(hcodes.begin(), hcodes.end());
    hcodes.resize(std::distance(hcodes.begin(), it));
    cout << "c Hashed " << cdata.size() << " clauses with "
         << it - hcodes.begin() << " (" << hcodes.size()
         << ") different hash codes\n";
  }
}

void Wcnf::feed_sat_solver_hards(SatSolver* const slv) {
  for (size_t i = 0; i < nHards(); i++) slv->addClause(getHard(i));
  if (slv->theory_is_unsat()) {
    unsat = true;
    log(1, "c WCNF found hards to be unsat\n");
  }
}

void Wcnf::mxBvars() {
  // modify the WCNF by finding at most one constraints among the bvars
  // and replacing all of these vars by one bvar.
  Bvars bvars{this};
  if (params.verbosity > 4) {
    cout << "BEFORE mxbvars\n";
    printFormula(bvars);
  }
  processMxs(mxFinder(bvars), bvars);
  if (params.verbosity > 4) {
    Bvars newbvars{this};
    cout << "AFTER mxbvars\n";
    printFormula(newbvars);
  }
}

void Wcnf::processMxs(const vector<vector<Lit>> mxs, Bvars& bvars) {
  // mxs should be disjoint collection of mx sets. Each set should be a
  // non empty set of blits all of which have the same weight. These
  // blits have the property that at most one of them can be true
  //(given the hard clauses).

  // If the blits are cores (making them true relaxes the soft clause),
  // then at most one of the corresponding soft clauses can be
  // falsified. If the blits are non-cores then at most one of the
  // corresponding soft clauses can be true.

  // Modify the WCNF to account for this mx.
  // Also store the mxes so that users of Wcnf can access them.

  if (unsat) return;

  // debug
  // cout << "Before processmx\n";
  // printFormula(bvars);

  vector<char> delMarks(nSofts(), 0);
  auto orig_nsofts = nSofts();
  vector<Lit> blits{};
  for (auto& mx : mxs) {
    // debug
    /*cout << "Processing " << (bvars.isCore(mx[0]) ? "Core Mx: [ " :
    "NonCore Mx: [ "); for(auto l : mx) cout << l << "(ci = " <<
    bvars.clsIndex(l) <<
    "),
    "; cout << "]\n"; if(mx.size() < 2) cout << "ERROR: Got mutual exclusion
    with less than 2 vars" << mx << "\n"; bool dc = bvars.isCore(mx[0]);
    for(auto l : mx) {
      if(dc && !bvars.isCore(l) || !dc && bvars.isCore(l))
        cout << "ERROR mx with mixed core/non-cores" << mx << "\n";
      if(bvars.wt(var(l)) != bvars.wt(var(mx[0])))
      cout << "ERROR mx with different wts" << mx << "\n";
      }*/
    // debug
    if (mx.empty()) {
      cout << "c WARNING. Mx finder returned empty mx\n";
      continue;
    }
    Weight unitWt = mx.size() > 0 ? bvars.wt(var(mx[0])) : 0.0;
    bool core = mx.size() > 0 ? bvars.isCore(mx[0]) : false;
    auto dvar = var_Undef;
    auto dlit = lit_Undef;
    blits.clear();
    // if the mx is a core-mx it means that of a set of soft clauses
    // at most one can be false. One could define a new variable d
    // that is true if one of these soft clauses if false, and then a
    // unit clause (-d) asking the solver to satisfy all of the soft
    // clauses.  However, having the solver assume -d leads to weaker
    // conflicts than having it assume (-b1, -b2, ..., -bk) (where
    // these are the soft clauses in the mutex. So we do not do the
    // transformation, but we do obtain a bvar for each soft and tell
    // the solver that these bvars are mutex (this info can be added
    // to the hitting set solvers.
    if (core) {
      for (auto l : mx) {
        auto ci = bvars.clsIndex(l);
        auto sftcls = getSoft(ci);
        if (sftcls.size() == 0) {
          cout << "c ERROR WCNF processMxs encountered zero length soft "
                  "clause\n";
          continue;
        } else if (sftcls.size() == 1)
          blits.push_back(~sftcls[0]);
        else {
          auto blit = mkLit(bvars.newBVar());
          blits.push_back(blit);
          sftcls.push_back(blit);
          delMarks[ci] = 1;
          _addHardClause(sftcls);
          _addSoftClause(~blit, unitWt);
        }
      }
      mutexes.emplace_back(blits, core, lit_Undef);
      // dvar = bvars.newBvar();
      // dlit = mkLit(dvar);
      // for (auto l : mx) {
      //   auto ci = bvars.clsIndex(l);
      //   auto sftcls = getSoft(ci);
      //   if (sftcls.size() == 0) {
      //     cout << "c ERROR WCNF processMxs encountered zero length soft "
      //             "clause\n";
      //     continue;
      //   }
      //   sftcls.push_back(dlit);
      //   delMarks[ci] = 1;
      //   _addHardClause(sftcls);
      // }
      // _addSoftClause(~dlit, unitWt);
    } else {
      // If the mutex is a non-core then there is a set of soft clauses at
      // most one of which can be TRUE. So we can add the wt of all but one
      // to base cost. In additon the transformation to a single new soft
      // clause
      // (-d) where -d == one is TRUE, (d == all are false) is beneficial.
      // As the solver assuming -d corresponds to assuming a disjunction
      // (-b1 \/ -b2 \/ ... \/ -bk) (where the bi are the softs in the
      // mutex. Assuming this disjunction leads to more powerful clauses
      // than picking an arbitrary -bi to make true (i.e., an arbitrary soft
      // clause to make true).
      //
      // We encode d as -d -> the disjunction of all literals in these
      // soft clauses. But TODO we should also check the approach of Bacchus
      // and Katesirelos CAV 2015 of introducing a b-var for each soft
      // clause and then have -d -> the disjunction of these bvars

      for (auto l : mx) {
        auto ci = bvars.clsIndex(l);
        auto sftcls = getSoft(ci);
        if (sftcls.size() == 0) {
          cout << "c ERROR WCNF processMxs encountered zero length soft "
                  "clause\n";
          continue;
        }
        for (auto x : sftcls) blits.push_back(x);  // union of softs
        delMarks[ci] = 1;
      }
      dvar = bvars.newBVar();
      dlit = mkLit(dvar);
      blits.push_back(dlit);
      _addHardClause(blits);
      base_cost += unitWt * (mx.size() - 1);
      nfalse_softs += (mx.size() - 1);
      _addSoftClause(~dlit, unitWt);

      // if (!orig_usoft.empty())
      //   mutexes.emplace_back(std::move(orig_usoft), coreMx, lit_Undef);
    }
  }

  // now rewite the softs
  MaxHS::Packed_vecs<Lit> tmp;
  size_t j{0};
  for (size_t i = 0; i < nSofts(); i++)
    if (i >= delMarks.size() || !delMarks[i]) {
      // delmarks don't extend to newly added softs
      tmp.addVec(getSoft(i));
      soft_clswts[j++] = soft_clswts[i];
    }
  soft_clswts.resize(j);
  soft_clswts.shrink_to_fit();
  soft_cls = std::move(tmp);

  computeWtInfo();

  if (params.verbosity) {
    cout << "c WCNF mutexes: original #softs " << orig_nsofts
         << " #softs after mx-transforms " << nSofts() << "\n";
    cout << "c WCNF mutexes: reduction in softs " << orig_nsofts - nSofts()
         << "\n";
  }

  if (params.verbosity > 2) {
    cout << "Process mx\n";
    cout << "mutexes\n";
    for (auto& mx : mutexes) cout << mx << "\n";
  }

  // debug
  // Bvars newbvars {this};
  // cout << "After processmx\n";
  // printFormula(newbvars);
}

class MXFinder {
  // helper class for finding mutually exclusive bvars
 public:
  MXFinder(Wcnf* f, Bvars& b)
      : bvars{b},
        theWcnf{f},
        sat_solver{new cadicalSolver},
        blitMarks(2 * bvars.n_vars(), 0) {}
  ~MXFinder() {
    for (auto p : blitMXes)
      if (p) delete p;
  }
  bool findMxs(vector<vector<Lit>>&);
  int nImpCalls{0};

 private:
  Bvars bvars;
  Wcnf* theWcnf;
  SatSolver_uniqp sat_solver;
  vector<uint8_t> blitMarks;
  const uint8_t inmx{1};
  const uint8_t in2s{2};
  bool fbeq();
  size_t getMXLitSize(Lit);
  const vector<Lit>* getMXLits(Lit);
  uint64_t totalMxMem{};
  vector<vector<Lit>*> blitMXes;
  void getMXRecomputeSizes(const vector<Lit>&);
  vector<Lit> growMx(Lit start);
  // debug
  void MXprintLit(Lit l) {
    cout << l << " (mkr=" << (int)blitMarks[toInt(l)]
         << (bvars.isCore(l) ? " C " : " NC ") << " wt = " << bvars.wt(var(l))
         << ") ";
  }
};

vector<vector<Lit>> Wcnf::mxFinder(Bvars& bvars) {
  // Return collection of mutally exclusive bvars. Where each set of
  // Bvars are cores or are non-cores and each is associated with a
  // soft clause of identical weight; finder does all the work: we use
  // a function wrapper to ensure MXFinder is deallocated after it has
  // done its work.
  MXFinder finder{this, bvars};
  vector<vector<Lit>> mxs;
  if (!finder.findMxs(mxs)) unsat = true;
  if (params.verbosity > 0)
    cout << "c WCNF mx finder used " << finder.nImpCalls
         << " calls to UP engine\n";
  return mxs;
}

bool MXFinder::findMxs(vector<vector<Lit>>& mxs) {
  // Top level findMxs method of MXFinder class.  find mx collections
  // of bvars and store in mxs. Return false if we found a
  // contradiction.

  // TODO: currently no attempt is made to reclaim
  // memory during the computation.  In particular, once a blit is in
  // an MX we shouldn't have to store its set of implicants. (BUT the
  // marking code is tricky so this is not necessarily an easy mod,
  // this needs a redesign..

  bool timedOut{false};
  double start_time = Minisat::cpuTime();

  // 1. Initialize solver with fbeq
  if (!fbeq()) {
    if (params.verbosity > 0)
      cout << "c WCNF detected input to be unsat during preprocessing\n";
    return false;
  }

  // Two stage. Note that absorbing a blit into a mutex blocks it and
  // its negation from being in any other mutex. So try to grow big
  // mutexes we delay the processing of mutexes of size 2. E.g., it
  // could be that (b1,b2) are mutex but so is (b2, b3, b4, b5) so we
  // don't want to absorbe b2 into (b1, b2).

  vector<Lit> toProcess;
  vector<Lit> twos;  // blits that might generate 2s. Process later.

  // toProcess will be used as a stack. Non core mx bump the base cost
  // so process those first (put at back of vector--top of stack)
  bool find_cores = params.mx_find_mxes == 3 || params.mx_find_mxes == 1;
  bool find_ncores = params.mx_find_mxes == 3 || params.mx_find_mxes == 2;
  if (find_cores)
    for (size_t i = 0; i < theWcnf->nSofts(); ++i) {
      toProcess.push_back(bvars.litOfCls(i));
    }
  if (find_ncores)
    for (size_t i = 0; i < theWcnf->nSofts(); ++i) {
      toProcess.push_back(~bvars.litOfCls(i));
    }

  int loops{0};
  while (!toProcess.empty()) {
    // in this loop toProcess.back() is either selected, or blits that
    // are mx with it are processed. Thus the size of
    // toProcess.back()'s mx lit set is decreased. So eventually toProcess
    // will become empty.
    ++loops;
    // check for time out and mem-limit.
    if (totalMxMem >=
            static_cast<uint64_t>(1024 * 1024) * params.mx_mem_limit ||
        (params.mx_cpu_lim > 0 && loops % 500 == 0 &&
         (Minisat::cpuTime() - start_time) > params.mx_cpu_lim)) {
      timedOut = true;
      if (totalMxMem >=
          static_cast<uint64_t>(1024 * 1024) * params.mx_mem_limit)
        cout << "c WCNF mx finder hit its memory limit. "
             << "Potentially more mxes could be found with -mx-mem-lim made "
                "larger\n";
      if ((Minisat::cpuTime() - start_time) > params.mx_cpu_lim)
        cout << "c WCNF mx finder hit its time limit. "
             << "Potentially more mxes could be found with -mx-cpu-lim made "
                "larger\n";
      break;
    }

    Lit blit = toProcess.back();
    if (blitMarks[toInt(blit)]) {
      // in an mx or in twos
      toProcess.pop_back();
      continue;
    }
    auto mx = getMXLits(blit);

    // Debug
    /*cout << "findmxs first loop unmarked processing ";
    MXprintLit(blit);
    cout << " mx = [";
    for(auto p : *mx) {
    MXprintLit(p);
    cout << " (" << getMXLitSize(p) << "), ";
    }
    cout << "]\n";*/

    if (mx->size() <= 1) {
      if (mx->size() == 1) {
        blitMarks[toInt(blit)] = in2s;
        twos.push_back(blit);
        // cout << "Pushed into twos\n";
      }
      // cout << "Got rid of\n";
      toProcess.pop_back();
      continue;
    }

    // Potential mx of size > 2 (but not guaranteed!)
    // cout << "Growing from start ";

    Lit start = blit;
    size_t size{mx->size()};

    // Find start.
    for (auto l : *mx) {
      size_t sz = getMXLitSize(l);
      if (sz > size) {
        size = sz;
        start = l;
      }
    }

    // cout << start << "\n";
    // grow from start
    auto tmp = growMx(start);

    // tmp might be small and it might not contain blit
    if (tmp.size() <= 2) {
      // cout << "too small\n";

      // no easy way to get rid of start from toProcess
      // so mark as being in twos
      blitMarks[toInt(blit)] = in2s;
      if (tmp.size() == 2) {
        // but only put there if it has potential
        twos.push_back(start);
        // cout << "put in twos\n";
      }
    }

    else {  // legit mx for this stage
      for (auto b : tmp) {
        // cout << "adding mx\n";
        blitMarks[toInt(b)] = inmx;
        blitMarks[toInt(~b)] = inmx;
      }
      mxs.push_back(std::move(tmp));
    }
  }

  if (!timedOut)
    while (!twos.empty()) {
      Lit blit = twos.back();
      // cout << "Processing twos " << blit << "\n";
      twos.pop_back();
      if (blitMarks[toInt(blit)] == inmx) {
        // cout << "Already marked\n";
        continue;
      }
      auto tmp = growMx(blit);
      if (tmp.size() > 1) {
        if (tmp.size() > 2) cout << "c WARNING. WCNF large mx got into twos\n";
        for (auto b : tmp) {
          blitMarks[toInt(b)] = inmx;
          blitMarks[toInt(~b)] = inmx;
        }
        mxs.push_back(std::move(tmp));
      }
    }

  if (params.verbosity > 0) {
    cout << "c WCNF mutexes: #mutexes found = " << mxs.size() << "\n";
    if (mxs.size() > 0) {
      double core_total_size{0};
      double ncore_total_size{0};
      int cores{0};
      int ncores{0};
      for (auto& mx : mxs) {
        if (bvars.isCore(mx[0])) {
          cores++;
          core_total_size += mx.size();
        } else {
          ncores++;
          ncore_total_size += mx.size();
        }
      }
      cout << "c WCNF mutexes: #cores mutexes = " << cores;
      if (cores) cout << " ave. size = " << fix4_fmt(core_total_size / cores);
      cout << "\n";
      cout << "c WCNF mutexes: #non-cores mutexes = " << ncores;
      if (ncores)
        cout << " ave. size = " << fix4_fmt(ncore_total_size / ncores);
      cout << "\n";
      cout << "c WCNF mutexes: time used = "
           << time_fmt(Minisat::cpuTime() - start_time) << "\n";
    }
  }
  return true;
}

vector<Lit> MXFinder::growMx(Lit start) {
  // Starting with blit start, grow a mx constraint (at most one)
  // The algorithm is as follows.
  // We have a set mx that is initially start, and a set of candidates
  // initially the negations of the appropriate implications of start
  //(returned by getMXLits)
  // Appropriate means that these negations are cores if start is a core
  // and they are non-cores if start is a non-core. Furthermore, none
  // of them is marked (as used in another mx constraint), and they all
  // correspond to soft clauses with identical weight as start's soft
  // clause. These restrictions are maintained by the subroutine getMXLits

  // The algorithm's invariant is that
  //(a) the lits in mx form an mx constraint
  //(b) every single lit in candidates is mx with each lit in mx
  //    Thus mx can be grown by any member of candidates.

  // 1. All l in candidates are MX (mutually exclusive) with all c in
  //   mx. So when we move l from candidates to mx (i.e., we select l)
  //   we must remove from candiates all l' that are not MX with l.
  //   Here we sort candidates by |mxset(l)\cap candidates|. So that
  //   if l is selected we can retain as many members of candidates as
  //   is possible. We only sort at the start. A greedy method would
  //   always recompute the size of this intersection since candidates
  //   is shrinking after every l is selected. But we don't do that
  //   to avoid this cost.

  vector<Lit> origCandidates(*getMXLits(start));  // need copy
  std::set<Lit> candidates(origCandidates.begin(), origCandidates.end());

  // compute |mxset(l) \cap candidates| or all l\in candidates
  vector<int> interSize;
  for (auto l : origCandidates) {
    int count{0};
    int ci{bvars.clsIndex(l)};
    for (auto l1 : *getMXLits(l))
      if (candidates.find(l1) != candidates.end()) ++count;
    if (static_cast<size_t>(ci) >= interSize.size())
      interSize.resize(ci + 1, 0);
    interSize[ci] = count;
  }

  auto mxSize = [this, &interSize](Lit l1, Lit l2) {
    return interSize[bvars.clsIndex(l1)] > interSize[bvars.clsIndex(l2)];
  };
  sort(origCandidates.begin(), origCandidates.end(), mxSize);

  // debug*****
  /*cout << "growMX on " << start << "\n";
    cout << start << " mx = " << *getMXLits(start) << "\n";
    cout << "candidiates\n";
    for(size_t i = 0; i < origCandidates.size(); i++) {
    auto l = origCandidates[i];
    auto ci = bvars.clsIndex(l);
    cout << l << " mx = "
    << *getMXLits(l) << " [" << interSize[ci] << "]\n";
    }*/
  // debug*****

  vector<Lit> mx{start};
  size_t cani{0};
  while (!candidates.empty()) {
    /*cout << "Loop #" << cani << " candidiates [";
      for(auto l : candidates)
      cout << l << ", ";
      cout << "]\n";*/

    auto l = origCandidates[cani++];  // select according to computed order
    // cout << " growing " << cani << " lit = " << l;
    auto it = candidates.find(l);
    if (it == candidates.end()) {
      // cout << " Not found\n";
      continue;
    }
    // cout << " Found\n";
    mx.push_back(l);
    candidates.erase(it);
    // remove from candidates all elements not mx with the newly selected l
    auto l_mx = getMXLits(l);
    std::set<Lit> newLitMx(l_mx->begin(), l_mx->end());
    for (auto p = candidates.begin(); p != candidates.end();) {
      // cout << "checking active candiate " << *p;
      if (newLitMx.find(*p) == newLitMx.end()) {
        // this candidate not mx with newly added lit---not a candidate
        // anymore; cout << " Pruning\n";
        candidates.erase(p++);  // sets! Erase invalidates p, so pass copy
                                // and increment here.
      } else {
        ++p;
        // cout << " Keeping\n";
      }
    }
  }

  // debug
  // cout << "Grow found mx " << mx << "\n";
  /*for(size_t i=0; i<mx.size(); i++)
    for(size_t j=i+1; j < mx.size(); j++)
    if(!solver.isMX(mx[i], mx[j])) {
        cout << "ERROR MX faulty (" << mx[i] << ", " << mx[j] << ")\n";
        break;
    }
    for(size_t i=0; i< mx.size(); i++)
    if(    (bvars.isCore(mx[0]) && bvars.isNonCore(mx[i]))
        || (bvars.isNonCore(mx[0]) && bvars.isCore(mx[i]))
    || (bvars.wt(var(mx[0])) != bvars.wt(var(mx[0])))   ) {
    cout << "ERROR MX faulty (" << mx[0] << ", " << mx[i] << ")\n";
    break;
    }*/

  return mx;
}

bool MXFinder::fbeq() {
  // Add fbeq to solver.
  for (size_t i = 0; i < theWcnf->nHards(); i++)
    sat_solver->addClause(theWcnf->getHard(i));
  if (sat_solver->theory_is_unsat()) return false;

  for (size_t i = 0; i < theWcnf->nSofts(); i++) {
    Lit blit = bvars.litOfCls(i);
    if (theWcnf->softSize(i) > 1) {
      vector<Lit> sftCls{theWcnf->getSoft(i)};
      sftCls.push_back(blit);
      sat_solver->addClause(sftCls);
      if (sat_solver->theory_is_unsat()) return false;

      vector<Lit> eqcls(2, lit_Undef);
      eqcls[1] = ~blit;
      for (auto l : theWcnf->softs()[i]) {
        eqcls[0] = ~l;
        sat_solver->addClause(eqcls);
        if (sat_solver->theory_is_unsat()) return false;
      }
    }
  }
  return true;
}

const vector<Lit>* MXFinder::getMXLits(Lit l) {
  // return *unmarked* literals of same type (same wt and same core status)
  // that are mutually exclusive with l. This is accomplished by asking the
  // sat solver to find the implications of l that have the same weight and
  // the opposite core status. If l --> l1, then l and ~l1 are mutually
  // exclusive. To make this more efficient we remember the compute sets in
  // blitMXex.
  //
  // Since literals are marked at various stages, we must prune these cached
  // sets to remove newly marked literals (??do we need to do this??)

  Weight lWT{0};  // rather than set up a class of template, we capture this
                  // local var in our lambda.
  auto not_coreMX = [this, &lWT](Lit l1) {
    return blitMarks[toInt(l1)] == inmx || !bvars.isNonCore(l1) ||
        bvars.wt(var(l1)) != lWT;
  };
  auto not_nonCoreMX = [this, &lWT](Lit l1) {
    return blitMarks[toInt(l1)] == inmx || !bvars.isCore(l1) ||
        bvars.wt(var(l1)) != lWT;
  };

  if (blitMXes.size() <= static_cast<size_t>(toInt(l)))
    blitMXes.resize(toInt(l) + 1, nullptr);

  vector<Lit> imps;
  if (!blitMXes[toInt(l)]) {  // not cached compute
    if (totalMxMem >=
        static_cast<uint64_t>(1024 * 1024) * params.mx_mem_limit) {
      // no more space for storing implications (2GB)...pretend that there
      // aren't any
      blitMXes[toInt(l)] = new vector<Lit>(std::move(imps));
      return blitMXes[toInt(l)];
    }

    lWT = bvars.wt(var(l));
    ++nImpCalls;
    sat_solver->findImplications(l, imps);

    if (bvars.isCore(l))
      imps.erase(std::remove_if(imps.begin(), imps.end(), not_coreMX),
                 imps.end());
    else
      imps.erase(std::remove_if(imps.begin(), imps.end(), not_nonCoreMX),
                 imps.end());

    // debug
    /*vector<Lit> t;
      t.push_back(l);
      vector<Lit> o;
      vector<Lit> pruned;
      sat_solver->findImplications(t, o);
      t.clear();
      for(size_t i=0; i < o.size(); ++i) {
      if(blitMarks[toInt(o[i])] != inmx &&
      ((bvars.isCore(l) && bvars.isNonCore(o[i]))
          || (bvars.isNonCore(l) && bvars.isCore(o[i])))
      && (bvars.wt(var(l)) == bvars.wt(var(o[i]))))
      pruned.push_back(o[i]);
      }

      std::sort(pruned.begin(), pruned.end());
      std::sort(imps.begin(), imps.end());
      if(pruned.size() != imps.size()) {
      cout << "ERROR regular interface returned different set of imps for ";
      cout << l << " (marks = " << (int) blitMarks[toInt(l)] <<
      (bvars.isCore(l) ? " Core " : " nonCore "); cout << " wt = " <<
      bvars.wt(var(l)) << ")\n";

      cout << " imps = [ ";
      for(auto x : imps)
      cout << x << " (marks = " << (int) blitMarks[toInt(x)] <<
      (bvars.isCore(x) ? " Core " : " nonCore ")
      << " wt = " << bvars.wt(var(x)) << "), ";
      cout << "] (" << imps.size() << ")\n";

      cout << " pruned = [ ";
      for(auto x : pruned)
      cout << x << " (marks = " << (int) blitMarks[toInt(x)] <<
      (bvars.isCore(x) ? " Core " : " nonCore ")
      << " wt = " << bvars.wt(var(x)) << "), ";
      cout << "] (" << pruned.size() << ")\n";
      }

      else {
      for(size_t i =0; i < pruned.size(); i++)
      if(pruned[i] != imps[i]) {
          cout << "ERROR regular interface returned different set of
      imps\n";

          cout << l << " (marks = " << (int) blitMarks[toInt(l)] <<
      (bvars.isCore(l) ? " Core " : " nonCore "); cout << " wt = " <<
      bvars.wt(var(l)) << ")\n";

          cout << " imps = [ ";
          for(auto x : imps)
      cout << x << " (marks = " << (int) blitMarks[toInt(x)] <<
      (bvars.isCore(x) ? " Core " : " nonCore ")
      << " wt = " << bvars.wt(var(x)) << "), ";
          cout << "] (" << imps.size() << ")\n";

          cout << " pruned = [ ";
          for(auto x : pruned)
      cout << x << " (marks = " << (int) blitMarks[toInt(x)] <<
      (bvars.isCore(x) ? " Core " : " nonCore ")
      << " wt = " << bvars.wt(var(x)) << "), ";
          cout << "] (" << pruned.size() << ")\n";

          break;
      }
      }*/
    // debug

    for (size_t i = 0; i < imps.size(); ++i)
      imps[i] = ~imps[i];  // convert from implication to MX

    totalMxMem += sizeof(Lit) * imps.size();
    blitMXes[toInt(l)] = new vector<Lit>(std::move(imps));
    return blitMXes[toInt(l)];
  }

  // otherwise prune for marks
  auto v = blitMXes[toInt(l)];
  size_t j{0};
  for (size_t i = 0; i < v->size(); ++i)
    if (blitMarks[toInt((*v)[i])] != inmx) (*v)[j++] = (*v)[i];

  v->resize(j);
  return v;
}

size_t MXFinder::getMXLitSize(Lit l) {
  // potentially we could aproximate this by returning
  // the size without pruning for marks and using getMXrecomputeSizes
  // on each found mx to get most of the sizes right.
  return getMXLits(l)->size();
}

void MXFinder::getMXRecomputeSizes(const vector<Lit>& newlyMarked) {
  auto nic = nImpCalls;
  for (auto l : newlyMarked) {  // note we don't need to update vectors of
                                // marked lits
    auto v = getMXLits(l);      // this is automatic from getMXLits. v is set of
                                // lits that need update
    for (auto x : *v) {
      auto vx = blitMXes[toInt(x)];
      if (!vx) continue;
      size_t j{0};
      for (size_t i = 0; i < vx->size(); ++i)
        if (blitMarks[toInt((*vx)[i])] == inmx) (*vx)[j++] = (*vx)[i];
      vx->resize(j);
    }
  }
  if (nImpCalls > nic && params.verbosity > 0)
    cout << "c WARNING getMXRecomputeSizes used some implication calls!\n";
}

// Internal variable number to input file numbering.

void Wcnf::remapVars() {
  vector<char> appears(nVars(), false);
  for (const auto cls : hard_cls)
    for (auto l : cls) {
      appears[var(l)] = true;
      if constexpr (prepro_DEBUG) {
        if (get_eqLit(l) != l)
          log(0, "c ERROR: lit of hards should have been equality replaced:", l,
              " = ", get_eqLit(l), " : ", ~l, " = ", get_eqLit(~l), '\n');
        if (fixed_value(l) != l_Undef)
          log(0, "c ERROR: lit of hard is fixed, should have been removed:", l,
              " = ", fixed_value(l), '\n');
      }
    }

  flipped_vars.resize(nVars(), 0);
  for (const auto cls : soft_cls) {
    for (auto l : cls) {
      appears[var(l)] = true;
      if constexpr (prepro_DEBUG) {
        if (get_eqLit(l) != l)
          log(0, "c ERROR: lit of softs should have been equality replaced:", l,
              " = ", get_eqLit(l), '\n');
        if (fixed_value(l) != l_Undef)
          log(0, "c ERROR: lit of softs is fixed, should have been removed:", l,
              " = ", fixed_value(l), '\n');
      }
    }
    // convert formula so that unit softs are of the form (-x) instead
    // of (x)---so making the 'blit' x true incurs the cost
    if (cls.size() == 1 && !sign(cls[0])) flipped_vars[var(cls[0])] = true;
  }

  Var nxtvar{0};
  ex2in.resize(nVars(), var_Undef);
  in2ex.resize(nVars(), var_Undef);
  for (Var v = 0; v != static_cast<Var>(nVars()); v++)
    if (appears[v]) {
      in2ex[nxtvar] = v;
      ex2in[v] = nxtvar;
      ++nxtvar;
    }
  maxvar = nxtvar - 1;

  MaxHS::Packed_vecs<Lit> tmp;
  vector<Lit> c;
  for (auto cls : hard_cls) {
    c.clear();
    for (auto l : cls) {
      assert(ex2in[var(l)] != var_Undef);
      c.push_back(
          mkLit(ex2in[var(l)], flipped_vars[var(l)] ? !sign(l) : sign(l)));
    }
    tmp.addVec(c);
  }
  hard_cls = std::move(tmp);
  tmp.clear();

  for (auto cls : soft_cls) {
    c.clear();
    for (auto l : cls) {
      assert(ex2in[var(l)] != var_Undef);
      c.push_back(
          mkLit(ex2in[var(l)], flipped_vars[var(l)] ? !sign(l) : sign(l)));
    }
    tmp.addVec(c);
  }
  soft_cls = std::move(tmp);

  for (auto& mx : mutexes) {
    for (auto& l : mx.soft_clause_lits_mod())
      l = mkLit(ex2in[var(l)], flipped_vars[var(l)] ? !sign(l) : sign(l));
    auto& el = mx.encoding_lit_mod();
    if (el != lit_Undef)
      el = mkLit(ex2in[var(el)], flipped_vars[var(el)] ? !sign(el) : sign(el));
  }
}

// INPUT PROBLEM/INTERNAL PROBLEM interface
// Take model found by solver and rewrite it into model of original formula
vector<lbool> Wcnf::rewriteModelToInput(const vector<lbool>& ubmodel) {
  // all original internal vars are preserved by solver. But
  // more vars might be added after.
  vector<lbool> ex_model(nOrigVars(), l_Undef);

  for (size_t i = 0; i < in2ex.size() && i < ubmodel.size(); ++i) {
    if (in2ex[i] != var_Undef && in2ex[i] < static_cast<Var>(nOrigVars())) {
      auto ex_var = in2ex[i];
      ex_model[ex_var] = flipped_vars[ex_var] ? ubmodel[i].neg() : ubmodel[i];
      if constexpr (prepro_DEBUG) {
        if (fixed_value(ex_var) != l_Undef)
          log(0,
              "c ERROR in rewriteModelToInput. fixed_values set for "
              "external "
              "variable\n");
        if (get_eqLit(mkLit(ex_var)) != mkLit(ex_var) ||
            get_eqLit(~mkLit(ex_var)) != ~mkLit(ex_var))
          log(0,
              "c ERROR in rewriteModelToInput. eqLits set for external "
              "variable\n");
      }
    }
  }

  // augment model with known fixed values. Note fixed_values have been
  // updated to account for equalities--so don't need to check value of
  // equivalent lit.
  for (size_t i = 0; i < ex_model.size(); i++) {
    if (ex_model[i] == l_Undef && fixed_value(i) != l_Undef) {
      ex_model[i] = fixed_value(i);
    }
  }

  // set values of variables removed by eq substitution
  // that do not have a fixed value.
  for (size_t i = 0; i < ex_model.size(); i++) {
    auto lt = mkLit(i);
    if (fixed_value(lt) == l_Undef && lt != get_eqLit(lt)) {
      if constexpr (prepro_DEBUG) {
        if (ex_model[i] != l_Undef)
          log(0,
              "c ERROR sat model sets variable that should have been "
              "equality "
              "substituted ",
              lt, " = ", get_eqLit(lt), '\n');
      }
      auto eq = get_eqLit(lt);
      if (ex_model[var(eq)] == l_Undef) ex_model[var(eq)] = l_False;
      if (sign(lt) == sign(eq))
        ex_model[i] = ex_model[var(eq)];
      else
        ex_model[i] = ex_model[var(eq)].neg();
    }
  }

  for (size_t i = 0; i < ex_model.size(); i++)
    if (ex_model[i] == l_Undef) ex_model[i] = l_False;

  if (using_maxpre) {
    vector<int> truelits;
    for (int i = 0; i < static_cast<int>(ex_model.size()); ++i) {
      int var = i + 1;
      if (ex_model[i] == l_True)
        truelits.push_back(var);
      else if (ex_model[i] == l_False)
        truelits.push_back(-var);
    }
    auto file_model = pif->reconstruct(truelits);
    ex_model.clear();
    for (auto tlit : file_model) {
      auto var = abs(tlit) - 1;
      lbool tv = (tlit > 0) ? l_True : l_False;
      if (static_cast<size_t>(var) >= ex_model.size())
        ex_model.resize(var + 1, l_False);
      ex_model[var] = tv;
    }
  }
  return ex_model;
}

// verify model against original formula.
Weight Wcnf::checkModel(const vector<lbool>& ubmodel, int& nfalseSofts) {
  assert(params.chkSoln);
  nfalseSofts = 0;
  bool hards_not_sat{false};
  vector<lbool> ex_model = rewriteModelToInput(ubmodel);
  for (auto hc : file_hards) {
    bool isSat = false;
    for (auto lt : hc) {
      auto var = abs(lt) - 1;
      if (((lt < 0) && ex_model[var] == l_False) ||
          ((lt > 0) && ex_model[var] == l_True)) {
        isSat = true;
        break;
      }
    }
    if (!isSat) hards_not_sat = true;
  }

  if (hards_not_sat) {
    cout << "c ERROR WCNF. Model does not satisfy the hards\n";
  }

  Weight w{0};
  int i = 0;
  for (auto sc : file_softs) {
    bool isSat = false;
    for (auto lt : sc) {
      auto var = abs(lt) - 1;
      if (((lt < 0) && ex_model[var] == l_False) ||
          ((lt > 0) && ex_model[var] == l_True)) {
        isSat = true;
        break;
      }
    }
    if (!isSat) {
      w += file_wts[i];
      ++nfalseSofts;
    }
    ++i;
  }
  return w;
}

// Stats and output
void Wcnf::computeWtInfo() {
  diffWts.clear();
  diffWtCounts.clear();
  transitionWts.clear();

  if (soft_clswts.size() == 0) {
    wt_min = wt_max = wt_mean = wt_var = 0;
    return;
  }

  wt_min = wt_max = 0;
  vector<Weight> wts(soft_clswts);
  std::sort(wts.begin(), wts.end());

  wt_min = wts.front();
  wt_max = wts.back();

  wt_mean = wt_var = 0;
  for (auto x : wts) wt_mean += x;
  total_cls_wt = wt_mean;

  wt_mean /= wts.size();
  for (auto x : wts) wt_var += (x - wt_mean) * (x - wt_mean);
  wt_var /= wts.size() - 1;

  size_t j{0};
  diffWts.push_back(wts[0]);
  diffWtCounts.push_back(0);
  for (size_t i = 0; i < wts.size(); i++)
    if (wts[j] == wts[i])
      diffWtCounts.back()++;
    else {
      j = i;
      diffWts.push_back(wts[i]);
      diffWtCounts.push_back(1);
    }

  Weight wtSoFar = diffWts[0] * diffWtCounts[0];
  for (size_t i = 1; i < diffWts.size(); i++) {
    if (diffWts[i] > wtSoFar) transitionWts.push_back(diffWts[i]);
    wtSoFar += diffWts[i] * diffWtCounts[i];
  }
}

void Wcnf::printFormulaStats() {
  auto n_units = nOrigUnits;
  cout << "c Instance: " << instance_file_name << "\n";
  cout << "c #vars: " << nVars() << "\n";
  cout << "c #Clauses: " << hard_cls.size() + n_units + soft_cls.size() << "\n";
  cout << "c Total Soft Cls Wt: " << wt_fmt(totalWt()) << '\n';
  cout << "c HARD: #Clauses = " << hard_cls.size() + n_units
       << ", Total Lits = " << hard_cls.total_size() + n_units << ", Ave Len = "
       << fix4_fmt((hard_cls.size() + n_units > 0)
                       ? (1.0 * hard_cls.total_size() + n_units) /
                           (hard_cls.size() + n_units)
                       : 0.0)
       << " #units = " << n_units << "\n";
  cout << "c SOFT: #Clauses = " << soft_cls.size()
       << ", Total Lits = " << soft_cls.total_size() << ", Ave Len = "
       << fix4_fmt(soft_cls.size()
                       ? static_cast<double>(soft_cls.total_size()) /
                           soft_cls.size()
                       : 0)
       << "\n";
  cout << "c Total Soft Clause Weight (+ basecost): " << wt_fmt(totalClsWt())
       << " (+ " << wt_fmt(baseCost()) << " from " << nfalse_softs
       << " false softs)\n";
  cout << "c SOFT%: "
       << fix4_fmt(soft_cls.size() + hard_cls.size() + n_units
                       ? (100.0 * soft_cls.size()) /
                           (soft_cls.size() + hard_cls.size() + n_units)
                       : 0)
       << "%\n";
  cout << "c #distinct weights: " << nDiffWts()
       << ", mean = " << fix4_fmt(aveSftWt())
       << ", std. dev = " << fix4_fmt(std::sqrt(varSftWt()))
       << ", min = " << wt_fmt(minSftWt()) << ", max = " << wt_fmt(maxSftWt())
       << "\n";
  cout << "c Parse time: " << fix4_fmt(parsing_time) << " secs.\n";
  cout << "c Wcnf Space Required: "
       << fix4_fmt(
              ((hard_cls.total_size() + soft_cls.total_size()) * sizeof(Lit) +
               soft_clswts.size() * sizeof(Weight)) /
              (1024 * 1024))
       << "MB\n";
  if (unsat) cout << "c Wcnf is UNSAT (hards are contradictory)\n";
  cout << "c ================================"
       << "\n";
}

void Wcnf::printSimpStats() {
  cout << "c After WCNF Simplification\n";
  cout << "c VARS: #vars = " << nVars() << '\n';
  cout << "c HARD: #Clauses = " << hard_cls.size()
       << ", Total Lits = " << hard_cls.total_size() << ", Ave Len = "
       << fix4_fmt((hard_cls.size() > 0)
                       ? (1.0 * hard_cls.total_size()) / hard_cls.size()
                       : 0.0)
       << "\n";
  cout << "c SOFT: #Clauses = " << soft_cls.size()
       << ", Total Lits = " << soft_cls.total_size() << ", Ave Len = "
       << fix4_fmt((1.0 * soft_cls.total_size()) / soft_cls.size()) << "\n";
  cout << "c Total Soft Clause Weight (+ basecost): " << wt_fmt(totalClsWt())
       << " (+ " << wt_fmt(baseCost()) << " from " << nfalse_softs
       << " false softs)\n";
  cout << "c #distinct weights: " << nDiffWts()
       << ", mean = " << fix4_fmt(aveSftWt())
       << ", std. dev = " << fix4_fmt(std::sqrt(varSftWt()))
       << ", min = " << wt_fmt(minSftWt()) << ", max = " << wt_fmt(maxSftWt())
       << "\n";
  cout << "c Total Clauses: " << hard_cls.size() + soft_cls.size() << "\n";
  cout << "c Wcnf Space Required: "
       << fix4_fmt(((hard_cls.total_size() + soft_cls.total_size() +
                     file_hards.total_size() + file_softs.total_size()) *
                        sizeof(Lit) +
                    (soft_clswts.size() + file_wts.size()) * sizeof(Weight)) /
                   1000000.0)
       << "MB\n";
  if (unsat) cout << "c Wcnf is UNSAT (hards are contradictory)\n";
  cout << "c ================================"
       << "\n";
}

void Wcnf::printFormula([[maybe_unused]] std::ostream& out) const {}

/*  // TODO modify to optionally output new DIMACS file.
  out << "c Wcnf---Print Formula\n";
  out << "c Dimacs (Vars, Clauses, TOP) = (" << dimacs_nvars << " ,"
      << dimacs_nclauses << " ," << dimacs_top << ")";
  out << " maxvar = " << nVars() << "\n";
  if (unsat) out << " formula is UNSAT\n";
  out << "c Hard Clauses # = " << hard_cls.size() + hard_units.size() <<
  "\n"; out << "c Soft Clauses, # = " << soft_cls.size() << "\n"; out << "c
  Base cost = " << wt_fmt(base_cost) << "\n"; out << "c HARD Units\n"; out
  << hard_units << "\n"; out << "c HARD SCCs\n"; out << all_scc << "\n"; out
  << "c HARDS\n"; out << hard_cls;

  out << "c SOFTS\n";
  for (size_t i = 0; i < soft_cls.size(); i++) {
    out << wt_fmt(soft_clswts[i]) << " ";
    for (auto& item : soft_cls[i]) out << item << " ";
    out << "0 \n";
  }
  }*/

void Wcnf::printFormula([[maybe_unused]] Bvars& bvars,
                        [[maybe_unused]] std::ostream& out) const {}
/*  // TODO modify to optionally output new DIMACS file.
  out << "c Wcnf---Print Formula\n";
  out << "c Dimacs (Vars, Clauses, TOP) = (" << dimacs_nvars << " ,"
      << dimacs_nclauses << " ," << wt_fmt(dimacs_top) << ")";
  out << " maxvar = " << nVars() << "\n";
  out << "c totalClsWt = " << wt_fmt(totalClsWt()) << "\n";
  if (unsat) out << " formula is UNSAT\n";
  out << "c Hard Clauses # = " << hard_cls.size() << "\n";
  out << "c Hard Units # = " << hard_units.size() << "\n";
  out << "c Hard SCC # = " << all_scc.size() << "\n";

  int nh{0};
  for (size_t i = 0; i < all_scc.size(); i++)
    out << "scc#" << nh++ << ": " << all_scc[i] << "\n";
  nh = 0;

  for (size_t i = 0; i < hard_units.size(); i++)
    out << "h#" << nh++ << ": " << hard_units[i] << "\n";

  for (size_t i = 0; i < nHards(); i++)
    out << "h#" << nh++ << ": " << getHard(i) << "\n";

  out << "c Soft Clauses, # = " << soft_cls.size() << "\n";
  out << "c Base cost = " << base_cost << "\n";

  int ns{0};
  for (size_t i = 0; i < nSofts(); i++)
    out << "c#" << ns++ << " blit = " << bvars.litOfCls(i)
        << " wt = " << wt_fmt(getWt(i)) << " : " << getSoft(i) << "\n";
        }*/

void Wcnf::printDimacs([[maybe_unused]] std::ostream& out) const {}
/*
//output wcnf in new standard dimacs
 out << "c maxhs preprocess output\n";
  if (unsat) {
    out << "c HARDS are UNSAT\n";
    out << "h -1 0\n";
    out << "h 1 0 \n";
    return;
  }

  if (baseCost() > 0) {
    out << wt_fmt(baseCost()) << " " << 1 << " 0\n";
    out << wt_fmt(baseCost()) << " " << -1 << " 0\n";
  }

  for (size_t i = 0; i < nSofts(); i++) {
    out << wt_fmt(soft_clswts[i]) << " ";
    for (auto l : soft_cls[i]) out << map_in2ex(l) << " ";
    out << "0\n";
  }

  if (!hard_units.empty())
    for (auto l : hard_units) out << l << " 0\n";

  for (size_t i = 0; i < nHards(); i++) {
    out << "h ";
    for (auto l : hard_cls[i]) out << map_in2ex(l) << " ";
    out << "0\n";
  }
  }*/
