/*****[lslsnsolver.cc]
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
#include <vector>

#include "maxhs/core/lslnssolver.h"
#include "maxhs/ifaces/cadicalsatsolver.h"
#include "maxhs/ifaces/satsolver.h"
#include "maxhs/utils/options.h"

#ifdef GLUCOSE
namespace Minisat = Glucose;
#endif
using namespace Minisat;

using namespace MaxHS;
using MaxHS_Iface::cadicalSolver;
using MaxHS_Iface::SatSolver_uniqp;
using std::vector;

static const char* maxhs_lns = " maxhslns";
Option<unsigned> opt_epochs(
    maxhs_lns, "epochs",
    "number of epochs where each epoch is a local search followed "
    "by a large neighbour search",
    100, 1);
Option<unsigned> opt_ls_flips(maxhs_lns, "flips",
                              "number of flips allowed in "
                              "each local search episode",
                              100000, 1);
Option<unsigned> opt_lns_rounds(
    maxhs_lns, "lns-rounds",
    "maximum number of lns rounds before anothe local search ", 10, 1);
Option<double> opt_init_lns_fix(
    maxhs_lns, "init-lns-fix",
    "Percentage of "
    "fixed literals for first LNS computation (default 0.7)",
    0.7, 0.0, 1.0);
Option<unsigned> opt_lns_no_improve_before_fix_shrink(
    maxhs_lns, "lns_no_improve",
    "number of lns attempts (with new random "
    "neighbourhoods) before shrinking the percentage "
    "of fixed literals",
    10, 1);
Option<double> opt_lns_shrink_fix_factor(
    maxhs_lns, "lns-fix-shrink",
    "Reduce percentage fixed literals by this factor when no improvement "
    "occurs ",
    0.1, 0.0, 1.0);

LSLNS_solver::LSLNS_solver(Wcnf* f) : theWcnf{f} {
  params.instance_file_name = theWcnf->fileName();
}

void LSLNS_solver::solve() {
  if (!get_initial_feasible_model()) {
    unsat = true;
    return;
  }
  for (int i = 0; i < opt_epochs; i++) {
    local_search();
    large_neighbourhood_search();
  }
}

bool LSLNS_solver::get_initial_feasible_model() {
  SatSolver_uniqp satsolver{new cadicalSolver};
  Var bvar = theWcnf->nVars();
  satsolver->set_option("phase", false);
  for (size_t i = 0; i < theWcnf->nHards(); i++)
    satsolver->addClause(theWcnf->getHard(i));
  for (size_t i = 0; i < theWcnf->nSofts(); i++) {
    std::vector<Lit> sftCls{theWcnf->getSoft(i)};
    sftCls.push_back(Minisat::mkLit(bvar++));
    satsolver->addClause(sftCls);
  }
  best_model.clear();
  if (satsolver->solve() == l_False)
    return false;
  else {
    best_model.resize(theWcnf->nVars(), l_Undef);
    for (size_t i = 0; i < theWcnf->nVars(); i++)
      best_model[i] = satsolver->modelValue(i);
    return true;
  }
}

void LSLNS_solver::local_search() { total_local_searches++; }
int LSLNS_solver::print_soln_stats() {
  if (best_model.empty()) {
    cout << "c no model found\n";
    cout << "s UNKNOWN\n";
    return 0;
  }
  if (unsat) {
    cout << "c formula has no feasible models\n";
    cout << "s UNSATISFIABLE\n";
    return 20;
  }
  std::vector<lbool> external_model = theWcnf->rewriteModelToInput(best_model);
  int nfalseSofts;
  auto w = theWcnf->checkModelFinal(best_model, nfalseSofts);
  cout << "o " << wt_fmt(w) << '\n';
  if (params.printSoln) {
    if (params.printNewFormat) {
      std::string vline(external_model.size() + 2, ' ');
      vline[0] = 'v', vline[1] = ' ';
      for (size_t i = 0; i < external_model.size(); i++)
        vline[i + 2] = (external_model[i] == l_True ? '1' : '0');
      cout << vline << std::endl;
    } else {
      cout << 'v';
      for (size_t i = 0; i < external_model.size(); i++) {
        cout << ' ' << (external_model[i] == l_True ? "" : "-") << i + 1;
      }
      cout << std::endl;
    }
  }
  return 10;
}
