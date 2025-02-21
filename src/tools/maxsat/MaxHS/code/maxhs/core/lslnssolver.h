/*****[lslsnsolver.h]
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

#ifndef LSLNSsolver_h
#define LSLNSsolver_h

#include <vector>

#ifdef GLUCOSE
#include "glucose/core/SolverTypes.h"
namespace Minisat = Glucose;
#else
#include "minisat/core/SolverTypes.h"
#endif

#include "maxhs/core/maxsolvertypes.h"
#include "maxhs/core/wcnf.h"

using Minisat::lbool;

namespace MaxHS {

class LSLNS_solver {
 public:
  LSLNS_solver(Wcnf* f);
  void                      solve();
  const std::vector<lbool>& get_best_model() { return best_model; }
  const Wcnf*               getWcnf() { return theWcnf; }
  int                       print_soln_stats();

 protected:
  Wcnf*              theWcnf;
  std::vector<lbool> best_model{};
  Weight             wt_best_soln{};
  bool               unsat{false};
  bool               get_initial_feasible_model();
  void               local_search();
  void               large_neighbourhood_search();
  uint64_t           total_flips{};
  uint64_t           total_epochs{};
  uint64_t           total_local_searches{};
  uint64_t           total_lns_searches{};
};
  
}  // namespace MaxHS
#endif
