/***********[maxhslns.cc]
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

#include <cerrno>
#include <csignal>
#include <iostream>

#ifdef GLUCOSE
#include "glucose/utils/System.h"
#else
#include "minisat/utils/System.h"
#endif

#include "maxhs/core/lslnssolver.h"
#include "maxhs/core/wcnf.h"
#include "maxhs/utils/options.h"
#include "maxhs/utils/params.h"

using std::cerr;
using std::cout;
using std::string;

#ifdef GLUCOSE
namespace Minisat = Glucose;
#endif
using namespace Minisat;

static MaxHS::LSLNS_solver* thesolver{};

static void SIGINT_exit([[maybe_unused]] int signum) {
  if (thesolver) thesolver->print_soln_stats();
  cout << std::endl;
  cerr << std::endl;
  // Note that '_exit()' rather than 'exit()' has to be used.
  // The reason is that 'exit()' calls destructors and may cause deadlocks
  // if a malloc/free function happens to be running (these functions are
  // guarded by
  //  locks for multithreaded use).
  _exit(0);
}

constexpr int majorVer{0};
constexpr int minorVer{1};
constexpr int update{0};

int main(int argc, char** argv) {
  const string usage{"USAGE: "};
  const string maxhs_lns_cat{"A: General maxhs_lns"};
  Options::set_usage_help(usage + argv[0] +
                          " [options] <input-file>\n  where input may be "
                          "either in plain or gzipped DIMACS.\n");
  Option<rlim_t> opt_cpu_lim(maxhs_lns_cat, "cpu-lim",
                         "Limit on CPU time allowed in seconds.\n",
                         Options::Unsigned_No_Limit);
  Option<rlim_t> opt_mem_lim(maxhs_lns_cat, "mem-lim",
                         "Limit on memory usage in megabytes.\n", INT32_MAX,
                         IntRange(0, INT32_MAX));
  BoolOption opt_version(maxhs_lns_cat, "version", "Print version number and exit",
                     false);

  if(!Options:parse_options(argc, argv) return 1;
  if (opt_version) {
    cout << "maxhs_lns " << majorVer << "." << minorVer << "." << update
         << "\n";
    return (0);
  }
  cout << "c maxhs_lns " << majorVer << "." << minorVer << "." << update
       << "\n";
  cout << "c Instance: " << argv[argc - 1] << "\n";
  if (params.printOptions) {
    cout << "c Parameter Settings\n";
    cout << "c ============================================\n";
    Options::print_option_settings(cout, "c ");
    cout << "c ============================================\n";
    cout << "c\n";
  }

  signal(SIGINT, SIGINT_exit);
  signal(SIGXCPU, SIGINT_exit);
  signal(SIGSEGV, SIGINT_exit);
  signal(SIGTERM, SIGINT_exit);
  signal(SIGABRT, SIGINT_exit);

  if (opt_cpu_lim != Options::Unsigned_No_Limit) {
    rlimit rl;
    getrlimit(RLIMIT_CPU, &rl);
    if (rl.rlim_max == RLIM_INFINITY || opt_cpu_lim < rl.rlim_max) {
      rl.rlim_cur = opt_cpu_lim;
      if (setrlimit(RLIMIT_CPU, &rl) == -1)
        cout << "c WARNING! Could not set resource limit: "
                "CPU-time.\n";
    }
  }

  if (opt_mem_lim != Options::Unsigned_No_Limit) {
    rlim_t new_mem_lim = opt_mem_lim * 1024 * 1024;
    rlimit rl;
    getrlimit(RLIMIT_AS, &rl);
    if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max) {
      rl.rlim_cur = new_mem_lim;
      if (setrlimit(RLIMIT_AS, &rl) == -1)
        cout << "c WARNING! Could not set resource limit: "
                "Virtual memory.\n";
    }
  }

  if (argc < 2) {
    cout << "c ERROR: no input file specfied:\n"
            "USAGE: %s [options] <input-file>\n  where input "
            "may be either "
            "in plain or gzipped DIMACS.\n";
    exit(0);
  }

  int ret_val{};
  try {
    Wcnf theFormula{};
    if (!theFormula.inputDimacs(argv[1])) return 1;
    MaxHS::LSLNS_solver S(&theFormula);
    thesolver = &S;
    thesolver->solve();
    ret_val = thesolver->print_soln_stats();
  } catch (const std::bad_alloc&) {
    cout << "c Memory Exceeded\n";
    ret_val = thesolver->print_soln_stats();
  } catch (...) {
    cout << "c Unknown exception probably memory.\n";
    ret_val = thesolver->print_soln_stats();
  }
  cout << std::endl;
  cerr << std::endl;
  return ret_val;
}
