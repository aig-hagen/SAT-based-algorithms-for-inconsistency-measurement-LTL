/***********[maxhs.cc]
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
#include <string>

#ifdef GLUCOSE
#include "glucose/utils/System.h"
#else
#include "minisat/utils/System.h"
#endif

#include "maxhs/core/maxsolver.h"
#include "maxhs/core/wcnf.h"
#include "maxhs/utils/options.h"
#include "maxhs/utils/params.h"

using std::cout;
using std::string;

#ifdef GLUCOSE
namespace Minisat = Glucose;
#endif
using namespace Minisat;
using namespace MaxHS;

static MaxSolver* thesolver{};

static void SIGINT_exit(int signum) {
  if (thesolver) {
    thesolver->printStatsAndExit(signum, 1);
  } else {
    fflush(stdout);
    fflush(stderr);
    std::cout << std::flush;
    std::cerr << std::flush;
    // Note that '_exit()' rather than 'exit()' has to be used.
    // The reason is that 'exit()' calls destructors and may cause deadlocks
    // if a malloc/free function happens to be running (these functions are
    // guarded by
    //  locks for multithreaded use).
    _exit(0);
  }
}

constexpr int majorVer{5};
constexpr int minorVer{0};
constexpr int update{0};

int main(int argc, char** argv) {
  if (!params.sverbosity && !params.cplex_output) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.tie(0);
    std::cerr.tie(0);
  }
  std::cout << std::setprecision(0);
  if (params.flush_io) std::cout << std::unitbuf;

  const string usage{"USAGE: "s};
  const string maxhs_cat{"A: General MaxHS"s};
  Options::set_usage_help(
      usage + argv[0] +
      " [options] <input-file>,"
      "\n\twhere input is DIMACS instance either plain or gzipped text."
      "\n\tIf input missing read instance from <stdin>."
      "\n\tFor options use -<option>=<value> with equal sign and no spaces."
      "\n\t-option or -no-<option> for Boolean options.");

  Option<rlim_t> opt_cpu_lim(maxhs_cat, "cpu-lim",
                             "Limit on CPU time allowed in seconds.",
                             Options::Unsigned_No_Limit);
  Option<rlim_t> opt_mem_lim(maxhs_cat, "mem-lim",
                             "Limit on memory usage in megabytes.",
                             Options::Unsigned_No_Limit);
  Option<bool> opt_version(maxhs_cat, "version",
                           "Print version number and exit", false);
  Option<string> opt_filename(
      maxhs_cat, "filename",
      "specify instance filename (useful when reading from stdin)", "<stdin>");

  if (!Options::parse_options(argc, argv)) return 1;
  if (argc > 2) {
    cout << "Error. Found extra arguments on command line. Expected only "
            "filename after options\n";
  }
  params.readOptions();
  cout << "c MaxHS " << majorVer << "." << minorVer << "." << update << "\n";
  if (opt_version) { return (0); }
  cout << "c Instance: ";
  if (argc == 2)
    cout << argv[1] << "\n";
  else
    cout << string{opt_filename} << "\n";

  if (params.printOptions) {
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

  try {
    Wcnf theFormula{};
    if (argc == 2) {
      if (!theFormula.inputDimacs(argv[1], argv[1], true)) return 1;
    } else {
      if (!theFormula.inputDimacs("", opt_filename, true)) return 1;
    }

    /*cout << "Different wts " << theFormula.nDiffWts() << " =
    " << theFormula.getDiffWts() << "\n"; cout << "Mutexes ";
    for(size_t i=0; i < theFormula.getMxs().size(); i++) {
    cout << "#" << i << ". "; if(theFormula.isCoreMx(i)) cout
    << "Core Mx     "; else cout << "Non-Core Mx "; cout <<
    theFormula.getMxs()[i] << "\n";
    }
    theFormula.printFormula();*/

    MaxHS::MaxSolver S(&theFormula);
    thesolver = &S;
    S.solve();
    S.printStatsAndExit(-1, 0);
  } catch (const std::bad_alloc&) {
    cout << "c Memory Exceeded\n" << std::flush;
    thesolver->printStatsAndExit(100, 1);
  } catch (...) {
    cout << "c Unknown exception probably memory.\n" << std::flush;
    thesolver->printStatsAndExit(200, 1);
  }
  std::cout << std::flush;
  std::cerr << std::flush;
  fflush(stdout);
  fflush(stderr);
  return 0;
}
