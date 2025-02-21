/*****[params.cc]
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
#include <iostream>
#include <limits>
using std::cout;

#include "maxhs/utils/options.h"
#include "maxhs/utils/params.h"

using namespace MaxHS;

static const char* maxhs = "A: General MaxHS";
static const char* abstract = "B: MaxHS with abstractions";
static const char* specialSolvers = "C: MaxHS special solvers";
static const char* disjoint = "D: Disjoint Phase";
static const char* seed = "E: Seeding";
static const char* seqOfSat = "F: Sequence of Sat";
static const char* muser = "G: Core Minimization";
static const char* cplex = "H: CPLEX";
static const char* pop = "I: CPLEX Solution Pool and Populate";
static const char* pre = "J: Preprocessing";
static const char* debug = "K: Debugging";

// General Controls
Option<unsigned> opt_verb(maxhs, "verb",
                          "Verbosity level (0=silent, 1=some, 2=more, "
                          "3=debugging output, 4=more debugging output).",
                          1, 0, 5);
Option<unsigned> opt_sverb(
    maxhs, "sverb",
    "Sat solver verbosity level (0=silent, 1=some, 2=more,3=debugging output, "
    "4=more debugging output).",
    0, 0, 4);

Option<bool> opt_fbeq(maxhs, "fbeq", "Use Fbeq theory", false);

Option<bool> opt_printOptions(maxhs, "printOptions", "Print paramater settings",
                              false);
Option<bool> opt_printBstSoln(maxhs, "printBstSoln",
                              "Print best solution found", false);
Option<bool> opt_printSoln(maxhs, "printSoln", "Print solution", false);
Option<bool> opt_printOldFormat(maxhs, "printSoln-old-format",
                                "Print solution in old format", false);
Option<bool> opt_chkSoln(
    maxhs, "chkSoln",
    "Check that returned solution is feasible using original input formula. "
    "This requires storing an unchanged copy of the input instance--for larger "
    "instances this can occupy up to 1GB",
    true);

Option<double> opt_tolerance(
    maxhs, "tolerance",
    "For floating point weights only: return solution when when |soln-cost - "
    "lower bound| <= tolerance",
    1e-6, 0.0);
Option<bool> opt_flush_io(maxhs, "flush-io",
                          "Always flush cout (for debugging)", false);

// Muser Controls
Option<unsigned> opt_mintype(muser, "mintype",
                             "JD: 0 = no minimization of "
                             "constraints/cores,  1 = Use Muser",
                             1, 0, 1);
Option<double> opt_mus_cpu_lim(
    muser, "mus-cpu-lim",
    "FB: CPU time limit for minimizing each core (-1 == no limit).", 2.5,
    Options::No_Limit);
Option<double> opt_mus_min_red(
    muser, "mus-min-red",
    "FB: Run muser only if on average it can remove at least this "
    "fraction of a core (-1 == no limit). (eventually the muser is turned off)",
    0.10, Options::No_Limit, 1.0);
Option<unsigned> opt_muser_verb(
    muser, "mverb",
    "Muser verbosity level (0=silent, 1=some, 2=more,3=debugging output, "
    "4=more debugging output).",
    0, 0, 4);

Option<bool> opt_abstract(abstract, "abstract", "JB: abstract cores", true);
Option<double> opt_abstract_max_ave_size(
    abstract, "abstract-max-ave-size",
    "Don't do abstractions if average core size is greater than this limit "
    "(-1==no limit)",
    100.0, Options::No_Limit);
Option<unsigned> opt_abstract_cplex_cores(
    abstract, "abstract-cplex_cores",
    "Generate cores from abstracted cplex solution (0=cores from non-abstract "
    "soln only, 1=cores "
    "from abstract soln only, 2=cores from both abstract and non-abstract soln",
    2, 0, 2);
Option<unsigned> opt_abstract_greedy_cores(
    abstract, "abstract-greedy_cores",
    "Generate cores from abstracted greedy solution (0=cores from non-abstract "
    "soln only, 1=cores "
    "from abstract soln only, 2=cores from both abstract and non-abstract soln",
    2, 0, 2);
Option<unsigned> opt_abstract_min_size(
    abstract, "abstract-minsize",
    "JB: minimum-size of summation before adding them", 2, 1);
Option<unsigned> opt_abstract_max_core_size(
    abstract, "abstract-max_core_size",
    "JB: max-size of core to consider when abstracting", 1000, 1);
Option<unsigned> opt_abstract_min_cores(
    abstract, "abstract-min-cores",
    "Only allow softs to be clustered into abstractions when they appear in "
    "this minimum number of cores",
    2);

Option<double> opt_abs_cpu(abstract, "abs-cpu",
                           "CPU limit for abstraction round", 256.0,
                           Options::No_Limit);

Option<double> opt_cpu_per_exhaust(
    abstract, "exhaust-cpu-lim",
    "JB: CPU time limit for exhausting summations (-1 == no limit).", 3.0,
    Options::No_Limit);
Option<double> opt_abstract_gap(
    abstract, "abstract-gap",
    "If the lp-relaxation does not improve by this we consider doing "
    "abstraction ",
    1.0);

Option<double> opt_initial_abstract_gap(
    abstract, "1st-abstract-gap",
    "If seeding and initial disjoint does not improve the lp-relaxation gap by "
    "this amount we consider doing abstraction",
    5.0, 0.0);

// special solver controls
Option<bool> opt_use_hs_solver(
    specialSolvers, "use-hs-solver",
    "Use a specialized solver for hitting set problems", true);
Option<double> opt_hs_cpu(specialSolvers, "hs-cpu",
                          "When solving a hitting set problem this is the cpu "
                          "limit for the special solver",
                          Options::No_Limit, Options::No_Limit);

Option<bool> opt_use_lsu_solver(
    specialSolvers, "use-lsu-solver",
    "Try lsu on problems with small numbers of softs", true);
Option<double> opt_unwt_lsu_cpu(
    specialSolvers, "unwt-lsu-cpu",
    "When solving an unweighted lsu problem this is the cpu limit for the special solver",
    50.0, Options::No_Limit);
Option<double> opt_wt_lsu_cpu(
    specialSolvers, "wt-lsu-cpu",
    "When solving an unweighted lsu problem this is the cpu limit for the special solver",
    1500.0, Options::No_Limit);
Option<int> opt_lsu_var_limit(
    specialSolvers, "lsu-max-vars",
    "Try lsu when problems no more than this number of softs", 256);

// Disjoint Phase Controls
Option<bool> opt_dsjnt(disjoint, "dsjnt",
                       "JD: Find disjoint cores in a first phase.", true);

Option<double> opt_dsjnt_cpu_for_hs(
    disjoint, "dsjnt",
    "When solving a hitting set problem this is the maximum cpu time spent on "
    "finding disjoint cores",
    250, Options::No_Limit);

Option<double> opt_dsjnt_cpu_per_core(
    disjoint, "dsjnt-cpu-lim",
    "FB: CPU time limit for finding each disjoint cores (-1 == no limit).",
    30.0, Options::No_Limit);

Option<double> opt_dsjnt_mus_cpu_lim(
    disjoint, "dsjnt-mus-cpu-lim",
    "FB: CPU time limit for minimizing each *disjoint* core (-1 == no limit).",
    10.0, Options::No_Limit);

// Noncore and Seeding Options
Option<bool> opt_seed_learnts(
    seed, "seed-learnts",
    "FB: seed any learnts available when seeding is performed.", true);

Option<unsigned> opt_coretype(
    maxhs, "coretype",
    "JD: Type of constraints to learn and"
    " feed to CPLEX (0 = core constraints only) (1 = mixed constraints).",
    0, 0, 1);
Option<unsigned> opt_seedtype(
    seed, "seedtype",
    "FB: Type of seeded constraints allowed, 0 = no seeding, 1 = cores only, 2 "
    "= also allow non-cores, 3 = also allow mixed constraints",
    3, 0, 3);
Option<int> opt_maxseeds(seed, "seed-max",
                         "FB: maximum number of seeded constraints",
                         1024 * 512);

Option<unsigned> opt_seed_all_limit(
    seed, "seed-all-limit",
    "If the total number of variables is <= this limit and the total number of "
    "clauses <= 64* this limit) then seed all clauses into "
    "CPLEX (subject to \"seed-max\" limit...CPLEX will try to solve "
    "but SAT might also be used",
    256 * 2);

Option<double> opt_seed_all_cpu_before_cplex(
    seed, "seed_cpu_before_cplex",
    "CPU time limit before calling cplex when all clauses seeded", 200.0,
    Options::No_Limit);

Option<double> opt_all_seeded_first_cplex_cpu(
    seed, "all-seeded-1st-cplex-cpu",
    "CPU limit for first cplex solve when all clauses seeded", 100.0,
    Options::No_Limit);

Option<double> opt_all_seeded_first_abs_cpu(
    seed, "all-seeded-1st-abs-cpu",
    "CPU limit first abstraction when all clauses seeded", 60.0,
    Options::No_Limit);

Option<double> opt_all_seeded_2nd_abs_cpu(
    seed, "all-seeded_2nd_abs_cpu",
    "CPU limit second abstraction when all clauses seeded", 240.0,
    Options::No_Limit);

// Populate and Soln Pool
//   limit
Option<unsigned> opt_cplex_solnpool_cap(
    pop, "cplex-solnpool-cap", "Set the capacity of cplex solution pool", 256);

Option<unsigned> opt_cplex_pop_nsoln(pop, "cplex-pop-nsoln",
                                     "Set the size of cplex population pool",
                                     512 / 2);
Option<double> opt_cplex_pop_cpu_lim(
    pop, "cplextime-pop-cpu-lim",
    "CPU time limit on cplex populate (-1 == no limit)", 7.5,
    Options::No_Limit);

Option<unsigned> opt_trypop(pop, "cplex-populate",
                            "Use cplex populate to obtain more solutions "
                            "(0=never) (1=when potentially useful) (2=always)",
                            1, 0, 2);

Option<unsigned> opt_conflicts_from_ub(
    pop, "ub-conflicts",
    "FB: Generate conflicts from upper bound (0=neve) (1=when potentially "
    "useful) (2=always)",
    1, 0, 2);

// Sequence of Sat Options
Option<double> opt_optcores_cpu_per(
    seqOfSat, "optcores-cpu-lim",
    "FB: CPU time limit for finding each additional core (-1 == no limit).", 10,
    Options::No_Limit);

Option<unsigned> opt_nonopt(
    seqOfSat, "nonopt",
    "JD: Method for relaxing clauses of current "
    "core (0 = pick a random clause, 1 = pick clause appearing in most cores"
    ", 2 = relax a fraction of each core (set fraction with \"relaxfrac\" "
    "parameter), 3 = remove all clauses in core making next core disjoint.",
    3, 0, 3);

Option<unsigned> opt_abstract_assumps(
    seqOfSat, "abstract-assumps",
    "Method for relaxing abstract assumptions (0 = remove summation outputs "
    "like ordinary b-vars, 1 = relax summations to be next output, 2 = relax "
    "only one summation output at a time",
    1, 0, 2);

Option<double> opt_relaxfrac(
    seqOfSat, "relaxfrac",
    "FB: After accumulating frac-rampup-end clauses "
    "relax this fraction of current core, picking clauses most frequently "
    "occuring in cores (must have \"nonopt=2\").",
    0.3, 0.0, 1.0);
Option<unsigned> opt_frac_rampup_start(
    seqOfSat, "frac-rampup-start",
    "FB: When nonopt = 2 (relax a fraction) relax only"
    " one clause until this many cores accumulated",
    128);
Option<unsigned> opt_frac_rampup_end(
    seqOfSat, "frac-rampup-end",
    "FB: When nonopt = 2 (relax a fraction) increase "
    "fract of core relaxed linearly to reach final \"relaxfrac\"  after this "
    "many cores accumulated",
    512);
Option<unsigned> opt_max_cores_before_cplex(
    seqOfSat, "max-cores-before-cplex",
    "FB: Force a call to Cplex after this many constraints", 300);

Option<double> opt_max_cpu_before_cplex(
    seqOfSat, "max-cpu-before-cplex",
    "FB: Force a call to Cplex after this many CPU seconds (-1 == no limit)",
    200.0, Options::No_Limit);

Option<bool> opt_b_m_s(
    seqOfSat, "use-ub-mipstart",
    "FB: Use current Sat solver upper bound model as cplex start. This entails "
    "deleting all other starts",
    true);

Option<unsigned> opt_sort_assumps(seqOfSat, "sort-assumps",
                                  "FB: (0=don't sort, 1=place best softs to "
                                  "relax at top of trail, 2 reverse of 1)",
                                  0, 0, 2);

Option<bool> opt_find_forced(
    seqOfSat, "find-forced",
    "Check for forced variables by UP or by the upper bound", false);

Option<bool> opt_lp_harden(
    seqOfSat, "lp-harden",
    "Use LP version of CPLEX model to force soft clauses", true);

// Other Controls
// CPLEX Solver Options
Option<unsigned> opt_cplex_threads(
    cplex, "cplex-threads",
    "Allow cplex to use this many threads (1 = sequential processing)", 1, 1,
    124);
Option<bool> opt_cplex_tune(
    cplex, "cplex-tune",
    "Use cplex parameter setting recommended by cplex-tune", false);
Option<double> opt_cplex_min_ticks(
    cplex, "cplex-min-ticks",
    "Run CPLEX for at least this 1000's of its deterministic ticks can allow "
    "CPLEX to find better feasible (non-optimal) solutions",
    4.0, 1.0);

// //Preprocessing
Option<bool> opt_preprocess(pre, "preprocess", "Use minisat preprocessor",
                            true);
Option<bool> opt_prepro_wcnf_eqs(pre, "wcnf-eqs",
                                 "Find and reduce equalities in wcnf", true);
Option<bool> opt_prepro_wcnf_units(pre, "wcnf-units",
                                   "Reduce wcnf by hard units", true);
Option<bool> opt_prepro_wcnf_harden(
    pre, "wcnf-harden", "Try to harden soft clauses by satisfiability tests",
    true);

// //Mutexes
Option<unsigned> opt_prepro_mx_find_mxes(
    pre, "mx-find-mxes",
    "Detect mutually exclusive soft clauses in the input formula (0=don't, 1= "
    "find at most one false (core-mxes), 2= find at most one true "
    "(non-core-mxes), 3=1&2)",
    2, 0, 3);

Option<unsigned> opt_prepro_mx_mem_lim(
    pre, "mx-mem-lim", "Limit on memory usage in megabytes of the mutex finder",
    512 * 3);

Option<bool> opt_prepro_use_maxpre(pre, "use-maxpre",
                                   "Try to simplify formula by running maxpre",
                                   true);

Option<double> opt_prepro_maxpre_time(pre, "maxpre-time",
                                      "Time to allocate to maxpre", 250, 0);

Option<bool> opt_prepro_simplify_and_exit(
    pre, "simplify-only",
    "Write simplified WCNF file with new suffix then exit. If "
    "mx-exit-if-no-mutexes we exit before writing if no mutexes found",
    false);

Option<bool> opt_prepro_mx_seed_originals(
    pre, "mx-seed-mxes",
    "Allow original softs clauses in mutexes to be seeded to CPLEX when "
    "formula is transformed",
    true);  // true

Option<bool> opt_prepro_mx_constrain_hs(
    pre, "mx-constrain-hs",
    "Ensure that computed hitting sets satisfy the discovered soft clause "
    "mutexes",
    true);  // true

Option<double> opt_prepro_mx_max_cpu(
    pre, "mx-cpu-lim", "Max time to spend on mx detection (-1 == no limit)",
    15.0, Options::No_Limit);

// //Debugging
Option<bool> opt_cplex_data_chk(debug, "cplex-data-chk",
                                "Run cplex data checker on its input", true);
Option<bool> opt_cplex_write_model(debug, "cplex-wrt-model",
                                   "Make cplex write out each of its models",
                                   false);
Option<bool> opt_cplex_output(debug, "cplex-output", "Turn on cplex output",
                              false);

void Params::readOptions() {
  verbosity = opt_verb;
  sverbosity = opt_sverb;
  mverbosity = opt_muser_verb;

  printOptions = opt_printOptions;
  printBstSoln = opt_printBstSoln;
  printOldFormat = opt_printOldFormat;
  printSoln = opt_printSoln;
  flush_io = opt_flush_io;
  if (printOldFormat) printSoln = true;
  chkSoln = opt_chkSoln;
  tolerance = opt_tolerance;
  min_type = opt_mintype;
  mus_cpu_lim = (opt_mus_cpu_lim > 0) ? opt_mus_cpu_lim : noLimit;
  mus_min_red = (opt_mus_min_red > 0) ? opt_mus_min_red : noLimit;

  use_hs_slv = opt_use_hs_solver;
  hs_cpu = opt_hs_cpu;
  use_lsu_slv = opt_use_lsu_solver;
  lsu_unwt_cpu = opt_unwt_lsu_cpu;
  lsu_wt_cpu = opt_wt_lsu_cpu;
  lsu_max_vars = opt_lsu_var_limit;

  dsjnt_phase = opt_dsjnt;
  dsjnt_cpu_for_hs =
      (opt_dsjnt_cpu_for_hs > 0) ? opt_dsjnt_cpu_for_hs : noLimit;
  dsjnt_cpu_per_core =
      (opt_dsjnt_cpu_per_core > 0) ? opt_dsjnt_cpu_per_core : noLimit;
  dsjnt_mus_cpu_lim =
      (opt_dsjnt_mus_cpu_lim > 0) ? opt_dsjnt_mus_cpu_lim : noLimit;
  optcores_cpu_per =
      (opt_optcores_cpu_per > 0) ? opt_optcores_cpu_per : noLimit;
  find_forced = opt_find_forced;

  seed_type = opt_seedtype;
  seed_max = opt_maxseeds;
  seed_learnts = opt_seed_learnts;
  seed_all_limit = opt_seed_all_limit;
  seed_all_cpu = opt_seed_all_cpu_before_cplex;
  frac_to_relax = opt_relaxfrac;
  frac_rampup_start = opt_frac_rampup_start;
  frac_rampup_end = opt_frac_rampup_end;
  max_cores_before_cplex = opt_max_cores_before_cplex;
  max_cpu_before_cplex = opt_max_cpu_before_cplex;
  //  max_cplex_calls_before_opt = opt_max_cplex_calls_before_optimum;
  lp_harden = opt_lp_harden;

  sort_assumps = opt_sort_assumps;
  bestmodel_mipstart = opt_b_m_s;
  fbeq = opt_fbeq;

  abstract_assumps = opt_abstract_assumps;
  int no = opt_nonopt;
  switch (no) {
    case 0:
      coreRelaxFn = CoreRelaxFn::rand;
      break;
    case 1:
      coreRelaxFn = CoreRelaxFn::maxoccur;
      break;
    case 2:
      coreRelaxFn = CoreRelaxFn::frac;
      break;
    case 3:
      coreRelaxFn = CoreRelaxFn::dsjn;
      break;
  }

  int ct = opt_coretype;
  switch (ct) {
    case 0:
      coreType = CoreType::cores;
      break;
    case 1:
      coreType = CoreType::mixed;
      break;
    default:
      coreType = CoreType::cores;
      break;
  }

  // cplex limits
  cplex_threads = opt_cplex_threads;
  cplex_tune = opt_cplex_tune;
  cplex_min_ticks = opt_cplex_min_ticks;
  cplex_data_chk = opt_cplex_data_chk;
  cplex_write_model = opt_cplex_write_model;
  cplex_output = opt_cplex_output;

  // cplex_solnpool_cap = opt_cplex_solnpool_cap;
  cplex_pop_nsoln = opt_cplex_pop_nsoln;
  cplex_pop_cpu_lim =
      (opt_cplex_pop_cpu_lim > 0) ? opt_cplex_pop_cpu_lim : noLimit;
  // trypop_cplextime_ub = (opt_trypop_cplextime_ub > 0) ?
  // opt_trypop_cplextime_ub : noLimit; trypop_feedtime_lb =
  // opt_trypop_feedtime_lb;
  //  trypop = opt_trypop;
  trypop = opt_trypop;
  conflicts_from_ub = opt_conflicts_from_ub;

  if (opt_cplex_pop_nsoln == 0) trypop = 0;

  // preprocess
  preprocess = opt_preprocess;
  wcnf_eqs = opt_prepro_wcnf_eqs;
  wcnf_harden = opt_prepro_wcnf_harden;
  wcnf_units = opt_prepro_wcnf_units;

  simplify_and_exit = opt_prepro_simplify_and_exit;
  mx_find_mxes = opt_prepro_mx_find_mxes;
  mx_mem_limit = opt_prepro_mx_mem_lim;
  use_maxpre = opt_prepro_use_maxpre;
  maxpre_time = opt_prepro_maxpre_time;
  mx_seed_originals = opt_prepro_mx_seed_originals;
  mx_constrain_hs = opt_prepro_mx_constrain_hs;
  mx_cpu_lim = (opt_prepro_mx_max_cpu > 0) ? opt_prepro_mx_max_cpu : noLimit;

  abstract = opt_abstract;
  abstract_max_ave_size = opt_abstract_max_ave_size;
  abstract_cplex_cores = opt_abstract_cplex_cores;
  abstract_greedy_cores = opt_abstract_greedy_cores;
  abstract_min_size = opt_abstract_min_size;
  abstract_max_core_size = opt_abstract_max_core_size;
  abstract_min_cores = opt_abstract_min_cores;
  all_seeded_first_cplex_cpu = opt_all_seeded_first_cplex_cpu;
  all_seeded_first_abs_cpu = opt_all_seeded_first_abs_cpu;
  all_seeded_2nd_abs_cpu = opt_all_seeded_2nd_abs_cpu;
  cpu_per_exhaust = opt_cpu_per_exhaust;
  abstract_gap = opt_abstract_gap;
  initial_abstract_gap = opt_initial_abstract_gap;
  abs_cpu = opt_abs_cpu;
}

Params params;
