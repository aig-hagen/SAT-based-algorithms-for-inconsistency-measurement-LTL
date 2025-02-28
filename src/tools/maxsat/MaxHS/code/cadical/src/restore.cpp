#include "internal.hpp"

namespace CaDiCaL {

/*------------------------------------------------------------------------*/

// In incremental solving after a first call to 'solve' has finished and
// before calling the internal 'solve' again incrementally we have to
// restore clauses which have the negation of a literal as a witness literal
// on the extension stack, which was added as original literal in a new
// clause or in an assumption.  This procedure has to be applied
// recursively, i.e., the literals of restored clauses are treated in the
// same way as literals of a new original clause.
//
// To figure out whether literals are such witnesses we have a 'witness'
// bit for each external literal, which is set in 'block', 'elim', and
// 'decompose' if a clause is pushed on the extension stack.  The witness
// bits are recomputed after restoring clauses.
//
// We further mark in the external solver newly internalized external
// literals in 'add' and 'assume' since the last call to 'solve' as tainted
// if they occur negated as a witness literal on the extension stack.  Then
// we go through the extension stack and restore all clauses which have a
// tainted literal (and its negation a marked as witness).
//
// Since the API contract disallows to call 'val' and 'failed' in an
// 'UNKNOWN' state. We do not have to internalize literals there.
//
// In order to have tainted literals accepted by the internal solver they
// have to be active and thus we might need to 'reactivate' them before
// restoring clauses if they are inactive. In case they have completely
// been eliminated and removed from the internal solver in 'compact', then
// we just use a new internal variable.  This is performed in 'internalize'
// during marking external literals as tainted.
//
// To check that this approach is correct the external solver can maintain a
// stack of original clauses and current assumptions both in terms of
// external literals.  Whenever 'solve' determines that the current
// incremental call is satisfiable we check that the (extended) witness does
// satisfy the saved original clauses, as well as all the assumptions. To
// enable these checks set 'opts.check' as well as 'opts.checkwitness' and
// 'opts.checkassumptions' all to 'true'.  The model based tester actually
// prefers to enable the 'opts.check' option and the other two are 'true' by
// default anyhow.
//
// See our SAT'19 paper [FazekasBiereScholl-SAT'19] for more details.

/*------------------------------------------------------------------------*/

void External::restore_clause (
  const vector<int>::const_iterator & begin,
  const vector<int>::const_iterator& end) {
  LOG(begin, end, "restoring external clause");
  for (auto p = begin; p != end; p++) {
    int ilit = internalize(*p);
    internal->add_original_lit(ilit);
    internal->stats.restoredlits++;
  }
  internal->add_original_lit(0);
  internal->stats.restored++;
}

/*------------------------------------------------------------------------*/

void External::restore_clauses() {

  assert(internal->opts.restoreall == 2 || !tainted.empty());

  START(restore);
  internal->stats.restorations++;

  struct { int64_t weakened, satisfied, restored, removed; } clauses;
  memset(&clauses, 0, sizeof clauses);

  if (internal->opts.restoreall && tainted.empty())
    PHASE("restore", internal->stats.restorations,
          "forced to restore all clauses");

  unsigned numtainted = 0;
  for (const auto b : tainted)
    if (b) numtainted++;

  PHASE("restore", internal->stats.restorations,
    "starting with %u tainted literals %.0f%%",
    numtainted, percent (numtainted, 2u*max_var));

  auto end_of_extension = extension.end();
  auto p = extension.begin(), q = p;

  // Go over all witness labelled clauses on the extension stack, restore
  // those necessary, remove restored and flush satisfied clauses.
  //
  while (p != end_of_extension) {

    clauses.weakened++;

    assert(!*p);
    const auto saved = q;  // Save old start.
    *q++ = *p++;           // Copy zero '0'.

    // Copy witness part and try to find a tainted witness literal in it.
    //
    int tlit = 0;  // Negation tainted.
    int elit;
    //
    assert(p != end_of_extension);
    //
    while ((elit = *q++ = *p++)) {

      if (marked(tainted, -elit)) {
        tlit = elit;
        LOG("negation of witness literal %d tainted", tlit);
      }

      assert(p != end_of_extension);
    }

    // Now find 'end_of_clause' (clause starts at 'p') and at the same time
    // figure out whether the clause is actually root level satisfied.
    //
    int satisfied = 0;
    auto end_of_clause = p;
    while (end_of_clause != end_of_extension && (elit = *end_of_clause)) {
      if (!satisfied && fixed(elit) > 0)
        satisfied = elit;
      end_of_clause++;
    }

    // Do not apply our 'FLUSH' rule to remove satisfied (implied) clauses
    // if the corresponding option is set simply by resetting 'satisfied'.
    //
    if (satisfied && !internal->opts.restoreflush) {
      LOG(p, end_of_clause,
        "forced to not remove %d satisfied", satisfied);
      satisfied = 0;
    }

    if (satisfied || tlit || internal->opts.restoreall) {

      if (satisfied) {
        LOG(p, end_of_clause,
            "flushing implied clause satisfied by %d from extension stack",
            satisfied);
        clauses.satisfied++;
      } else {
        restore_clause(p, end_of_clause);  // Might taint literals.
        clauses.restored++;
      }

      clauses.removed++;
      p = end_of_clause;
      q = saved;

    } else {

      LOG(p, end_of_clause, "keeping clause on extension stack");

      while (p != end_of_clause)  // Copy clause too.
        *q++ = *p++;
    }
  }

  extension.resize(q - extension.begin());
  shrink_vector(extension);

#ifndef QUIET
  if (clauses.satisfied)
    PHASE("restore", internal->stats.restorations,
          "removed %" PRId64 " satisfied %.0f%% of %" PRId64 " weakened clauses",
          clauses.satisfied,
          percent(clauses.satisfied, clauses.weakened),
          clauses.weakened);
  else
    PHASE("restore", internal->stats.restorations,
          "no satisfied clause removed out of %" PRId64 " weakened clauses",
          clauses.weakened);

  if (clauses.restored)
    PHASE("restore", internal->stats.restorations,
          "restored %" PRId64 " clauses %.0f%% out of %" PRId64 " weakened clauses",
          clauses.restored,
          percent(clauses.restored, clauses.weakened),
          clauses.weakened);
  else
    PHASE("restore", internal->stats.restorations,
          "no clause restored out of %" PRId64 " weakened clauses",
          clauses.weakened);
#endif

  numtainted = 0;
  for (const auto& b : tainted)
    if (b) numtainted++;

  PHASE("restore", internal->stats.restorations,
    "finishing with %u tainted literals %.0f%%",
    numtainted, percent(numtainted, 2u * max_var));

  LOG("extension stack clean");
  tainted.clear();

  // Finally recompute the witness bits.
  //
  witness.clear();
  const auto begin_of_extension = extension.begin();
  p = extension.end();
  while (p != begin_of_extension) {
    while (*--p)
      assert(p != begin_of_extension);
    int elit;
    assert(p != begin_of_extension);
    while ((elit = *--p)) {
      mark(witness, elit);
      assert(p != begin_of_extension);
    }
  }

  STOP(restore);
}

}
