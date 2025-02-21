/********************************************************************************************
Glucose3SolverProxy.h -- Copyright (c) 2020, Tobias Paxian

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************/

#ifndef GLUCOSE3SOLVERPROXY_H
#define GLUCOSE3SOLVERPROXY_H

#include <fstream>
#include "SATSolverProxy.h"
#include "glucose3/core/Solver.h"

class Glucose3SolverProxy : public SATSolverProxy {
 public:
  Glucose3SolverProxy();
  ~Glucose3SolverProxy();
  SATSolverType GetSATSolverType();
  unsigned int GetModel(int var);
  int NewVariable();
  void NewClause();
  void AddLiteral(unsigned int *lit);
  bool CommitClause();
  void ResetClause();
  void AddAssumption(unsigned int *lit);
  void ClearAssumption();
  unsigned int Solve();
  void Reset(void);
  bool AddClause(std::vector<unsigned int> &clause);
  unsigned int GetNumberOfClauses();
  unsigned int GetNumberOfVariables();
  //    void SetFrozen(int variable);

 protected:
  //    Glucose3::SimpSolver* _glucose3;
  Glucose3::Solver *_glucose3;
  Glucose3::vec<Glucose3::Lit> _currentClause;
  Glucose3::vec<Glucose3::Lit> _currentAssumptions;
};

#endif  // GLUCOSE3SOLVERPROXY_H
