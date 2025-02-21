/********************************************************************************************
CryptoMiniSat565SolverProxy.h -- Copyright (c) 2020, Tobias Paxian

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

#ifndef CRYPTOMINISATSOLVERPROXY_H
#define CRYPTOMINISATSOLVERPROXY_H

#include <fstream>  // std::ofstream
#include <vector>
#include "SATSolverProxy.h"
#include "cryptoMiniSat565/cryptominisat5/cryptominisat.h"

class CryptoMiniSATSolverProxy : public SATSolverProxy {
 public:
  CryptoMiniSATSolverProxy(unsigned int noOfThreads);
  ~CryptoMiniSATSolverProxy();

  // pure virtual functions - have to be implemented in the different
  // solverProxy's.
  SATSolverType GetSATSolverType(void);

  unsigned int GetModel(int var);
  int NewVariable();

  void NewClause();
  void AddLiteral(unsigned int *lit);
  bool CommitClause();
  void ResetClause();

  bool AddClause(std::vector<unsigned int> &clause);

  void AddAssumption(unsigned int *lit);
  void ClearAssumption();

  unsigned int Solve();

  void Reset(void);

  // the following functions can be implemented
  // if not it's either translated to the pure virtual functions
  // or it's just not necessary (some extended functionality of SAT solvers)
  void SetNumberOfThreads(unsigned int n);
  //    unsigned int GetNumberOfClauses();
  unsigned int GetNumberOfVariables();
  //  void AddLiteral(int* lit);
  //  void AddLiteral(int* lit, bool sign);
  //    void AddVariablePrio(unsigned int variable, unsigned int prio);
  //    virtual void SetFrozen(int variable);
  void NewVariables(unsigned number);

  // different solver configurations
  void SetReconf(unsigned reconf);

 protected:
  CMSat::SATSolver *_cmSAT;
  std::vector<CMSat::Lit> _currentClause;
  std::vector<CMSat::Lit> _currentAssumptions;
};

#endif  // CRYPTOMINISATSOLVERPROXY_H
