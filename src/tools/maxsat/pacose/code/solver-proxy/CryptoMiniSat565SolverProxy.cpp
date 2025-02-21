/********************************************************************************************
CryptoMiniSat565SolverProxy.cpp -- Copyright (c) 2020, Tobias Paxian

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

#include "CryptoMiniSat565SolverProxy.h"
#include <vector>
#include "SATSolverProxy.h"
#include "cryptoMiniSat565/cryptominisat5/cryptominisat.h"
#include "iostream"

CryptoMiniSATSolverProxy::CryptoMiniSATSolverProxy(unsigned int noOfThreads)
    : SATSolverProxy(noOfThreads), _cmSAT(new CMSat::SATSolver) {
  Reset();
}

CryptoMiniSATSolverProxy::~CryptoMiniSATSolverProxy() { delete _cmSAT; }

SATSolverType CryptoMiniSATSolverProxy::GetSATSolverType() {
  std::cout << CMSat::SATSolver::get_version();
  return SATSolverType::CRYPTOMINISAT;
}

unsigned int CryptoMiniSATSolverProxy::GetModel(int var) {
  CMSat::lbool varLbool = _cmSAT->get_model()[static_cast<unsigned int>(var)];

  unsigned int rv;
  //    std::cout << "varLbool " << varLbool << std::endl;
  if (varLbool == CMSat::l_True) {
    rv = static_cast<unsigned int>(var << 1);
  } else if (varLbool == CMSat::l_False) {
    rv = static_cast<unsigned int>(var << 1 ^ 1);
  } else if (varLbool == CMSat::l_Undef) {
    std::cout << "Strange Value, l_Undef for variable: " << var << std::endl;
    //         undefined means, both values are possible!
    rv = static_cast<unsigned int>(var << 1);
  } else {
    std::cout << "Strange value for varLbool: " << varLbool << std::endl;
    assert(false);
  }
  return rv;
}

int CryptoMiniSATSolverProxy::NewVariable() {
  // TODO: maybe to integrate new vars - to get multiple new vars at once!
  _cmSAT->new_var();
  return static_cast<int>(_cmSAT->nVars() - 1);
}

void CryptoMiniSATSolverProxy::NewClause() { assert(_currentClause.empty()); }

void CryptoMiniSATSolverProxy::AddLiteral(unsigned int *lit) {
  _currentClause.push_back(CMSat::Lit::toLit(*lit));
}

bool CryptoMiniSATSolverProxy::CommitClause() {
  return _cmSAT->add_clause(_currentClause);
}

void CryptoMiniSATSolverProxy::ResetClause() { _currentClause.clear(); }

bool CryptoMiniSATSolverProxy::AddClause(std::vector<unsigned int> &clause) {
  std::vector<CMSat::Lit> tmpClause;
  tmpClause.reserve(clause.size());
  for (unsigned int literal : clause) {
    tmpClause.push_back(CMSat::Lit::toLit(literal));
  }
  return _cmSAT->add_clause(tmpClause);
}

void CryptoMiniSATSolverProxy::AddAssumption(unsigned int *lit) {
  _currentAssumptions.push_back(CMSat::Lit::toLit(*lit));
}

void CryptoMiniSATSolverProxy::ClearAssumption() {
  _currentAssumptions.clear();
}

unsigned int CryptoMiniSATSolverProxy::Solve() {
  if (_cpuLimit > 0) {
    _cmSAT->set_timeout_all_calls(_cpuLimit);
  }

  EnableMemoryLimit();

  CMSat::lbool solution = _cmSAT->solve(&_currentAssumptions);
  unsigned int rv = 0;
  if (solution == CMSat::l_True) {
    rv = 10;
  } else if (solution == CMSat::l_False) {
    rv = 20;
  } else {
    rv = 0;
  }
  return rv;
}

void CryptoMiniSATSolverProxy::Reset() {
  delete _cmSAT;
  _cmSAT = new CMSat::SATSolver;
  SetNumberOfThreads(_noOfThreads);
  NewVariable();
  ClearAssumption();
  ResetClause();
  //  _cmSAT->set_allow_otf_gauss();
}

void CryptoMiniSATSolverProxy::SetNumberOfThreads(unsigned int n) {
  _noOfThreads = n;
  _cmSAT->set_num_threads(n);
}

// void CryptoMiniSATSolverProxy::SetFrozen(int variable)
//{
////    _cmSAT->_setFrozen(variable, true);
//}

unsigned int CryptoMiniSATSolverProxy::GetNumberOfVariables() {
  return static_cast<unsigned int>(_cmSAT->nVars());
}

void CryptoMiniSATSolverProxy::NewVariables(unsigned number) {
  _cmSAT->new_vars(number);
}

void CryptoMiniSATSolverProxy::SetReconf(unsigned reconf) {
  // possible values
  // 0,3,4,6,7,12,13,14,15,16

  if (!(reconf == 0 || reconf == 3 || reconf == 4 || reconf == 6 ||
        reconf == 7 || reconf == 12 || reconf == 13 || reconf == 14 ||
        reconf == 15 || reconf == 16)) {
    std::cout << "Wrong value for reconf, has to be out of "
                 "0,3,4,6,7,12,13,14,15,16. Actual value is: "
              << reconf << std::endl;
    exit(1);
  }
  _cmSAT->set_reconf_at(0);
  _cmSAT->set_reconf(reconf);
}

// unsigned int CryptoMiniSATSolverProxy::GetNumberOfClauses()
//{
////    return static_cast<unsigned int>(_cmSAT->nClauses());
//}
