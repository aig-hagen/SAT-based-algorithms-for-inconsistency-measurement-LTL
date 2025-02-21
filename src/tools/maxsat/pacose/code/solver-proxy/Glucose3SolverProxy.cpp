/********************************************************************************************
Glucose3SolverProxy.cpp -- Copyright (c) 2020, Tobias Paxian

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

#include "Glucose3SolverProxy.h"
#include <iostream>
#include "glucose3/core/Solver.h"

Glucose3SolverProxy::Glucose3SolverProxy()
    : SATSolverProxy(),
      _glucose3(new Glucose3::Solver()),
      _currentClause(),
      _currentAssumptions() {
  Reset();
}

Glucose3SolverProxy::~Glucose3SolverProxy() { delete _glucose3; }

SATSolverType Glucose3SolverProxy::GetSATSolverType() {
  return SATSolverType::GLUCOSE3;
}

unsigned int Glucose3SolverProxy::GetModel(int var) {
  int varLbool = toInt(_glucose3->modelValue(var));
  unsigned int rv;
  if (varLbool == 0) {
    rv = static_cast<unsigned int>(var << 1);
  } else if (varLbool == 1) {
    rv = static_cast<unsigned int>(var << 1 ^ 1);
  } else if (varLbool == 2) {
    std::cout << "Strange value undefined for variable: " << var << std::endl;
    rv = static_cast<uint32_t>(var);
  } else {
    std::cout << "Strange value for varLbool: " << varLbool << std::endl;
    assert(false);
  }
  return rv;
}

int Glucose3SolverProxy::NewVariable() { _glucose3->newVar(); }

void Glucose3SolverProxy::NewClause() {
  if (_currentClause.size() > 0) {
    exit(1);
  }
}

void Glucose3SolverProxy::AddLiteral(unsigned int *lit) {
  _currentClause.push(Glucose3::toLit(static_cast<int>(*lit)));
}

bool Glucose3SolverProxy::CommitClause() {
  return _glucose3->addClause(_currentClause);
}

void Glucose3SolverProxy::ResetClause() { _currentClause.clear(); }

bool Glucose3SolverProxy::AddClause(std::vector<unsigned int> &clause) {
  //    Glucose3::vec<Glucose3::Lit> tmpClause;
  ResetClause();
  for (unsigned int literal : clause) {
    //        tmpClause.push(Glucose3::toLit(static_cast<int>(literal)));
    AddLiteral(&literal);
  }
  bool rv = CommitClause();
  ResetClause();
  //    return _glucose3->addClause(tmpClause);
  return rv;
}

void Glucose3SolverProxy::AddAssumption(unsigned int *lit) {
  _currentAssumptions.push(Glucose3::toLit(static_cast<int>(*lit)));
}

void Glucose3SolverProxy::ClearAssumption() { _currentAssumptions.clear(); }

unsigned int Glucose3SolverProxy::Solve() {
  EnableTimeLimit();
  EnableMemoryLimit();

  //    std::cout << "Vars: " << _glucose3->nVars() << std::endl;
  //    std::cout << "Clauses: " << _glucose3->nClauses() << std::endl;
  unsigned int rv = _glucose3->solve(_currentAssumptions) ? 10 : 20;
  return rv;
}

void Glucose3SolverProxy::Reset() {
  // std::cout << "Going to delete _glucose3" << std::endl;
  delete _glucose3;
  // std::cout << "Going to make new SimpSolver pointer" << std::endl;
  _glucose3 = new Glucose3::Solver();
  // std::cout << "Going to call NewVariable" << std::endl;
  NewVariable();
  // std::cout << "Going to call ClearAssumption" << std::endl;
  ClearAssumption();
  // std::cout << "Going to call ResetClause" << std::endl;
  ResetClause();
  // std::cout << "Done" << std::endl;
}

// void Glucose3SolverProxy::SetFrozen(int variable)
//{
//    _glucose3->setFrozen(variable, true);
//}

unsigned int Glucose3SolverProxy::GetNumberOfVariables() {
  return static_cast<unsigned int>(_glucose3->nVars());
}

unsigned int Glucose3SolverProxy::GetNumberOfClauses() {
  return static_cast<unsigned int>(_glucose3->nClauses());
}
