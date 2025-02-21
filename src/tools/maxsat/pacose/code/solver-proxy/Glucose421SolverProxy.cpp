/********************************************************************************************
Glucose421SolverProxy.cpp -- Copyright (c) 2020, Tobias Paxian

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

#include "Glucose421SolverProxy.h"
#include "SATSolverProxy.h"
#include "glucose421/core/Solver.h"
#include "iostream"
//#include "glucose421/Solver.h"

Glucose421SolverProxy::Glucose421SolverProxy()
    : _glucose421(new Glucose::Solver) {
  Reset();
}

Glucose421SolverProxy::~Glucose421SolverProxy() { delete _glucose421; }

SATSolverType Glucose421SolverProxy::GetSATSolverType() {
  return SATSolverType::GLUCOSE421;
}

unsigned int Glucose421SolverProxy::GetModel(int var) {
  int varLbool = toInt(_glucose421->modelValue(var));
  unsigned int rv;
  //    std::cout << "varLbool " << varLbool << std::endl;
  if (varLbool == 0) {
    rv = static_cast<unsigned int>(var << 1);
  } else if (varLbool == 1) {
    rv = static_cast<unsigned int>(var << 1 ^ 1);
  } else if (varLbool == 2) {
    std::cout << "Strange Value undef for variable: " << var << std::endl;
    //         undefined means, both values are possible!
    rv = static_cast<uint32_t>(var);
  } else {
    std::cout << "Strange value vor varLbool: " << varLbool << std::endl;
    std::cout << "Value var: " << var << std::endl;
    assert(false);
  }
  return rv;
}

int Glucose421SolverProxy::NewVariable() { return _glucose421->newVar(); }

void Glucose421SolverProxy::NewClause() { _currentClause.clear(); }

void Glucose421SolverProxy::AddLiteral(unsigned int *lit) {
  // TODO: Ask why toLit is called when SATSolverProxy already converts them
  //  std::cout << "AddLit: " << *lit << "  CC.size()" << _currentClause.size()
  //            << std::endl;
  _currentClause.push(Glucose::toLit(static_cast<int>(*lit)));
}

bool Glucose421SolverProxy::CommitClause() {
  //  for (int i = 0; i < _currentClause.size(); i++) {
  //    std::cout << Glucose::toInt(_currentClause[i]) << ", ";
  //  }
  //  std::cout << std::endl;

  return _glucose421->addClause(_currentClause);
}

void Glucose421SolverProxy::ResetClause() { _currentClause.clear(); }

bool Glucose421SolverProxy::AddClause(std::vector<unsigned int> &clause) {
  Glucose::vec<Glucose::Lit> tmpClause;
  for (unsigned int literal : clause) {
    tmpClause.push(Glucose::toLit(static_cast<int>(literal)));
  }
  //    std::cout << Glucose::toInt(literal) << ", ";
  //  }
  //  std::cout << std::endl;
  return _glucose421->addClause(tmpClause);
}

void Glucose421SolverProxy::AddAssumption(unsigned int *lit) {
  //  std::cout << Glucose::toInt(*lit) << ";; ";
  _currentAssumptions.push(Glucose::toLit(static_cast<int>(*lit)));
}

void Glucose421SolverProxy::ClearAssumption() { _currentAssumptions.clear(); }

unsigned int Glucose421SolverProxy::Solve() {
  EnableTimeLimit();
  EnableMemoryLimit();

  //  for (int i = 0; i < _currentAssumptions.size(); i++) {
  //    std::cout << Glucose::toInt(_currentAssumptions[i]) << "; ";
  //  }
  //  std::cout << std::endl;

  return _glucose421->solve(_currentAssumptions) ? 10 : 20;
}

unsigned int Glucose421SolverProxy::Simplify() {
  return _glucose421->simplify() ? 10 : 20;
}

unsigned int Glucose421SolverProxy::SolveLimited() {
  EnableTimeLimit();
  EnableMemoryLimit();
  int cr = toInt(_glucose421->solveLimited(_currentAssumptions));
  //  _glucose421->printIncrementalStats();
  if (cr == 0)
    return 10;
  else if (cr == 1)
    return 20;
  else if (cr == 2)
    return 0;
  else
    return 99;
}

void Glucose421SolverProxy::SetPropagationBudget(int64_t propagationBudget) {
  _glucose421->setPropBudget(propagationBudget);
}

void Glucose421SolverProxy::Reset() {
  delete _glucose421;
  _glucose421 = new Glucose::Solver;
  //    _glucose421 = new Glucose::Solver;
  NewVariable();
  ClearAssumption();
  ResetClause();
}

// void Glucose421SolverProxy::ResetCounter() {
//  _variableCounter = 0;
//  _clauseCounter = 0;
//}

// unsigned Glucose421SolverProxy::GetVariableCounter() {
//  return _variableCounter;
//}
// unsigned Glucose421SolverProxy::GetClauseCounter() { return _clauseCounter; }

// void Glucose421SolverProxy::SetFrozen(int variable)
//{
//    _glucose421->setFrozen(variable, true);
//}

unsigned int Glucose421SolverProxy::GetNumberOfVariables() {
  return static_cast<unsigned int>(_glucose421->nVars());
}

unsigned int Glucose421SolverProxy::GetNumberOfClauses() {
  return static_cast<unsigned int>(_glucose421->nClauses());
}
