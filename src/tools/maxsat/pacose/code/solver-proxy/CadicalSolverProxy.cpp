/********************************************************************************************
CadicalSolverProxy.cpp -- Copyright (c) 2020, Tobias Paxian

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

#include "CadicalSolverProxy.h"
#include <iostream>

CadicalSolverProxy::CadicalSolverProxy()
    : _cadical(new CaDiCaL::Solver), _vars(0), _noClauses(0), _hasVars(false) {
  Reset();
}

CadicalSolverProxy::~CadicalSolverProxy() { delete _cadical; }

SATSolverType CadicalSolverProxy::GetSATSolverType() {
  return SATSolverType::CADICAL;
}

unsigned int CadicalSolverProxy::GetModel(int var) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;
  //    std::cout << var << std::endl;
  unsigned int uLit;
  //    int value = _cadical->val(var);
  if (_model.empty()) {
    //        std::cout << "Model is empty" << std::endl;
    return 0;
  }
  int value = _model[static_cast<unsigned>(var)];
  if (value < 0) {
    uLit = (static_cast<unsigned int>(var) << 1) ^ 1;
  } else {
    uLit = static_cast<unsigned int>(var << 1);
  }

  //    std::cout << "value: " << value << "  uLit: " << uLit << std::endl;
  return uLit;
}

int CadicalSolverProxy::NewVariable() {
  //  return _cadical->max();
  _vars++;
  return static_cast<int>(_vars);
}

void CadicalSolverProxy::NewVariables(unsigned number) {
  //  _cadical->init(number);
  if (number > _vars) {
    _vars = number;
  }
  _cadical->reserve(static_cast<int>(number));
}

void CadicalSolverProxy::NewClause() {}

void CadicalSolverProxy::AddLiteral(int *lit) {
  //  if (_cadical->state() == CaDiCaL::State::INVALID) {
  //    std::cout << "Cannot Add Literal - invalid state!!! - abort solving!"
  //              << std::endl;
  //    exit(1);
  //  }
  _cadical->add(*lit);
  _hasVars = true;
  //  std::cout << *lit << std::endl;
}

void CadicalSolverProxy::AddLiteral(unsigned *lit) {
  int literal;
  if ((*lit) & 1) {
    literal = static_cast<int>(-(*lit >> 1));
  } else {
    literal = static_cast<int>(*lit >> 1);
  }
  AddLiteral(&literal);
}

bool CadicalSolverProxy::CommitClause() {
  //  _cadical->add(0);
  if (!_hasVars) {
    std::cout << "empty clause - do not add it!!" << std::endl;
    return true;
  }

  int lastLit = 0;
  AddLiteral(&lastLit);
  //    std::cout << "nextCl" << std::endl;
  _noClauses++;
  return true;
}

void CadicalSolverProxy::ResetClause() { _hasVars = false; }

void CadicalSolverProxy::AddAssumption(unsigned int *lit) {
  //    if (_assumption == 0) {
  //        _assumption = NewVariable();
  //    }
  //    AddLiteral(lit);
  //    AddLiteral(&_assumption);
  //    CommitClause();

  int literal;
  if ((*lit) & 1) {
    literal = static_cast<int>(-(*lit >> 1));
  } else {
    literal = static_cast<int>(*lit >> 1);
  }

  _assumptions.push_back(literal);
  //    std::cout << "ass.back: " << _assumptions.back() << std::endl;
}

void CadicalSolverProxy::ClearAssumption() {
  //    AddLiteral(&_assumption);
  //    CommitClause();
  //    _assumption = 0;
  _assumptions.clear();
}

unsigned int CadicalSolverProxy::Solve() {
  // TODO: Use in-built methods when available
  EnableTimeLimit();
  EnableMemoryLimit();
  //  if (_assumption != 0) {
  //    _cadical->assume(-_assumption);
  //  }
  //    static_cast<unsigned int>(_cadical->solve());
  for (auto assumption : _assumptions) {
    //    if (_cadical->state() != CaDiCaL::State::READY) {
    //      std::cout << "Cannot Add Assumption - invalid state!!! - abort
    //      solving!"
    //                << std::endl;
    //      //      exit(1);
    //    }
    _cadical->assume(assumption);
  }
  ClearAssumption();
  //  std::cout << "Vars: " << GetNumberOfVariables() << std::endl;
  //  std::cout << "Clauses: " << GetNumberOfClauses() << "/" << _noClauses
  //            << std::endl;

  //  if (_cadical->state() != CaDiCaL::State::READY) {
  //    std::cout << "Cannot Solve - invalid state!!! - abort solving!"
  //              << std::endl;
  //    //    exit(1);
  //  }
  int rv = _cadical->solve();
  if (rv == 10) {
    SaveWholeModel();
  } else {
    _model.clear();
  }

  //    return static_cast<unsigned int>(_cadical->solve());
  return static_cast<unsigned int>(rv);
}

void CadicalSolverProxy::SetFrozen(int variable) { _cadical->freeze(variable); }

void CadicalSolverProxy::MeltFrozen(int variable) { _cadical->melt(variable); }

void CadicalSolverProxy::Reset() {
  delete _cadical;
  _cadical = new CaDiCaL::Solver;
}

unsigned int CadicalSolverProxy::GetNumberOfVariables() {
  return static_cast<unsigned>(_cadical->vars());
}

unsigned int CadicalSolverProxy::GetNumberOfClauses() {
  return static_cast<unsigned>(_cadical->irredundant());
}

void CadicalSolverProxy::SaveWholeModel() {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;
  _model.clear();
  _model.push_back(0);

  for (int var = 1; var <= static_cast<int>(GetNumberOfVariables()); var++) {
    _model.push_back(_cadical->val(var));
  }
}
