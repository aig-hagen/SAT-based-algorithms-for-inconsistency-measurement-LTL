/********************************************************************************************
AntomSolverProxy.cpp -- Copyright (c) 2020, Tobias Paxian

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


#include "AntomSolverProxy.h"
#include "SATSolverProxy.h"
#include "antom/antom.h"
#include "antom/antombase.h"

//#include <algorithm>

AntomSolverProxy::AntomSolverProxy()
    : SATSolverProxy()
    , _antom(new antom::Antom())
    , _currentClause()
    , _currentAssumptions()
    , _howOftenAddClauseCalled(0)
{
    // Skip first literal to assert static signals !
    Reset();
}

AntomSolverProxy::~AntomSolverProxy()
{
    delete _antom;
}

SATSolverType AntomSolverProxy::GetSATSolverType()
{
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    return SATSolverType::ANTOM;
}

void AntomSolverProxy::Reset()
{
    ClearAssumption();
    ResetClause();
    _antom->InstanceReset();
    NewClause();
    unsigned int uLit = static_cast<unsigned int>(NewVariable() * 2);
    AddLiteral(&uLit);
    CommitClause();
}

bool AntomSolverProxy::AddClause(std::vector<unsigned int> &clause)
{
    if (clause.size() == 1) {
        return _antom->AddClause(clause, 1);
    } else {
        return _antom->AddClause(clause);
    }
}

unsigned int AntomSolverProxy::Solve()
{
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    //    std::cout << _howOftenAddClauseCalled << " antom addClause is called." << std::endl;
    //    std::cout << GetNumberOfClauses() << " AntomClauseDBSize." << std::endl;
    //    std::cout << std::dec << std::endl;
    //    _antom->DumpCNF();
    //    assert(false);

    if (_cpuLimit > 0) {
        _antom->SetCPULimit(_cpuLimit);
    }
    if (_memoryLimit > 0) {
        _antom->SetMemoryLimit(_memoryLimit);
    }
    return _antom->Solve(_currentAssumptions);
}

unsigned int AntomSolverProxy::GetModel(int var)
{
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    //    return _antom->Model().at(static_cast<unsigned int>(var));
    unsigned int vari = static_cast<unsigned int>(abs(var));
    return _antom->Model()[vari];
}

int AntomSolverProxy::NewVariable()
{
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    return static_cast<int>(_antom->NewVariable());
}

void AntomSolverProxy::ResetClause()
{
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    _currentClause.clear();
}

void AntomSolverProxy::AddAssumption(unsigned int *lit)
{
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    _currentAssumptions.push_back(*lit);
}

void AntomSolverProxy::ClearAssumption()
{
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    _currentAssumptions.clear();
}

bool AntomSolverProxy::CommitClause()
{
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    _howOftenAddClauseCalled++;
    //    std::sort(_currentClause.begin(), _currentClause.end());
    //    std::cout << std::dec;
    //    for(unsigned int lit : _currentClause)
    //    {
    //        std::cout << lit << " ";
    //    }
    //    std::cout << std::endl;
    bool success = _antom->AddClause(_currentClause, 1);
    _currentClause.clear();
    return success;
}

void AntomSolverProxy::AddLiteral(unsigned int *lit)
{
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    _currentClause.push_back(*lit);
}

void AntomSolverProxy::NewClause()
{
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    //    ResetClause();
    assert(_currentClause.empty());
}

unsigned int AntomSolverProxy::GetNumberOfVariables()
{
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    return _antom->Variables();
}

unsigned int AntomSolverProxy::GetNumberOfClauses()
{
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    return _antom->Clauses();
}

void AntomSolverProxy::AddVariablePrio(unsigned int variable, unsigned int prio)
{
    _antom->setVariablePriority(variable, prio);
}
