/********************************************************************************************
MapleGlucoseSolverProxy.cpp -- Copyright (c) 2020, Tobias Paxian

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

#include "MapleGlucoseSolverProxy.h"
#include "SATSolverProxy.h"
#include "mapleLCMDistChronoBT/SimpSolver.h"
//#include "MapleLCMDistChronoBT/solver.h"
#include "mapleLCMDistChronoBT/core/SolverTypes.h"

MapleGlucoseSolverProxy::MapleGlucoseSolverProxy(void)
    : SATSolverProxy()
    , _maple(new Minisat::SimpSolver)
    ,
    //    _maple(new Minisat::Solver),
    _currentClause()
    , _currentAssumptions()
{
    Reset();
}

MapleGlucoseSolverProxy::~MapleGlucoseSolverProxy(void)
{
    delete _maple;
}

SATSolverType MapleGlucoseSolverProxy::GetSATSolverType()
{
    return SATSolverType::MAPLEGLUCOSE;
}

unsigned int MapleGlucoseSolverProxy::GetModel(int var)
{
    int varLbool = toInt(_maple->modelValue(var));
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
        assert(false);
    }
    return rv;
}

int MapleGlucoseSolverProxy::NewVariable()
{
    return _maple->newVar();
}

void MapleGlucoseSolverProxy::NewClause()
{
    assert(_currentClause.size() == 0);
}

void MapleGlucoseSolverProxy::AddLiteral(unsigned int *lit)
{
    _currentClause.push(Minisat::toLit(static_cast<int>(*lit)));
}

bool MapleGlucoseSolverProxy::CommitClause()
{
    return _maple->addClause(_currentClause);
}

void MapleGlucoseSolverProxy::ResetClause()
{
    _currentClause.clear();
}

bool MapleGlucoseSolverProxy::AddClause(std::vector<unsigned int> &clause)
{
    Minisat::vec<Minisat::Lit> tmpClause;
    for (unsigned int literal : clause) {
        tmpClause.push(Minisat::toLit(static_cast<int>(literal)));
    }
    return _maple->addClause(tmpClause);
}

void MapleGlucoseSolverProxy::AddAssumption(unsigned int *lit)
{
    _currentAssumptions.push(Minisat::toLit(static_cast<int>(*lit)));
}

void MapleGlucoseSolverProxy::ClearAssumption()
{
    _currentAssumptions.clear();
}

unsigned int MapleGlucoseSolverProxy::Solve()
{
    EnableTimeLimit();
    EnableMemoryLimit();

    return _maple->solve(_currentAssumptions) ? 10 : 20;
}

void MapleGlucoseSolverProxy::Reset()
{
    delete _maple;
    _maple = new Minisat::SimpSolver;
    //    _maple = new Minisat::Solver;
    NewVariable();
    ClearAssumption();
    ResetClause();
}

void MapleGlucoseSolverProxy::SetFrozen(int variable)
{
    // I'm not sure if needed! Because no iterative calls are possible!
    _maple->setFrozen(variable, true);
}

unsigned int MapleGlucoseSolverProxy::GetNumberOfVariables()
{
    return static_cast<unsigned int>(_maple->nVars());
}

unsigned int MapleGlucoseSolverProxy::GetNumberOfClauses()
{
    return static_cast<unsigned int>(_maple->nClauses());
}
