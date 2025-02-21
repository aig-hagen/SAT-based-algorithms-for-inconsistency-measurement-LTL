/********************************************************************************************
Copyright (c) 2014-2016, Sven Reimer,
Copyright (c) 2017-2020, Tobias Paxian

dPermission is hereby granted, free of charge, to any person obtaining a copy of
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

#include <assert.h>
#include <cmath>
#include <iomanip>

#include "Pacose.h"
#include "SATSolverProxy.h"
#include "Softclause.h"
#include "cascade.h"
#include "dgpw.h"
#include "multiplecascade.h"
#include "sorter.h"
#include "timemeasurement.h"
#include "timevariables.h"

#include <algorithm>
#include <chrono>
//#include <cmath>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

namespace DGPW {
DGPW::DGPW(Pacose *pacose)
    : _sumOfSoftWeights(),
      _softClauses(),
      _sorterTree(),
      _solver(pacose->_satSolver),
      _dgpwSetting(&pacose->_settings),
      _pacose(pacose),
      _maxSorterDepth(0),
      _maxsatResult(UNKNOWN),
      _optimum(nullptr),
      _timeVariables(nullptr),
      _mainCascade(nullptr),
      _mainMultipleCascade(nullptr),
      _greatestCommonDivisor(1),
      _currentBucketForTare(0),
      _satWeight(0),
      _moreThanTwoWeights(false),
      //    _topWeight(pacose->_top),
      _minWeight(1),
      _maxWeight(0),
      _resultUnknown(false),
      _clausesBefore(pacose->_nbClauses),
      _variablesBefore(pacose->_nbVars),
      _binaryClausesBefore(0),
      _ternaryClausesBefore(0),
      _addedClauses(0),
      _addedVariables(0),
      _addedBinaryClauses(0),
      _addedTernaryClauses(0),
      _binaryTopClauses(0),
      _ternaryTopClauses(0),
      _binaryBottomClauses(0),
      _ternaryBottomClauses(0),
      _solverCalls(0),
      _approxSolverCalls(0),
      _highestBucketMultiplicator(0),
      _currentBucketMultiplicator(0),
      _estimatedSatWeight(0),
      _diffToSatWeight(0),
      //    _collectedAssumptions(),
      _externalAssumptions({}),
      _fixedAssumptions({}),
      _lastResult(0),
      _previousPosition(1 - (_dgpwSetting->greedyPPSATPercentage *
                             0.01)),  //  - intitially 0% = 0.0
      _previousLowerBound(0.0),       // - intitially 100% = 1.0
      _noClauses(-1),  // number of clauses added, each clause forcing at least
                       // one SC to be true
      _noVars(-1),  // number of vars in one clause, one SC out of noVars has to
                    // be satisfied.
      _unsatisfiableSCs(0),
      _greedyPrepro(_dgpwSetting->greedyPrepro)

{
  _sorterTree.resize(1);
}

DGPW::~DGPW(void) {
  for (uint32_t i = 0; i != _sorterTree.size(); ++i) {
    for (uint32_t j = 0; j != _sorterTree[i].size(); ++j) {
      delete _sorterTree[i][j];
    }
  }
  //    for (uint32_t i = 0; i != _softClauses.size(); ++i)
  //      {
  //        delete _softClauses[i];
  //      }
  //    delete _dgpwSetting;
  // TODO: check, if _timeVariables, _mainCascade, and _mainMultipleCascade have
  // to be deleted here!
  //    delete _timeVariables;
}

/********** Interface to Glucose and reimplemented functions from AntomBase
 * **********/

// Creates and returns a fresh variable index not used so far.
uint32_t DGPW::NewVariable(void) {
  _addedVariables++;
  return static_cast<uint32_t>(_solver->NewVariable());
  //    Glucose::Var newVar = _glucose->newVar(true, true);
  //    return static_cast<uint32_t>(newVar) + 1;
}

// TODO: optimize -> change interface to Glucose data structures
bool DGPW::AddClause(std::vector<uint32_t> &clause, uint32_t lbd) {
  //      for(auto literal : clause) {
  //        _solver->SetFrozen(static_cast<int>(literal << 1));
  //      }
  _addedClauses++;
  return _solver->AddClause(clause);
  //    Glucose::vec<Glucose::Lit> newClause;
  //    for (uint32_t lit : clause )
  //      {
  //        assert(lit > 1);
  //        Glucose::Lit newLit = Glucose::toLit(static_cast<int>(lit - 2));
  //        newClause.push(newLit);
  //      }
  //    return _glucose->addClause_(newClause);
}

//  bool DGPW::AddClause(Glucose::vec<Glucose::Lit>& clause)
//  {
//    return _glucose->addClause_(clause);
//  }

unsigned DGPW::GetLastResult() { return _lastResult; }

bool DGPW::AddUnit(uint32_t lit) {
  std::vector<uint32_t> unitclause;
  unitclause.push_back(lit);
  //    std::cout << "                                              UNIT CLAUSE:
  //    " << lit << std::endl;

  return AddClause(unitclause);
}

uint32_t DGPW::Solve(void) {
  //      std::cout << "SOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOLVE()!" <<
  //      std::endl;
  _solver->ClearAssumption();

  if (!_fixedAssumptions.empty()) {
    _solver->AddAssumptions(_fixedAssumptions);
  }
  _lastResult = _solver->Solve();
  _solverCalls++;
  return _lastResult;
}

uint32_t DGPW::Solve(std::vector<uint32_t> &assumptions) {
  //      std::cout <<
  //      "SOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOLVE(assumptions)!" <<
  //      std::endl;
  _solver->ClearAssumption();
  if (!_fixedAssumptions.empty()) {
    _solver->AddAssumptions(_fixedAssumptions);
  }
  _solver->AddAssumptions(assumptions);
  _lastResult = _solver->Solve();
  _solverCalls++;
  return _lastResult;
  //    std::cout << "result: " << rv << std::endl;
  //    return rv;
}

// TODO: optimize!
// copies whole model vector
std::vector<uint32_t> DGPW::Model(void) const {
  //    std::cout << __PRETTY_FUNCTION__ << "  x  " << std::endl;

  //    _solver->ClearAssumptions();
  //    std::cout << _solver->Solve() << std::endl;

  std::vector<uint32_t> model;

  for (int32_t i = 0; i < static_cast<int32_t>(_solver->GetNumberOfVariables());
       ++i) {
    model.push_back(_solver->GetModel(i));
  }
  return model;
}

// TODO: Check this for correctness!
// TODO: Try to use this interface whenever possible
uint32_t DGPW::Model(uint32_t var) const {
  //      std::cout << __PRETTY_FUNCTION__ << "  y  " << std::endl;
  //    assert(var>0);
  //      std::cout << "_solver->GetLastResult(): " << _solver->GetLastResult()
  //      << std::endl;
  if (_lastResult != 10) return 0;

  assert(var <= Variables());
  unsigned value = _solver->GetModel(static_cast<int>(var));
  assert(value <= 2 * Variables() + 1);
  assert(value <= 2 * Variables() + 1);
  return value;
}

uint32_t DGPW::Variables() const { return _solver->GetNumberOfVariables(); }

uint32_t DGPW::Clauses() const { return _solver->GetNumberOfClauses(); }

uint32_t DGPW::StaticClauses() const {
  //    return static_cast<uint32_t>(_glucose->nClauses());
  return Clauses();
}

uint32_t DGPW::CurrentBinaryClauses() const {
  // TODO: implement me!
  return 0;
}
uint32_t DGPW::CurrentTernaryClauses() const {
  // TODO: implement me!
  return 0;
}

Settings *DGPW::GetSetting() { return _dgpwSetting; }

/********** End: Interface to Glucose and reimplemented functions from AntomBase
 * **********/

bool DGPW::AddWeightedSoftClause(std::vector<uint32_t> &clause,
                                 uint64_t weight) {
  SoftClause *sclaus = new SoftClause(0, clause, weight);
  _softClauses.push_back(sclaus);

  return true;
}

void DGPW::SetSoftClauseVector(std::vector<SoftClause *> *softClauses) {
  _softClauses = *softClauses;
}

// If sorter == nullptr, we want to count for all soft clauses
uint64_t DGPW::CountSatisfiedSoftClauses(Sorter *sorter,
                                         const std::vector<uint32_t> &model) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;

  std::vector<SoftClause *> softclauses;
  if (sorter == nullptr) {
    softclauses = _softClauses;
  } else {
    softclauses = sorter->GetSoftClauses();
  }

  uint64_t result = CountSatisfiedSoftClauses(softclauses, model);

  if (sorter != nullptr) {
    sorter->SetMinSatisfied(result);
  }
  return result;
}

uint64_t DGPW::CountSatisfiedSoftClauses(std::vector<SoftClause *> softclauses,
                                         const std::vector<uint32_t> &model) {
  if (_dgpwSetting->verbosity > 4) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << "SC.size: " << softclauses.size() << std::endl;
    //        std::cout << "Clauses.size: " << Clauses() << std::endl;
  }
  uint64_t result(0);
  // Proceed all soft clauses
  for (uint32_t i = 0; i != softclauses.size(); ++i) {
    uint32_t relaxlit = softclauses[i]->relaxationLit;

    //    if( model[relaxlit>>1] == relaxlit )
    if (Model(relaxlit >> 1) == relaxlit) {
      std::vector<uint32_t> clause(softclauses[i]->clause);
      uint32_t pos = 0;
      for (; pos != clause.size(); ++pos) {
        // clause satisfied without trigger?
        if (Model(clause[pos] >> 1) == clause[pos]) {
          result += softclauses[i]->weight;
          softclauses[i]->lastassignment = relaxlit ^ 1;
          //            std::cout << "++" << result;
          break;
        }
      }
      if (pos == clause.size()) {
        softclauses[i]->lastassignment = relaxlit;
      }
    } else if (Model(relaxlit >> 1) != 0) {
      assert(Model(relaxlit >> 1) == (relaxlit ^ 1));
      softclauses[i]->lastassignment = relaxlit ^ 1;
      result += softclauses[i]->weight;
      //        std::cout << " ++relax sat " << result;
    }
  }
  return result;
}

// Weigthed maxsat
// ____________________________________________________________________________________________________________________
uint32_t DGPW::MaxSolveWeightedPartial(int64_t &optimum) {
  //      _dgpwSetting->verbosity = 0;
  //      _dgpwSetting->solveAtFirst = true;
  //      // DGPW Standard settings!
  //      _dgpwSetting->encodeStrategy = ENCODEONLYIFNEEDED;
  //      _dgpwSetting->partitionStrategy = GROUPBYWEIGHT;
  //      _dgpwSetting->lastPos1 = true;
  //      _dgpwSetting->weightPlusOne = true;
  //      _dgpwSetting->encode01 = false;

  //      _dgpwSetting->solveAsWeighted = true;
  //      _dgpwSetting->analyze = true;

  //      _dgpwSetting->adderCaching = true;
  //      _dgpwSetting->coneOfInfluence = true;
  //      _dgpwSetting->exactBounding = true;
  //      _dgpwSetting->dGPW = true;
  //      _dgpwSetting->plainVariant = false;

  std::vector<uint32_t> externalAssumptions = {};
  uint32_t result = MaxSolveWeightedPartial(externalAssumptions, optimum);
  //  _timeVariables->DumpVariables(_dgpwSetting->currentCascade.iteration);
  //    std::cout << "result: " << result << std::endl;
  return result;
}

void DGPW::wbSortAndFilter(uint64_t UnSATWeight) {
  unsigned long noSCsBefore = _softClauses.size();

  //      for(unsigned i = 0; i < _softClauses.size(); i++)
  //      {
  //          if (_softClauses[i]->weight >= UnSATWeight) {
  //              // not satisfiable anymore - add negated unit clause of
  //              blockings AddUnit(_softClauses[i]->relaxationLit^1);
  //              _softClauses.erase(_softClauses.begin() + i);
  //          } /*else {
  //              // if weight is still fulfillable push it
  ////              tmpSoftClauses.push_back(_softClauses[i]);
  //          } */
  //      }

  //      _softClauses.assign(tmpSoftClauses.begin(), tmpSoftClauses.end());

  //    tmpSoftClauses.moveTo(_softClauses);
  if (noSCsBefore > _softClauses.size()) {
    printf("c A number of soft clauses remained = %zu\n", _softClauses.size());
  }
}

// ____________________________________________________________________________________________________________________
uint32_t DGPW::MaxSolveWeightedPartial(
    std::vector<uint32_t> &externalAssumptions, int64_t &optimum) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;
  //    _dgpwSetting->Print();

  // to change optimum value with another function called from cascade and
  // bucket
  _optimum = &optimum;
  //    _dgpwSetting->application = WEIGHTEDMAXSAT;

  _timeVariables = new TimeVariables();
  TimeMeasurement totalTime(&_timeVariables->total, true);
  //    struct rusage resources;
  //    getrusage(RUSAGE_SELF, &resources);
  //    double timeS = static_cast<double>(resources.ru_utime.tv_sec + 1.e-6 ) *
  //    static_cast<double>( resources.ru_utime.tv_usec ); double timeC = timeS;
  //    _control->SetStartTime(timeS);
  //    _control->SetSumTime(true);
  uint32_t currentresult(UNKNOWN);
  if (_dgpwSetting->checkIfSolutionIsUnique &&
      _dgpwSetting->currentCascade.iteration > 0) {
    if (_dgpwSetting->currentCascade._onlyWithAssumptions) {
      std::cout << "c BAD PARAMETER COMBINATION!!!" << std::endl;
      exit(1);
    }
    _dgpwSetting->currentCascade._onlyWithAssumptions = true;
  }

  // uint32_t averageBucketEntries(0);

  //#ifndef NDEBUG
  uint64_t softWeights = 0;
  std::for_each(_softClauses.begin(), _softClauses.end(),
                [&](SoftClause *SC) { softWeights += SC->weight; });
  _sumOfSoftWeights = softWeights;

  //	std::cout << "soft weights: " << softWeights << " sumofweights: " <<
  //_sumOfSoftWeights << std::endl;
  //    assert(softWeights == _sumOfSoftWeights);
  //#endif
  if (_sumOfSoftWeights == 1) {
  }

  //  std::cout << "c clauses before.........: " << _clausesBefore << std::endl;
  if (_dgpwSetting->verbosity > 0) {
    std::cout << "c #softclauses...........: " << _softClauses.size()
              << std::endl;
    std::cout << "c #sum of SoftWeights....: "
              << _sumOfSoftWeights * _greatestCommonDivisor << std::endl;
  }
  //    if (!_dgpwSetting->solveAtFirst)
  //      std::cout << "c #sum of SoftWeights....: " << _sumOfSoftWeights *
  //      _greatestCommonDivisor << std::endl;

  if (_dgpwSetting->onlyByTares ||
      _dgpwSetting->mcDivideStrategy != SOLVEINNORMALCASCADEMODE) {
    _dgpwSetting->encodeStrategy = ENCODEONLYIFNEEDED;
  }

  if (_dgpwSetting->solveAtFirst) {
    if (false) {
      _approxSolverCalls = _solverCalls;
      TimeMeasurement solvedFirst(&_timeVariables->solvedFirst);
      unsigned maxGreedyTime =
          static_cast<unsigned>(_dgpwSetting->greedyPPTimeLimit);
      unsigned maxRounds =
          static_cast<unsigned>(ceil(pow(_softClauses.size(), 0.3)) * 100);

      if (_softClauses.size() < _dgpwSetting->greedyMinSizeOfSet) {
        _greedyPrepro = 0;
      }

      if (_greedyPrepro != 1) {
        // old variant
        currentresult =
            GreedyMaxInitSATWeight(_greedyPrepro, maxGreedyTime, maxRounds);
      } else {
        // binary search
        currentresult =
            GreedyMaxInitSATWeightV2(_greedyPrepro, maxGreedyTime, maxRounds);
      }
      _approxSolverCalls = _solverCalls - _approxSolverCalls;
      _addedClauses = 0;
      _addedVariables = 0;
      if (!_fixedAssumptions.empty()) {
        return SATISFIABLE;
      }
      CalculateOverallOptimum(_satWeight, true);
      // now done in pacose
      // wbSortAndFilter(_sumOfSoftWeights - _satWeight);
    } else {
      currentresult = SATISFIABLE;
    }

    if (currentresult == UNSAT) {
      _satWeight = 0;
      std::cout << "c first Solver call is UNSAT!" << std::endl;
      std::cout << "UNSAT" << std::endl;
      if (!_dgpwSetting->currentCascade._onlyWithAssumptions) {
        return UNSAT;
      }
    } else if (currentresult == UNKNOWN) {
      std::cout << "c UNKNOWN!!!!!! shouldn't be!" << std::endl;
      return UNKNOWN;
    }

    if (_softClauses.size() == 1 && _satWeight == 0 &&
        _sumOfSoftWeights == _softClauses[0]->weight) {
      // try to solve the only remaining softclause!
      std::vector<uint32_t> assumptions;
      assumptions.push_back(_softClauses[0]->relaxationLit ^ 1);
      currentresult = Solve(assumptions);
      if (currentresult == UNSAT) {
        std::cout << "c The only SC is NOT satisfiable!" << std::endl;
        AddUnit(_softClauses[0]->relaxationLit);
        if (Solve() != 10) {
          std::cout << "Strange solver result, shouldn't be possible!"
                    << std::endl;
        }
        CalculateOverallOptimum(_satWeight, true);
        return SATISFIABLE;
      } else if (currentresult == SATISFIABLE) {
        std::cout << "c The only SC is SATISFIABLE!" << std::endl;
        AddUnit(_softClauses[0]->relaxationLit ^ 1);
        CalculateOverallOptimum(_satWeight, true);
        return SATISFIABLE;
      }
    }

    //    std::cout << "c SATWeight solved first.: " << _satWeight << std::endl;
    if (_softClauses.size() == 0) {
      return SATISFIABLE;
    }
  }

  if (_satWeight == _sumOfSoftWeights &&
      (!_dgpwSetting->currentCascade._onlyWithAssumptions ||
       _dgpwSetting->checkIfSolutionIsUnique)) {
    // add relaxation literals as UnitClauses!
    for (auto sc : _softClauses) {
      AddUnit(sc->relaxationLit ^ 1);
    }
    std::cout << "c All SoftClauses are Satisfiable!" << std::endl;
    return SATISFIABLE;
  }

  // calc and set all trivially conflicting softclauses
  // TOASK: maybe better seperately for each bucket - to know boundaries of
  // bucket! together with mode of solving bucket parts to get max
  // CheckAllWeightedConflictingSoftclauses();
  if (_dgpwSetting->mcDivideStrategy != SOLVEINNORMALCASCADEMODE) {
    _dgpwSetting->encodeStrategy = ENCODEONLYIFNEEDED;
    _mainMultipleCascade = new MultipleCascade(
        this, _dgpwSetting->onlyByTares, _dgpwSetting->tareCascadeOnlyByTares,
        _dgpwSetting->cascadeDivider, _dgpwSetting->maxBucketSize,
        _dgpwSetting->nOfCasc, _dgpwSetting->interimResult, _sumOfSoftWeights);
    //        std::cout << "Cascade Divider: " << _cascadeDivider << std::endl;
    if (_mainMultipleCascade->DivideAndConnectCascades(
            _dgpwSetting->mcDivideStrategy, &_softClauses)) {
      currentresult = _mainMultipleCascade->Solve();

      if (currentresult != UNKNOWN) currentresult = SATISFIABLE;
    }
  } else {
    _dgpwSetting->cascadeDivider = 0;
  }

  // no else if, because if _cascadeDivider returns only one cascade then
  // _cascadeDivider == 0
  if (_dgpwSetting->cascadeDivider == 0) {
    _clausesBefore = StaticClauses();
    _mainMultipleCascade = nullptr;
    //        std::cout << "_dgpwSetting->onlyByTares: " <<
    //        _dgpwSetting->onlyByTares << std::endl;
    _mainCascade = new Cascade(this, nullptr, _dgpwSetting->onlyByTares);

    //        std::cout << "PartitionStrategy: " <<
    //        _dgpwSetting->partitionStrategy << std::endl;

    _mainCascade->Fill(&_softClauses, _dgpwSetting->partitionStrategy,
                       _dgpwSetting->encodeStrategy);

    // if (_antomSetting->verbosity > 2 && _antomSetting->base > 2)
    //{
    // std::cout << "How often is Softclause in Bucket, because of base " <<
    // _baseMode << std::endl; DumpBucketStructure();
    //}

    if (_dgpwSetting->cascadeDivider == 0) {
      if (_dgpwSetting->encodeStrategy == ENCODEALL) {
        _clausesBefore = StaticClauses();
      }
      if (!_mainCascade->Encode()) {
        return UNKNOWN;
      }
      //      if (_dgpwSetting->encodeStrategy == ENCODEALL) {
      //        std::cout << "c #clauses of coding.....: "
      //                  << StaticClauses() - _clausesBefore << std::endl;
      //        std::cout << "c #variables of coding...: "
      //                  << Variables() - _variablesBefore << std::endl;
      //      }
    }

    // TOBI: SetDecisionStrategiesForCascadeModel!!!!
    // SetDecisionStrategiesForMaxSAT();
    currentresult =
        _mainCascade->Solve(_dgpwSetting->currentCascade._onlyWithAssumptions,
                            _dgpwSetting->currentCascade._solveTares);
  }

  //  if (_dgpwSetting->encodeStrategy != ENCODEALL &&
  //  _dgpwSetting->onlyByTares) {
  //    std::cout << "c #clauses of coding.....: " << _addedClauses <<
  //    std::endl; std::cout << "c #variables of coding...: " << _addedVariables
  //    << std::endl;
  //  }
  if (currentresult == UNKNOWN) {
    return UNKNOWN;
  }

  // if there is more than one cascade, we have to calculate
  // the OverallOptimum again to sum up the weights.
  //    CalculateOverallOptimum(_satWeight, true);
  optimum = optimum * _greatestCommonDivisor;

  if (_dgpwSetting->checkIfSolutionIsUnique &&
      _dgpwSetting->currentCascade.iteration > 0) {
    if (optimum == 0) {
      // special case 1 - all SCs are satisfied!
      for (auto SC : _softClauses) {
        AddUnit(SC->relaxationLit ^ 1);
      }
      //            currentresult = Solve();
      //            if (currentresult != SATISFIABLE) {
      //                std::cout << "c Strange result - shouldn't be possible!"
      //                << std::endl; exit(1);
      //            }
    } else if (_satWeight == 0) {
      // special case 2 - no SC could be satisfied
      for (auto SC : _softClauses) {
        AddUnit(SC->relaxationLit);
      }
      //            currentresult = Solve();
      //            if (currentresult != SATISFIABLE) {
      //                std::cout << "c Strange result - shouldn't be possible!"
      //                << std::endl; exit(1);
      //            }
    } else {
      std::vector<uint32_t> clause;
      for (auto SC : _softClauses) {
        clause.push_back(Model(SC->relaxationLit >> 1) ^ 1);
      }
      uint32_t relaxLit = NewVariable() << 1;
      clause.push_back(relaxLit);
      AddClause(clause);

      std::vector<uint32_t> assumptions = GetLastAssumptions();

      if (assumptions.size() == 0) {
        std::cout << "c NO ASSUMPTIONS - ERROR!!!" << std::endl;
        exit(1);
      }
      assumptions.push_back(relaxLit ^ 1);

      uint32_t currentresult = Solve(assumptions);

      if (currentresult != SATISFIABLE) {
        std::cout << "c SOLUTION IS UNIQUE! Add this solution as hard clauses "
                     "to the problem!"
                  << std::endl;
        AddUnit(relaxLit);
        clause.pop_back();
        // add solution as unit clauses
        for (auto lit : clause) {
          AddUnit(lit ^ 1);
        }
        // we don't need the encoding anymore - throw it away!!
        for (auto ass : assumptions) {
          AddUnit(ass ^ 1);
        }
        currentresult = Solve();
        if (currentresult != SATISFIABLE) {
          std::cout << "c Strange result - shouldn't be possible!" << std::endl;
          exit(1);
        }
      } else {
        std::cout << "c THERE ARE MORE THAN ONE SOLUTION!" << std::endl;
        assumptions.pop_back();
        AddUnit(relaxLit);
        for (auto ass : assumptions) {
          std::cout << ass << ", ";
          AddUnit(ass);
        }
      }
    }
  }

  if (_dgpwSetting->verbosity > 0) {
    if (!_dgpwSetting->formulaIsDivided) {
      std::cout << "s SATISFIABLE" << std::endl;
      //    std::cout << "o " << optimum << std::endl;
    } else {
      if (currentresult == SATISFIABLE) {
        std::cout << "c currently SAT" << std::endl;
        //      std::cout << "c local o " << optimum << std::endl;
      } else if (currentresult == 20) {
        std::cout << "c currently UNSAT" << std::endl;
        //      std::cout << "c local o " << optimum << std::endl;
      }
    }
  }

  if (_satWeight == _sumOfSoftWeights &&
      !_dgpwSetting->currentCascade._onlyWithAssumptions) {
    // add relaxation literals as UnitClauses!
    for (auto sc : _softClauses) {
      AddUnit(sc->relaxationLit ^ 1);
    }
    if (_dgpwSetting->verbosity > 0)
      std::cout << "c All SoftClauses are Satisfiable!" << std::endl;
    return SATISFIABLE;
  }
  if (_dgpwSetting->verbosity > 0)
    std::cout << "c number of variables: " << Variables() << std::endl;
  return currentresult;
}  // namespace DGPW

void DGPW::SetInitialAssumptions(std::vector<uint32_t> assumptions) {
  _externalAssumptions = assumptions;
}

unsigned DGPW::TryToSolveMoreThanOneWeight(
    std::vector<unsigned> &nextAssumptions, std::vector<SoftClause *> &UNSATSCs,
    unsigned round) {
  if (_dgpwSetting->verbosity > 5)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  unsigned maxSize = 0;

  if (round == 0) {
    // can we satisfy one of every divider elements by just dividing the SC
    // vector at every n'th element
    maxSize = 10;
  } else if (round == 1) {
    // can we satisfy at least 10 elements by just dividing the SC vector into
    // divider
    maxSize = static_cast<unsigned>(ceil(UNSATSCs.size() / 10));
  } else if (round == 2) {
    maxSize = static_cast<unsigned>(ceil(UNSATSCs.size() / 3));
  }

  unsigned currentresult = 0;
  std::vector<std::vector<unsigned>> clauses;
  std::vector<unsigned> clause;
  clauses.push_back(clause);

  for (unsigned i = 0; i < UNSATSCs.size(); ++i) {
    if (clauses.back().size() == maxSize && (UNSATSCs.size() - i >= maxSize)) {
      assert(clause.size() == 0);
      clauses.push_back(clause);
    }
    clauses.back().push_back(UNSATSCs[i]->relaxationLit ^ 1);
  }

  if (_dgpwSetting->verbosity > 1) {
    std::cout << "c maxSize = " << maxSize << std::endl;
    std::cout << "c UNSATSCs.size(): " << UNSATSCs.size() << std::endl;
    std::cout << "c number of added assumption clauses = " << clauses.size()
              << std::endl;
  }

  for (size_t index = 0; index < clauses.size(); ++index) {
    assert(clauses[index].size() > 1);
    // add relaxLit to clause and clause to the CNF, deactivate SC afterwards!
    clauses[index].push_back(NewVariable() << 1);
    // activate clause as assumption
    nextAssumptions.push_back(clauses[index].back() ^ 1);
    AddClause(clauses[index]);
  }

  currentresult = Solve(nextAssumptions);

  for (size_t index = clauses.size(); index > 0; --index) {
    //          std::cout << nextAssumptions.back() << " == " <<
    //          clauses[index-1].back() << std::endl; std::cout <<
    //          "clauses.back().back() " << clauses.back().back() << std::endl;
    assert(nextAssumptions.back() == (clauses[index - 1].back() ^ 1));
    AddUnit(clauses[index - 1].back());
    assert(!nextAssumptions.empty());
    nextAssumptions.pop_back();
  }
  assert(Solve(nextAssumptions) == 10);

  return currentresult;
}

unsigned DGPW::TryToSolveMoreThanOneWeight2(
    std::vector<unsigned> &nextAssumptions, std::vector<SoftClause *> &UNSATSCs,
    unsigned round) {
  if (_dgpwSetting->verbosity > 5)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  unsigned maxSize = _lastSATSizeTryToSolveMoreThanOneWeight;

  unsigned currentresult = 0;
  std::vector<std::vector<unsigned>> clauses;
  std::vector<unsigned> clause;
  clauses.push_back(clause);

  // generate randomly sorted vector of indices
  //    std::cout << "UNSATSCs.size(): " << UNSATSCs.size() << std::endl;
  std::vector<unsigned> indices(UNSATSCs.size());
  std::size_t n(0);
  std::generate(std::begin(indices),
                std::begin(indices) + static_cast<unsigned>(UNSATSCs.size()),
                [&] { return n++; });

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle(indices.begin(), indices.end(),
               std::default_random_engine(seed));

  for (unsigned i = 0; i < UNSATSCs.size(); ++i) {
    if (clauses.back().size() == maxSize &&
        (UNSATSCs.size() - i > (maxSize / 2))) {
      assert(clause.size() == 0);
      //              std::cout << "maxSize: " << maxSize << std::endl;
      assert(clauses.back().size() > 1);
      clauses.push_back(clause);
    }
    // randomly picked next element
    clauses.back().push_back(UNSATSCs[indices[i]]->relaxationLit ^ 1);
    //        next element in row
    //        clauses.back().push_back(UNSATSCs[i]->relaxationLit ^ 1);
  }

  if (_dgpwSetting->verbosity > 1) {
    std::cout << "c lastSATSizeTryToSolveMoreThanOneWeight: "
              << _lastSATSizeTryToSolveMoreThanOneWeight << std::endl;
    std::cout << "c maxSize = " << maxSize << std::endl;
    std::cout << "c UNSATSCs.size(): " << UNSATSCs.size() << std::endl;
    //        std::cout << "c number of added assumption clauses = " <<
    //        clauses.size() << std::endl;
  }

  std::cout << "c Force " << clauses.size() << " soft clauses to be satisfied!"
            << std::endl;

  for (size_t index = 0; index < clauses.size(); ++index) {
    assert(clauses[index].size() > 1);
    // add relaxLit to clause and clause to the CNF, deactivate SC afterwards!
    clauses[index].push_back(NewVariable() << 1);
    // activate clause as assumption
    nextAssumptions.push_back(clauses[index].back() ^ 1);
    AddClause(clauses[index]);
  }

  //      _solver->SetMaxPropagations(_preproPropagationLimit);
  currentresult = Solve(nextAssumptions);

  for (size_t index = clauses.size(); index > 0; --index) {
    //          std::cout << nextAssumptions.back() << " == " <<
    //          clauses[index-1].back() << std::endl; std::cout <<
    //          "clauses.back().back() " << clauses.back().back() << std::endl;
    assert(nextAssumptions.back() == (clauses[index - 1].back() ^ 1));
    AddUnit(clauses[index - 1].back());
    assert(!nextAssumptions.empty());
    nextAssumptions.pop_back();
  }
  //      assert(Solve(nextAssumptions) == 10);
  //    std::cout << "c _lastSATSizeTryToSolveMoreThanOneWeight: "
  //              << _lastSATSizeTryToSolveMoreThanOneWeight << std::endl;

  if (currentresult == SATISFIABLE) {
    std::cout << "c                 SAT" << std::endl;
    _lastSATSizeTryToSolveMoreThanOneWeight = maxSize;
    if (maxSize > 20) {
      _lastSATSizeTryToSolveMoreThanOneWeight =
          _lastSATSizeTryToSolveMoreThanOneWeight / 1.25;
    } else if (maxSize > 3) {
      _lastSATSizeTryToSolveMoreThanOneWeight--;
    }
  } else if (currentresult == UNSAT) {
    std::cout << "c                 UNSAT" << std::endl;
    _lastSATSizeTryToSolveMoreThanOneWeight =
        1.5 * _lastSATSizeTryToSolveMoreThanOneWeight;
    //                  std::cout << "Sizes: " << UNSATSCs.size() << " " <<
    //                  neverSATSCs.size() << std::endl;
  } else {
    std::cout << "c               UNKNOWN" << std::endl;
    _lastSATSizeTryToSolveMoreThanOneWeight =
        2 * _lastSATSizeTryToSolveMoreThanOneWeight;
  }

  return currentresult;
}

unsigned DGPW::BinarySearchSatisfySCs(std::vector<unsigned> &nextAssumptions,
                                      std::vector<SoftClause *> *unsatSCs) {
  if (_dgpwSetting->verbosity > 2) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    std::cout << "c greedy strategy: kind of binary search to force optimal "
                 "number of SC to be true!"
              << std::endl;
  }
  if (unsatSCs->empty()) {
    if (_dgpwSetting->verbosity > 1) {
      std::cout << "No more SCs to satisfy!" << std::endl;
    }
    return SATISFIABLE;
  }

  //  unsatSCs = &_softClauses;

  assert(_previousPosition <= 1);
  assert(_previousPosition > 0);
  if (_dgpwSetting->verbosity > 0) {
    std::cout << std::setw(30) << "_previousPosition: " << std::setw(10)
              << _previousPosition << std::endl;
    std::cout << std::setw(30) << "_previousLowerBound: " << std::setw(10)
              << _previousLowerBound << std::endl;
  }
  unsigned noUnsatSCs = unsatSCs->size();
  //  noUnsatSCs = 25;
  //  _previousPosition = 0.18;

  double sqrt = std::sqrt(noUnsatSCs);
  unsigned numberOfPositions = (floor(sqrt * 2)) - 1;
  if (noUnsatSCs == 2) {
    numberOfPositions++;
  }

  unsigned newPosition = 0;

  if (_dgpwSetting->verbosity > 5) {
    std::cout << std::setw(30) << "newPosition: " << std::setw(10)
              << newPosition << std::endl;

    std::cout << std::setw(30) << "noUnsatSCs: " << std::setw(10) << noUnsatSCs
              << std::endl;
    std::cout << std::setw(30) << "numberOfPositions: " << std::setw(10)
              << numberOfPositions << std::endl;
    std::cout << std::setw(30) << "sqrt: " << std::setw(10) << sqrt
              << std::endl;
    std::cout << "If we have n clauses, at least n SC have to be satisfied!"
              << std::endl;
    std::cout
        << "If we have n vars, at least one out of n SC have to be satisfied!"
        << std::endl;
    for (unsigned i = 1; i < numberOfPositions + 1; i++) {
      newPosition = i;
      if (newPosition < sqrt) {
        //        noClauses = ceil(noUnsatSCs / newPosition);
        _noClauses = noUnsatSCs / newPosition;
      } else {
        _noClauses = numberOfPositions - newPosition + 1;
      }
      _noVars = noUnsatSCs / _noClauses;
      std::cout << "Position, noClauses, noVars: " << std::setw(5) << i
                << std::setw(8) << _noClauses << std::setw(8)
                << std::round(_noVars * 100) / 100 << std::endl;
    }
  }

  newPosition = ceil(_previousPosition * numberOfPositions);
  if (newPosition < sqrt) {
    //    noClauses = ceil(noUnsatSCs / newPosition);
    assert(newPosition != 0);
    //    if (newPosition == 0) newPosition = 1;
    _noClauses = noUnsatSCs / newPosition;
  } else {
    _noClauses = numberOfPositions - newPosition + 1;
  }
  _noVars = (noUnsatSCs / _noClauses);
  if (_dgpwSetting->verbosity > 2) {
    std::cout << std::endl
              << "Position, noClauses, noVars: " << std::setw(5) << newPosition
              << std::setw(8) << _noClauses << std::setw(8)
              << std::round(_noVars * 100) / 100 << std::endl;
  }

  unsigned upperBound = noUnsatSCs % static_cast<int>(floor(_noVars));
  if (_dgpwSetting->verbosity > 3) {
    std::cout << "Number of clauses  X  Clause Size: " << std::setw(10)
              << upperBound << " X " << ceil(_noVars) << std::endl;
    std::cout << "Number of clauses  X  Clause Size: " << std::setw(10)
              << _noClauses - upperBound << "  X  " << floor(_noVars)
              << std::endl;
  }

  //  std::vector<unsigned> &nextAssumptions;

  //  unsigned round;
  //  exit(1);

  unsigned relaxlit = NewVariable() << 1;
  std::vector<unsigned> clause{relaxlit};
  std::vector<std::vector<unsigned>> clauses(_noClauses, clause);

  // generate randomly sorted vector of indices
  std::vector<unsigned> indices(unsatSCs->size());
  std::size_t n(0);
  std::generate(std::begin(indices),
                std::begin(indices) + static_cast<unsigned>(unsatSCs->size()),
                [&] { return n++; });

  //  unsigned seed =
  //  std::chrono::system_clock::now().time_since_epoch().count();
  unsigned seed = 1234567890;
  std::shuffle(indices.begin(), indices.end(),
               std::default_random_engine(seed));

  int varsPerClause = ceil(_noVars);
  //  std::cout << "unsatSCs->size(): " << unsatSCs->size()
  //            << " clauses.size(): " << clauses.size() << std::endl;

  int i = 0;
  for (unsigned j = 0; j < clauses.size(); j++) {
    if (j == upperBound) varsPerClause = floor(_noVars);
    for (int k = 0; k < varsPerClause; k++) {
      //      std::cout << "i: " << i << " varsPerClause: " << varsPerClause
      //                << " j: " << j << " UB: " << upperBound << "; ";
      if (_dgpwSetting->verbosity > 4)
        std::cout << ((*unsatSCs)[indices[i]]->relaxationLit ^ 1) << " ";
      //                << std::endl;
      clauses[j].push_back((*unsatSCs)[indices[i]]->relaxationLit ^ 1);
      i++;
    }
    AddClause(clauses[j]);
    if (_dgpwSetting->verbosity > 4) std::cout << std::endl;
  }

  assert(i == unsatSCs->size());
  //  unsigned nextAssumptions nextAssumptions.push_back(clauses[index].back() ^
  //  1);
  //    AddClause(clauses[index]);
  //  }

  //  //      _solver->SetMaxPropagations(_preproPropagationLimit);

  nextAssumptions.push_back(relaxlit ^ 1);

  unsigned currentresult = Solve(nextAssumptions);
  nextAssumptions.pop_back();

  //  double noClauses;  // number of clauses added, each clause forcing at
  //  least
  // one SC to be true
  //  double noVars;  // number of vars in one clause, one SC out of noVars has
  //  to
  // be satisfied.

  //  greedyPPSATPercentage;
  //  greedyPPUNSATPercentage;
  //  greedyMinSizeOfSet;

  if (currentresult == SATISFIABLE) {
    if (_dgpwSetting->verbosity > 0) {
      std::cout << std::setw(100) << "SAT" << std::endl;
      if (_noVars == 1) {
        std::cout << std::setw(100) << "c ALL REMAINING SCs COULD BE SATISFIED!"
                  << std::endl;
      }
    }
    _previousPosition -= (_previousPosition - _previousLowerBound) /
                         (1 / (_dgpwSetting->greedyPPSATPercentage * 0.01));
    if (_previousPosition == 0) _previousPosition = 0.000001;
    return currentresult;
  } else if (currentresult == UNSAT) {
    if (_dgpwSetting->verbosity > 0)
      std::cout << std::setw(100) << "UNSAT" << std::endl;
    _previousLowerBound = _previousPosition;
    _previousPosition += (1 - _previousLowerBound) /
                         (1 / (_dgpwSetting->greedyPPUNSATPercentage * 0.01));
    if (_previousPosition > 1) _previousPosition = 1;

    // no more soft clauses can be satisfied
    if (_noClauses == 1) {
      if (_dgpwSetting->verbosity > 0)
        std::cout
            << std::setw(100)
            << "c NO ADDITIONAL SOFTCLAUSES FROM THIS SET ARE SATISFIABLE!"
            << std::endl;

      return currentresult;
    }

    return BinarySearchSatisfySCs(nextAssumptions, unsatSCs);
  } else {
    if (_dgpwSetting->verbosity > 0)
      std::cout << std::setw(100) << "STRANGE RESULT" << std::endl;
    return 0;
  }

  return 0;
}

unsigned DGPW::TestNextSoftClauses(std::vector<unsigned> &nextAssumptions,
                                   std::vector<SoftClause *> &UNSATSCs,
                                   bool isMaxSAT,
                                   bool tryToSolveMoreThanOneWeight) {
  if (_dgpwSetting->verbosity > 2)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  unsigned currentresult = 0;
  assert(_lastSATSizeTryToSolveMoreThanOneWeight >= 3);

  //    tryToSolveMoreThanOneWeight = true;
  // STRATEGY 2!!!
  //      std::cout << "c TTSMTOW= " << tryToSolveMoreThanOneWeight << "
  //      UNSATSCs: " << UNSATSCs.size() << "
  //      _lastSATSizeTryToSolveMoreThanOneWeight: " <<
  //      _lastSATSizeTryToSolveMoreThanOneWeight << std::endl;
  if (tryToSolveMoreThanOneWeight) {
    std::cout << "c strategy 2: solve more than one weight!" << std::endl;
    for (unsigned j = 0;; ++j) {
      if (3 * _lastSATSizeTryToSolveMoreThanOneWeight > UNSATSCs.size()) {
        j = 3;
        currentresult =
            TestIfNextHighestWeightsAreSAT(&nextAssumptions, &UNSATSCs, j);
        if (_lastSATSizeTryToSolveMoreThanOneWeight > 10) {
          _lastSATSizeTryToSolveMoreThanOneWeight /= 2;
        }
        if (_dgpwSetting->verbosity > 1) {
          if (currentresult == SATISFIABLE) {
            std::cout << "c                   SAT" << std::endl;
          } else if (currentresult == UNSAT) {
            std::cout << "c                 UNSAT" << std::endl;
          } else {
            std::cout << "c               UNKNOWN" << std::endl;
          }
        }
        return currentresult;
      }

      currentresult =
          TryToSolveMoreThanOneWeight2(nextAssumptions, UNSATSCs, j);

      if (currentresult == SATISFIABLE) {
        return currentresult;
      }
    }
  }

  // STRATEGY 1!!!
  std::cout << "c strategy 1: at least one more weight has to be solved!!"
            << std::endl;
  for (unsigned j = 0; j < 3; ++j) {
    if (isMaxSAT) {
      // is maxSAT instance
      j = 3;
    }

    // test if one of the next (to - from) highest relaxLits can be satisfied,
    // if next assumption is satisfied too
    currentresult =
        TestIfNextHighestWeightsAreSAT(&nextAssumptions, &UNSATSCs, j);

    if (currentresult == SATISFIABLE) {
      if (_dgpwSetting->verbosity > 1)
        std::cout << "c                   SAT" << std::endl;
      break;
    } else if (currentresult == UNSAT) {
      if (_dgpwSetting->verbosity > 1)
        std::cout << "c                 UNSAT" << std::endl;
      //                  std::cout << "Sizes: " << UNSATSCs.size() << " " <<
      //                  neverSATSCs.size() << std::endl;
    } else {
      if (_dgpwSetting->verbosity > 1)
        std::cout << "c               UNKNOWN" << std::endl;
    }
    if (isMaxSAT) {
      // is pure maxSAT instance - test only once if one of all weights weight
      // is satisfiable!
      break;
    }
  }
  //      std::cout << "RETURN" << std::endl;
  return currentresult;
}

unsigned DGPW::GreedyMaxInitSATWeightV2(int greedyPrepro, unsigned maxTime,
                                        unsigned maxRounds) {
  uint64_t softWeights = 0;
  std::for_each(_softClauses.begin(), _softClauses.end(),
                [&](SoftClause *SC) { softWeights += SC->weight; });
  _sumOfSoftWeights = softWeights;

  std::cout << "c _sumOfSoftWeights : " << _sumOfSoftWeights << std::endl;

  if (_dgpwSetting->verbosity > 2)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  std::cout << "c maxTime: " << maxTime << std::endl;
  std::cout << "c maxRounds: " << maxRounds << std::endl;

  TimeMeasurement timeSolvedFirst(&_timeVariables->solvedFirst, true);

  unsigned currentresult = 0;
  std::vector<unsigned> nextAssumptions;
  std::vector<unsigned> highestSATAssignment;
  std::vector<SoftClause *> UNSATSCs = {};
  std::vector<SoftClause *> neverSATSCs = {};
  uint64_t actualSATWeight = 0;
  uint64_t firstSATWeight = 0;
  std::vector<unsigned> sortedSCIndices(_softClauses.size());
  unsigned localMax = 0;
  unsigned previousMax = 0;
  uint64_t unsatisfiableSCWeight = 0;

  unsigned lastNewMaxRound = 0;
  bool isMaxSat = false;

  bool addedUnitClausesWithoutGettingNewMax = false;
  unsigned onlyOneSCSatByChance = 0;

  double timeLimit = 5 * pow(_softClauses.size(), 0.35);
  _preproPropagationLimit = static_cast<unsigned long>(timeLimit * 2500000);

  double lastNewMaxTime = 0;
  double lastSCDeletedTime = 0;
  int64_t opti = -1;

  _lastSATSizeTryToSolveMoreThanOneWeight = 4;

  // to get the right order of weights to test them!
  std::size_t n(0);
  std::generate(
      std::begin(sortedSCIndices),
      std::begin(sortedSCIndices) + static_cast<unsigned>(_softClauses.size()),
      [&] { return n++; });
  std::stable_sort(
      std::begin(sortedSCIndices), std::end(sortedSCIndices),
      [&](std::size_t i1, std::size_t i2) {
        return (_softClauses[i2]->weight > _softClauses[i1]->weight);
      });
  if (_softClauses[sortedSCIndices[0]]->weight ==
      _softClauses[sortedSCIndices.back()]->weight) {
    isMaxSat = true;
    std::cout << "c it is unweighted MaxSAT!!!" << std::endl;
  }

  // first solve, with external assumptions, if given
  // could be the last satisfying model of another DGPW
  // for DivideDGPW mode!
  currentresult = Solve(_externalAssumptions);
  if (greedyPrepro == 0) {
    if (currentresult == SATISFIABLE) {
      CalculateOverallOptimum(0, true);
    }
    return currentresult;
  }

  std::cout << "c PrePro TimeLimit: " << timeLimit << std::endl;
  std::cout << "c PrePro propagation limit: " << _preproPropagationLimit
            << std::endl;

  if (_dgpwSetting->verbosity > 1)
    std::cout << "c start solving model!" << std::endl;
  for (unsigned rounds = 0; rounds < maxRounds; ++rounds) {
    // informations...
    if (_dgpwSetting->verbosity > 0) {
      std::cout << "c time of greedy prepro: "
                << timeSolvedFirst.CurrentTimeDiff() << std::endl;
      if (_dgpwSetting->verbosity > 4)
        std::cout << "c local Max: " << localMax
                  << "  rounds / 2: " << rounds / 2
                  << "  rounds - lastNewMaxRound: " << rounds - lastNewMaxRound
                  << "  ceil(maxRounds/40): " << ceil(maxRounds / 40)
                  << std::endl;
      if (_dgpwSetting->verbosity > 1)
        std::cout << "c rounds: " << rounds << std::endl;
    }

    //        std::vector<unsigned int> indices;
    //        for (unsigned int i = 0; i < UNSATSCs.size(); i++) {
    //            indices.push_back(i);
    //        }

    assert(currentresult == SATISFIABLE);
    size_t unsatSizeBefore = UNSATSCs.size();
    // get Info about the actual SC Model, the return values are sorted
    // highest weight first!
    tie(nextAssumptions, UNSATSCs, actualSATWeight) =
        SatisfiedSCsInfo(&sortedSCIndices);

    if (rounds == 0) {
      //      std::cout << "firstSATWeight: " << actualSATWeight << std::endl;
      //      std::cout << "_SATWeight: " << _satWeight << std::endl;
      firstSATWeight = actualSATWeight;
    }

    if (_softClauses.size() == UNSATSCs.size() - 1) {
      if (_dgpwSetting->verbosity > 1) {
        std::cout << "c clauses satisfied by chance: "
                  << unsatSizeBefore - UNSATSCs.size() << std::endl;
        std::cout << "c !!!!!!!!!!!!!!!!!!!!!!!!THIS SC CAN ONLY BE SATISFIED "
                     "ALONE! IF "
                     "WEIGHT OF THIS SC < ACTUAL FOUND MAX, REMOVE IT!!!"
                  << std::endl;
      }
      onlyOneSCSatByChance++;
    }

    //        std::cout << "actualSATWeight: " << actualSATWeight << "
    //        _satWeight: " << _satWeight
    //                  << std::endl;

    if (actualSATWeight > _satWeight) {
      _satWeight = actualSATWeight;
      opti = static_cast<int64_t>(_sumOfSoftWeights - _satWeight);
      addedUnitClausesWithoutGettingNewMax = false;
      _pacose->CalculateSATWeight();
      highestSATAssignment = GetLastSatisfiableAssignment();
      if (_dgpwSetting->verbosity > 0) {
        std::cout << "c                                        NEW MAX FOUND: "
                  << _satWeight * _greatestCommonDivisor << std::endl;
        std::cout << "c local o " << opti * _greatestCommonDivisor << std::endl;
      }
      lastNewMaxTime = timeSolvedFirst.CurrentTimeDiff();
      lastNewMaxRound = rounds;
      if (_satWeight == _sumOfSoftWeights) {
        if (_dgpwSetting->verbosity > 0)
          std::cout << "c All weights could be satisfied!!!" << std::endl;
        break;
      }
    }

    if ((timeLimit < timeSolvedFirst.CurrentTimeDiff() - lastNewMaxTime &&
         timeLimit < timeSolvedFirst.CurrentTimeDiff() - lastSCDeletedTime) ||
        (timeSolvedFirst.CurrentTimeDiff() > maxTime)) {
      if (_dgpwSetting->verbosity > 0)
        std::cout << "c timeout in greedy preprocessor!" << std::endl;
      break;
    }

    //    currentresult = TestNextSoftClauses(nextAssumptions, UNSATSCs,
    //    isMaxSat,
    //                                        tryToSolveMoreThanOneWeight);

    if (previousMax != localMax) {
      //      _previousPosition = 1.0;
      _previousPosition = (_dgpwSetting->greedyPPSATPercentage * 0.01);
      _previousLowerBound = 0.0;
    }
    if (_dgpwSetting->greedyPPFixSCs == 0 ||
        (_dgpwSetting->greedyPPFixSCs == 1 && localMax > 0)) {
      nextAssumptions.clear();
      neverSATSCs = BuildOrderedIntersection(&UNSATSCs, &neverSATSCs);
      currentresult = BinarySearchSatisfySCs(nextAssumptions, &neverSATSCs);
    } else {
      currentresult = BinarySearchSatisfySCs(nextAssumptions, &UNSATSCs);
    }
    previousMax = localMax;

    if (currentresult != SATISFIABLE) {
      localMax++;
      if (_dgpwSetting->verbosity > 1)
        std::cout << "c                                              The "
                  << localMax << " th local Max is found: " << actualSATWeight
                  << std::endl;

      neverSATSCs = BuildOrderedIntersection(&UNSATSCs, &neverSATSCs);
      if (neverSATSCs.empty()) {
        if (_dgpwSetting->verbosity > 1)
          std::cout << "c All clauses were at least once satisfied!"
                    << std::endl;
        break;
      }
      if (_dgpwSetting->verbosity > 1)
        std::cout << "c neverSATSCs.size(): " << neverSATSCs.size()
                  << std::endl;

      //      if (_noClauses != 1 || localMax == 1) {
      if ((_dgpwSetting->greedyPPFixSCs == 0 && _noClauses != 1) ||
          (_dgpwSetting->greedyPPFixSCs == 1 &&
           (localMax == 1 || _noClauses != 1))) {
        nextAssumptions.clear();
        currentresult = Solve(nextAssumptions);
        continue;
      }

      for (unsigned iter = 0; iter < neverSATSCs.size(); ++iter) {
        if (_dgpwSetting->greedyPPFixSCs == 2) {
          // no additional SCs can be satisfied!
          // all remaining SCs are unsat!

          nextAssumptions.clear();
          nextAssumptions.push_back(neverSATSCs[iter]->relaxationLit ^ 1);
          //                  _solver->DeactivateLimits();
          currentresult = Solve(nextAssumptions);
          if (currentresult == SATISFIABLE) {
            if (_dgpwSetting->verbosity > 0) {
              std::cout << "c                   SAT" << std::endl;
              std::cout << "c weight: " << neverSATSCs[iter]->weight
                        << "  relaxLit: " << neverSATSCs[iter]->relaxationLit
                        << std::endl;
            }
            break;
          }
        }
        // remove Softclause from set
        if (currentresult == UNSAT) {
          //                      std::cout << "c                   SAT" <<
          //                      std::endl;
          _unsatisfiableSCs++;
          unsatisfiableSCWeight += neverSATSCs[iter]->weight;

          if (!_fixedAssumptions.empty()) {
            _fixedAssumptions.push_back(neverSATSCs[iter]->relaxationLit);
            if (_dgpwSetting->verbosity > 0) {
              std::cout << "c Softclause with weight "
                        << neverSATSCs[iter]->weight
                        << " couldn't be satisfied! - Add fixed assumption!"
                        << std::endl;
            }
            continue;
          }
          if (_dgpwSetting->verbosity > 0) {
            std::cout << "c Softclause with weight "
                      << neverSATSCs[iter]->weight
                      << " couldn't be satisfied! - Remove this softclause! - "
                         "new SC.size: "
                      << _softClauses.size() - 1 << "  neverSATSCs.size() "
                      << neverSATSCs.size() << "  UNSATSCs.size() "
                      << neverSATSCs.size() << std::endl;
          }
          addedUnitClausesWithoutGettingNewMax = true;
          lastSCDeletedTime = timeSolvedFirst.CurrentTimeDiff();
          AddUnit(neverSATSCs[iter]->relaxationLit);

          for (unsigned j = 0; j < _softClauses.size(); ++j) {
            if (_softClauses[j]->relaxationLit ==
                neverSATSCs[iter]->relaxationLit) {
              _sumOfSoftWeights -= _softClauses[j]->weight;
              _softClauses.erase(_softClauses.begin() + j);
            }
          }

          neverSATSCs.erase(neverSATSCs.begin() + iter);

          // to get the right order of weights to test them!
          sortedSCIndices.pop_back();
          std::size_t n(0);
          std::generate(std::begin(sortedSCIndices),
                        std::begin(sortedSCIndices) +
                            static_cast<unsigned>(_softClauses.size()),
                        [&] { return n++; });
          std::stable_sort(
              std::begin(sortedSCIndices), std::end(sortedSCIndices),
              [&](std::size_t i1, std::size_t i2) {
                return (_softClauses[i2]->weight > _softClauses[i1]->weight);
              });

          iter--;
        }
        if (timeSolvedFirst.CurrentTimeDiff() > maxTime) {
          if (_dgpwSetting->verbosity > 0)
            std::cout << "c timeout in greedy preprocessor!" << std::endl;
          break;
        }
      }

      if (currentresult != SATISFIABLE) {
        if (_dgpwSetting->verbosity > 0)
          std::cout << "c no more softclauses can be satisfieda!" << std::endl;
        break;
      }
    }

    //    if (neverSATSCs.empty() && localMax > 0) {

    if (neverSATSCs.empty() &&
        (_dgpwSetting->greedyPPFixSCs == 0 ||
         (_dgpwSetting->greedyPPFixSCs != 1 && localMax > 0))) {
      if (_dgpwSetting->verbosity > 0)
        std::cout << "c no more softclauses are satisfiable!" << std::endl;

      break;
    }
  }

  //      for (auto rl : highestSATAssignment)
  //        std::cout << rl << " ";
  //      std::cout << std::endl;
  if (_optimum != nullptr) {
    *_optimum = opti;
  }

  if (_unsatisfiableSCs > 0) {
    std::cout << "c unsatisfiable SCs......: " << _unsatisfiableSCs
              << std::endl;
    std::cout << "c remaining SCs..........: " << _softClauses.size()
              << std::endl;
    std::cout << "c weight of unsat SCs....: " << unsatisfiableSCWeight
              << std::endl;
  }
  if (onlyOneSCSatByChance > 0)
    std::cout << "c only 1 SC SAT by chance: " << onlyOneSCSatByChance
              << std::endl;

  //      _solver->DeactivateLimits();
  // TODO!! maybe faster without this line!
  if (addedUnitClausesWithoutGettingNewMax) {
    currentresult = Solve(highestSATAssignment);
    if (currentresult != 10) {
      std::cout << "c Strange result UNSAT after greedy Prepro, shouldn't be!"
                << std::endl;
      assert(false);
    }
    _pacose->CalculateSATWeight();
  }

  std::cout << "c max SatWeight Prepro...: " << _satWeight << std::endl;
  std::cout << "c add to 1st solved......: " << _satWeight - firstSATWeight
            << std::endl;
  std::cout << "c time greedyPrepro......: "
            << timeSolvedFirst.CurrentTimeDiff() << std::endl;

  return 10;
  //      return SATISFIABLE;
}

unsigned DGPW::GreedyMaxInitSATWeight(int greedyPrepro, unsigned maxTime,
                                      unsigned maxRounds) {
  uint64_t softWeights = 0;
  std::for_each(_softClauses.begin(), _softClauses.end(),
                [&](SoftClause *SC) { softWeights += SC->weight; });
  _sumOfSoftWeights = softWeights;

  std::cout << "c _sumOfSoftWeights : " << _sumOfSoftWeights << std::endl;

  if (_dgpwSetting->verbosity > 2)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  std::cout << "c maxTime: " << maxTime << std::endl;
  std::cout << "c maxRounds: " << maxRounds << std::endl;

  TimeMeasurement timeSolvedFirst(&_timeVariables->solvedFirst, true);

  unsigned currentresult = 0;
  std::vector<unsigned> nextAssumptions;
  std::vector<unsigned> highestSATAssignment;
  std::vector<SoftClause *> UNSATSCs = {};
  std::vector<SoftClause *> neverSATSCs = {};
  uint64_t actualSATWeight = 0;
  std::vector<unsigned> sortedSCIndices(_softClauses.size());
  unsigned localMax = 0;
  uint64_t unsatisfiableSCWeight = 0;
  unsigned zeroSuccByChanceSATClauses = 0;
  unsigned lastNewMaxRound = 0;
  bool isMaxSat = false;
  bool tryToSolveMoreThanOneWeight = false;
  bool addedUnitClausesWithoutGettingNewMax = false;
  unsigned onlyOneSCSatByChance = 0;

  double timeLimit = 5 * pow(_softClauses.size(), 0.35);
  _preproPropagationLimit = static_cast<unsigned long>(timeLimit * 2500000);

  double lastNewMaxTime = 0;
  double lastSCDeletedTime = 0;
  int64_t opti = -1;

  _lastSATSizeTryToSolveMoreThanOneWeight = 4;

  // to get the right order of weights to test them!
  std::size_t n(0);
  std::generate(
      std::begin(sortedSCIndices),
      std::begin(sortedSCIndices) + static_cast<unsigned>(_softClauses.size()),
      [&] { return n++; });
  std::stable_sort(
      std::begin(sortedSCIndices), std::end(sortedSCIndices),
      [&](std::size_t i1, std::size_t i2) {
        return (_softClauses[i2]->weight > _softClauses[i1]->weight);
      });
  if (_softClauses[sortedSCIndices[0]]->weight ==
      _softClauses[sortedSCIndices.back()]->weight) {
    isMaxSat = true;
    std::cout << "c it is unweighted MaxSAT!!!" << std::endl;
  }

  // first solve, with external assumptions, if given
  // could be the last satisfying model of another DGPW
  // for DivideDGPW mode!
  currentresult = Solve(_externalAssumptions);
  if (greedyPrepro == 0) {
    if (currentresult == SATISFIABLE) {
      CalculateOverallOptimum(0, true);
    }
    return currentresult;
  }

  std::cout << "c PrePro TimeLimit: " << timeLimit << std::endl;
  std::cout << "c PrePro propagation limit: " << _preproPropagationLimit
            << std::endl;

  if (_dgpwSetting->verbosity > 1)
    std::cout << "c start solving model!" << std::endl;
  for (unsigned rounds = 0; rounds < maxRounds; ++rounds) {
    // informations...
    if (_dgpwSetting->verbosity > 0) {
      std::cout << "c time of greedy prepro: "
                << timeSolvedFirst.CurrentTimeDiff() << std::endl;
      if (_dgpwSetting->verbosity > 4)
        std::cout << "c local Max: " << localMax
                  << "  rounds / 2: " << rounds / 2
                  << "  rounds - lastNewMaxRound: " << rounds - lastNewMaxRound
                  << "  ceil(maxRounds/40): " << ceil(maxRounds / 40)
                  << std::endl;
      if (_dgpwSetting->verbosity > 1)
        std::cout << "c rounds: " << rounds << std::endl;
    }

    //        std::vector<unsigned int> indices;
    //        for (unsigned int i = 0; i < UNSATSCs.size(); i++) {
    //            indices.push_back(i);
    //        }

    assert(currentresult == SATISFIABLE);
    size_t unsatSizeBefore = UNSATSCs.size();
    // get Info about the actual SC Model, the return values are sorted
    // highest weight first!
    tie(nextAssumptions, UNSATSCs, actualSATWeight) =
        SatisfiedSCsInfo(&sortedSCIndices);

    if (unsatSizeBefore - 1 == UNSATSCs.size()) {
      zeroSuccByChanceSATClauses++;
      if (_dgpwSetting->verbosity > 0)
        std::cout << "c zeroSuccByChanceSATClauses: "
                  << zeroSuccByChanceSATClauses << std::endl;
      if (zeroSuccByChanceSATClauses > 15 || isMaxSat) {
        if (!tryToSolveMoreThanOneWeight &&
            static_cast<unsigned>(UNSATSCs.size()) / 5 > 3) {
          _lastSATSizeTryToSolveMoreThanOneWeight =
              static_cast<unsigned>(UNSATSCs.size()) / 5;
        }
        tryToSolveMoreThanOneWeight = true;
      }
    } else {
      zeroSuccByChanceSATClauses /= 1.5;
    }

    //        std::cout << "actualSATWeight: " << actualSATWeight << "
    //        _satWeight: " << _satWeight
    //                  << std::endl;

    if (actualSATWeight > _satWeight) {
      _satWeight = actualSATWeight;
      opti = static_cast<int64_t>(_sumOfSoftWeights - _satWeight);
      addedUnitClausesWithoutGettingNewMax = false;
      _pacose->CalculateSATWeight();
      highestSATAssignment = GetLastSatisfiableAssignment();
      if (_dgpwSetting->verbosity > 0) {
        std::cout << "c                                        NEW MAX FOUND: "
                  << _satWeight * _greatestCommonDivisor << std::endl;
        std::cout << "c local o " << opti * _greatestCommonDivisor << std::endl;
      }
      lastNewMaxTime = timeSolvedFirst.CurrentTimeDiff();
      lastNewMaxRound = rounds;
      if (_satWeight == _sumOfSoftWeights) {
        if (_dgpwSetting->verbosity > 0)
          std::cout << "c All weights could be satisfied!!!" << std::endl;
        break;
      }
    }

    if ((timeLimit < timeSolvedFirst.CurrentTimeDiff() - lastNewMaxTime &&
         timeLimit < timeSolvedFirst.CurrentTimeDiff() - lastSCDeletedTime) ||
        (timeSolvedFirst.CurrentTimeDiff() > maxTime)) {
      if (_dgpwSetting->verbosity > 0)
        std::cout << "c timeout in greedy preprocessor!" << std::endl;
      break;
    }

    if (zeroSuccByChanceSATClauses > 35 &&
        timeLimit / 5 < timeSolvedFirst.CurrentTimeDiff() - lastNewMaxTime &&
        timeLimit / 5 < timeSolvedFirst.CurrentTimeDiff() - lastSCDeletedTime) {
      if (_dgpwSetting->verbosity > 0)
        std::cout << "c Break greedy prepro too few by chance satisfied "
                     "softclauses: "
                  << zeroSuccByChanceSATClauses << " rounds." << std::endl;
      break;
    }

    currentresult = TestNextSoftClauses(nextAssumptions, UNSATSCs, isMaxSat,
                                        tryToSolveMoreThanOneWeight);

    if (currentresult != SATISFIABLE) {
      if (!isMaxSat) {
        tryToSolveMoreThanOneWeight = false;
      }
      localMax++;
      if (_dgpwSetting->verbosity > 1)
        std::cout << "c                                              The "
                  << localMax << " th local Max is found: " << actualSATWeight
                  << std::endl;
      neverSATSCs = BuildOrderedIntersection(&UNSATSCs, &neverSATSCs);
      if (neverSATSCs.size() == 0) {
        if (_dgpwSetting->verbosity > 1)
          std::cout << "c All clauses were at least once satisfied!"
                    << std::endl;
        break;
      }
      //              for ( auto SC : neverSATSCs ) {

      assert(Solve() == SATISFIABLE);
      for (unsigned iter = 0; iter < neverSATSCs.size(); ++iter) {
        nextAssumptions.clear();
        nextAssumptions.push_back(neverSATSCs[iter]->relaxationLit ^ 1);
        //                  _solver->DeactivateLimits();
        currentresult = Solve(nextAssumptions);
        if (currentresult == SATISFIABLE) {
          if (_dgpwSetting->verbosity > 0) {
            std::cout << "c                   SAT" << std::endl;
            std::cout << "c weight: " << neverSATSCs[iter]->weight
                      << "  relaxLit: " << neverSATSCs[iter]->relaxationLit
                      << std::endl;
          }
          break;
        }
        if (currentresult == UNSAT) {
          //                      std::cout << "c                   SAT" <<
          //                      std::endl;
          _unsatisfiableSCs++;
          unsatisfiableSCWeight += neverSATSCs[iter]->weight;

          if (!_fixedAssumptions.empty()) {
            _fixedAssumptions.push_back(neverSATSCs[iter]->relaxationLit);
            if (_dgpwSetting->verbosity > 0) {
              std::cout << "c Softclause with weight "
                        << neverSATSCs[iter]->weight
                        << " couldn't be satisfied! - Add fixed assumption!"
                        << std::endl;
            }
            continue;
          }
          if (_dgpwSetting->verbosity > 0) {
            std::cout << "c Softclause with weight "
                      << neverSATSCs[iter]->weight
                      << " couldn't be satisfied! - Remove this softclause! - "
                         "new SC.size: "
                      << _softClauses.size() - 1 << "  neverSATSCs.size() "
                      << neverSATSCs.size() << "  UNSATSCs.size() "
                      << neverSATSCs.size() << std::endl;
          }
          addedUnitClausesWithoutGettingNewMax = true;
          lastSCDeletedTime = timeSolvedFirst.CurrentTimeDiff();
          AddUnit(neverSATSCs[iter]->relaxationLit);

          for (unsigned j = 0; j < _softClauses.size(); ++j) {
            if (_softClauses[j]->relaxationLit ==
                neverSATSCs[iter]->relaxationLit) {
              _sumOfSoftWeights -= _softClauses[j]->weight;
              _softClauses.erase(_softClauses.begin() + j);
            }
          }

          neverSATSCs.erase(neverSATSCs.begin() + iter);

          // to get the right order of weights to test them!
          sortedSCIndices.pop_back();
          std::size_t n(0);
          std::generate(std::begin(sortedSCIndices),
                        std::begin(sortedSCIndices) +
                            static_cast<unsigned>(_softClauses.size()),
                        [&] { return n++; });
          std::stable_sort(
              std::begin(sortedSCIndices), std::end(sortedSCIndices),
              [&](std::size_t i1, std::size_t i2) {
                return (_softClauses[i2]->weight > _softClauses[i1]->weight);
              });

          iter--;
        }
        if (timeSolvedFirst.CurrentTimeDiff() > maxTime) {
          if (_dgpwSetting->verbosity > 0)
            std::cout << "c timeout in greedy preprocessor!" << std::endl;
          break;
        }
      }
      if (currentresult != SATISFIABLE) {
        if (_dgpwSetting->verbosity > 0)
          std::cout << "c no more softclauses can be satisfied!" << std::endl;
        break;
      }
    }
  }

  //      for (auto rl : highestSATAssignment)
  //        std::cout << rl << " ";
  //      std::cout << std::endl;
  if (_optimum != nullptr) {
    *_optimum = opti;
  }

  if (_unsatisfiableSCs > 0) {
    std::cout << "c unsatisfiable SCs......: " << _unsatisfiableSCs
              << std::endl;
    std::cout << "c remaining SCs..........: " << _softClauses.size()
              << std::endl;
    std::cout << "c weight of unsat SCs....: " << unsatisfiableSCWeight
              << std::endl;
    std::cout << "c max SatWeight Prepro...: " << _satWeight << std::endl;
  }
  if (onlyOneSCSatByChance > 0)
    std::cout << "c only 1 SC SAT by chance: " << onlyOneSCSatByChance
              << std::endl;

  //      _solver->DeactivateLimits();
  // TODO!! maybe faster without this line!
  if (addedUnitClausesWithoutGettingNewMax) {
    currentresult = Solve(highestSATAssignment);
    if (currentresult != 10) {
      std::cout << "c Strange result UNSAT after greedy Prepro, shouldn't be!"
                << std::endl;
      assert(false);
    }
    _pacose->CalculateSATWeight();
  }
  std::cout << "c time greedyPrepro......: "
            << timeSolvedFirst.CurrentTimeDiff() << std::endl;

  return 10;
  //      return SATISFIABLE;
}  // namespace DGPW

std::vector<SoftClause *> DGPW::BuildOrderedIntersection(
    std::vector<SoftClause *> *UNSATSCs,
    std::vector<SoftClause *> *neverSATSCs) {
  if (_dgpwSetting->verbosity > 2)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  if ((*neverSATSCs).size() == 0) {
    return *UNSATSCs;
  }

  std::vector<SoftClause *> v_intersection;

  std::sort((*neverSATSCs).begin(), (*neverSATSCs).end(), SoftClause::ids);
  std::sort((*UNSATSCs).begin(), (*UNSATSCs).end(), SoftClause::ids);

  //              std::cout << "Build intersection" << std::endl;
  std::set_intersection((*neverSATSCs).begin(), (*neverSATSCs).end(),
                        (*UNSATSCs).begin(), (*UNSATSCs).end(),
                        std::back_inserter(v_intersection), SoftClause::ids);

  //  std::cout << "Intersection with size " << v_intersection.size() << ": "
  //            << " neverSatSCs.size(): " << (*neverSATSCs).size() <<
  //            std::endl;

  // assertion draussen, da die funktion 2x aufgerufen wurde
  //  assert(v_intersection.size() < (*neverSATSCs).size());

  //      std::cout << "v_intersection before ordering: " << std::endl;
  //      for(auto n : v_intersection)
  //          std::cout << n->weight << ' ';
  //      std::cout << std::endl << std::endl;

  std::sort(v_intersection.begin(), v_intersection.end(), SoftClause::bigger);

  //      std::cout << "v_intersection after ordering: " << std::endl;
  //      for(auto n : v_intersection)
  //          std::cout << n->weight << ' ';
  //      std::cout << std::endl << std::endl;

  return v_intersection;
}

unsigned DGPW::TestIfNextHighestWeightsAreSAT(
    std::vector<unsigned> *nextAssumptions, std::vector<SoftClause *> *UNSATSCs,
    unsigned round) {
  //      round = 3;
  unsigned from = 0;
  unsigned to = 1;
  if (round == 1) {
    from = 1;
    to = 10;
  } else if (round == 2) {
    from = 11;
    to = static_cast<unsigned>((*UNSATSCs).size());
  } else if (round == 3) {
    from = 1;
    to = static_cast<unsigned>((*UNSATSCs).size());
  }
  if (_dgpwSetting->verbosity > 2)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (to >= UNSATSCs->size()) {
    to = static_cast<unsigned>(UNSATSCs->size() - 1);
  }
  if (UNSATSCs->size() < from || from >= to) {
    return UNKNOWN;
  }
  assert(from < to);

  unsigned currentresult = 0;
  std::vector<unsigned> clause;

  //      std::cout << "UNSATSCs.size: " << UNSATSCs->size() << std::endl;
  //      std::cout << "current round: " << round << std::endl;

  for (unsigned i = from; i < to; ++i) {
    clause.push_back((*UNSATSCs)[i]->relaxationLit ^ 1);
  }

  if (clause.size() == 1) {
    nextAssumptions->push_back(clause[0]);
  } else {
    assert(clause.size() > 1);
    // add relaxLit to clause and clause to the CNF, deactivate SC afterwards!
    clause.push_back(NewVariable() << 1);
    // activate clause as assumption
    nextAssumptions->push_back(clause.back() ^ 1);
    AddClause(clause);
  }
  //      _solver->SetMaxPropagations(_preproPropagationLimit);
  currentresult = Solve(*nextAssumptions);

  if (clause.size() > 1) {
    // deactivate clause - with unit clause
    AddUnit(clause.back());
  }
  nextAssumptions->pop_back();

  return currentresult;
}

std::tuple<std::vector<unsigned>, std::vector<SoftClause *>, uint64_t>
DGPW::SatisfiedSCsInfo(std::vector<unsigned> *sortedSCIndices) {
  if (_dgpwSetting->verbosity > 2)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  std::vector<unsigned> nextAssumptions;
  std::vector<SoftClause *> UNSATSCs = {};
  uint64_t actualWeight = 0;
  uint64_t unsatWeight = 0;

  for (auto i : *sortedSCIndices) {
    // highest not yet looked at weight!
    SoftClause *SC = _softClauses[i];
    //          std::cout << "SCWEIGHT: " << SC->weight << std::endl;
    //          std::cout << "NEXT RELAXLIT: " << SC->relaxationLit <<
    //          std::endl; std::cout << "Number of vars: " << Variables() <<
    //          std::endl; std::cout << "NEXT MODEL OF RELAXLIT: " <<
    //          Model(SC->relaxationLit>>1) << std::endl;

    if (Model(SC->relaxationLit >> 1) == SC->relaxationLit) {
      //            std::cout << "RelaxLit is NOT satisfiable" << std::endl;
      for (auto lit : SC->clause) {
        if (Model(lit >> 1) == lit) {
          // clause is satisfied without relaxLit
          nextAssumptions.push_back(SC->relaxationLit ^ 1);
          SC->lastassignment = SC->relaxationLit ^ 1;
          actualWeight += SC->weight;
          break;
        }
      }
    } else {
      //            std::cout << "RelaxLit IS satisfiable" << std::endl;
      nextAssumptions.push_back(SC->relaxationLit ^ 1);
      SC->lastassignment = SC->relaxationLit ^ 1;
      actualWeight += SC->weight;
      //            std::cout << "RelaxLit IS satisfiable" << std::endl;
    }
    if (nextAssumptions.empty() ||
        nextAssumptions.back() != (SC->relaxationLit ^ 1)) {
      UNSATSCs.push_back(SC);
      SC->lastassignment = SC->relaxationLit;
      unsatWeight += SC->weight;
    }
  }
  //    std::cout << "DONE!" << std::endl;
  //    std::cout << unsatWeight << " " << actualWeight << " " <<
  //    _sumOfSoftWeights << std::endl;

  assert(unsatWeight + actualWeight == _sumOfSoftWeights);
  //    std::cout << "DONE!" << std::endl;
  return std::make_tuple(nextAssumptions, UNSATSCs, actualWeight);
}

std::vector<unsigned> DGPW::GetLastSatisfiableAssignment() {
  std::vector<unsigned> SCModel;
  for (auto SC : _softClauses) {
    SCModel.push_back(SC->lastassignment);
  }
  return SCModel;
}

StructureInfo DGPW::AnalyzeandConvertStructure() {
  std::cout << __FUNCTION__ << std::endl;
  StructureInfo structure = AnalyzeStructure();

  switch (structure) {
    case ISSAT:
      return ISSAT;
    case ISMAXSAT:
      if (_minWeight > 1) {
        DivideAllSoftClausesByFactor(_minWeight);
        std::cout << "c DivideByDivisormW: " << _minWeight << std::endl;
      }
      break;
    case CONVERTTOMAXSAT:
      ConvertFormulaToMaxSAT(_maxWeight);
      // TOASK: better - TOTELL
      // This line makes it faster with a high factor!!
      // But the value of _minWeight has to be multiplied to optimum!
      if (_minWeight > 1) {
        DivideAllSoftClausesByFactor(_minWeight);
        std::cout << "c DivideByDivisormW: " << _minWeight << std::endl;
      }

      structure = ISMAXSAT;
      break;
    case DIVWEIGHTSBYDIVISOR:
      DivideAllSoftClausesByFactor(_greatestCommonDivisor);
      structure = ISWEIGHTEDMAXSAT;
      break;
    case ISWEIGHTEDMAXSAT:
      break;
  }
  // Add information to SoftClauses
  // AddSoftClauseVector();
  return structure;
}

StructureInfo DGPW::AnalyzeStructure() {
  // assert(_sumOfSoftWeights > _topWeight);
  bool maxIsMin = (_maxWeight == _minWeight);
  bool withHC = Clauses() != 0;
  bool minIsSet = _minWeight != ((uint64_t)-1);
  bool maxIsSet = _maxWeight != 0;
  // is SAT formula
  bool onlyHC = !minIsSet && !maxIsSet && withHC;
  // is MaxSAT formula -- NOT PARTIAL
  bool onlySCwithOneWeight = !_moreThanTwoWeights && maxIsMin && !withHC;
  // is partial MaxSAT formula - Set _minWeight to one
  bool oneHCWeightAndSCWeight = maxIsMin && withHC;
  // if sum (_minWeight) < _maxWeight --> is partial MaxSAT formula
  bool onlySCwithTwoWeights =
      !_moreThanTwoWeights && !maxIsMin && minIsSet && maxIsSet && !withHC;
  // Calculate the greatest common divisor of all softclauseweights.
  uint64_t greatestCommonDivisor = _minWeight;

  for (uint32_t ind = 0; ind < _softClauses.size(); ++ind) {
    if (greatestCommonDivisor == 1) break;
    greatestCommonDivisor =
        GreatestCommonDivisor(greatestCommonDivisor, _softClauses[ind]->weight);
  }
  _greatestCommonDivisor = greatestCommonDivisor;

  if (_dgpwSetting->verbosity > 2) {
    std::cout << std::setw(30) << "_topWeight " << _topWeight << std::endl;
    std::cout << std::setw(30) << "_minWeight " << _minWeight << std::endl;
    std::cout << std::setw(30) << "_maxWeight " << _maxWeight << std::endl;
    std::cout << std::setw(30) << "_moreThanTwoWeights " << _moreThanTwoWeights
              << std::endl;
    std::cout << std::setw(30) << "maxIsMin " << maxIsMin << std::endl;
    std::cout << std::setw(30) << "withHC " << withHC << std::endl;
    std::cout << std::setw(30) << "minIsSet " << minIsSet << std::endl;
    std::cout << std::setw(30) << "maxIsSet " << maxIsSet << std::endl;
    std::cout << std::setw(30) << "onlyHC " << onlyHC << std::endl;
    std::cout << std::setw(30) << "onlySCwithOneWeight " << onlySCwithOneWeight
              << std::endl;
    std::cout << std::setw(30) << "oneHCWeightAndSCWeight "
              << oneHCWeightAndSCWeight << std::endl;
    std::cout << std::setw(30) << "onlySCwithTwoWeights "
              << onlySCwithTwoWeights << std::endl;
    std::cout << std::setw(30) << "_sumOfSoftWeights " << _sumOfSoftWeights
              << std::endl;
  }

  if (onlyHC) {
    return ISSAT;
  }
  if (onlySCwithTwoWeights) {
    uint32_t sumOfMinWeights(0);

    // the sum of all _minWeights has to be smaller than _maxWeight
    // then it can be solved like MaxSAT formula.
    for (uint32_t ind = 0; ind < _softClauses.size(); ++ind) {
      if (_softClauses[ind]->weight != _minWeight) continue;

      sumOfMinWeights += _softClauses[ind]->weight;
    }

    if (sumOfMinWeights < _maxWeight) {
      _sumOfSoftWeights = sumOfMinWeights;
      return CONVERTTOMAXSAT;
    }
  }
  if (oneHCWeightAndSCWeight || onlySCwithOneWeight || oneHCWeightAndSCWeight) {
    return ISMAXSAT;
  }
  //    if (greatestCommonDivisor > 1)
  //      {
  //        _greatestCommonDivisor = greatestCommonDivisor;
  //        std::cout << "c greatest common divisor: " <<
  //        greatestCommonDivisor
  //        << std::endl; return dgpw::DIVWEIGHTSBYDIVISOR;
  //      }
  return ISWEIGHTEDMAXSAT;
}

uint64_t DGPW::GreatestCommonDivisor(uint64_t a, uint64_t b) {
  uint64_t temp;
  while (b > 0) {
    temp = b;
    b = a % b;
    a = temp;
  }
  return a;
}

void DGPW::ConvertFormulaToMaxSAT(uint64_t maxWeight) {
  // configure all _maxWeights as hardclauses
  std::vector<SoftClause *> newSoftClauseVector;
  for (uint32_t ind = 0; ind < _softClauses.size(); ++ind) {
    if (_softClauses[ind]->weight != maxWeight) {
      newSoftClauseVector.push_back(_softClauses[ind]);
      continue;
    }

    AddClause(_softClauses[ind]->clause);
  }

  _softClauses = newSoftClauseVector;
}

void DGPW::DivideAllSoftClausesByFactor(uint64_t factor) {
  for (uint32_t ind = 0; ind < _softClauses.size(); ++ind) {
    _softClauses[ind]->weight /= factor;
  }
  _sumOfSoftWeights /= factor;
}

uint64_t DGPW::CalcGlobalOpt() { return _pacose->CalculateSATWeight(); }

uint64_t DGPW::CalculateOverallOptimum(uint64_t satWeight, bool countAgain) {
  if (_dgpwSetting->verbosity > 3)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  // struct rusage resources;
  // Cascade and Bucket counts the satisfied SC only for their SC's
  // for an overall result count again.
  if (countAgain) {
    // maybe _maxSatWeight - because it could be a local optima of this Solver
    // call
    satWeight = CountSatisfiedSoftClauses(nullptr, {});
  }

  if (*_optimum > static_cast<int64_t>(_sumOfSoftWeights - satWeight) ||
      *_optimum == -1) {
    /*
      if (_antomSetting->verbosity > 0)
      {
      getrusage(RUSAGE_SELF, &resources);
      double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double)
      resources.ru_utime.tv_usec; std::cout << "c " << (timeC -
      _control->GetStartTime()) << "s" << std::endl;
      }
    */
    // better actualize satWeight once at the end.
    _satWeight = satWeight;
    //    std::cout << "_sumOfSoftWeights: " << _sumOfSoftWeights
    //              << "  satWeight: " << satWeight << std::endl;
    *_optimum = static_cast<int64_t>(_sumOfSoftWeights - satWeight);

    if (!_dgpwSetting->formulaIsDivided) {
      std::cout << "o " << *_optimum * _greatestCommonDivisor << std::endl;
    } else {
      _pacose->CalculateSATWeight();
    }

    if (_dgpwSetting->verbosity > 0)
      std::cout << std::setw(50) << "Calculated Global SATWeight: " << satWeight
                << std::endl;

  } else if (_dgpwSetting->verbosity > 2) {
    /*
      if (_antomSetting->verbosity > 3)
      {
      getrusage(RUSAGE_SELF, &resources);
      double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double)
      resources.ru_utime.tv_usec;
      //std::cout << "c " << (timeC - _control->GetStartTime()) << "s" <<
      std::endl;
      }
    */
    if (!_dgpwSetting->formulaIsDivided)
      std::cout << "o " << *_optimum * _greatestCommonDivisor << std::endl;
    std::cout << std::setw(50) << "Calculated Global SATWeight: " << satWeight
              << std::endl;
  }

  return _satWeight;
}

void DGPW::SetMoreThanTwoWeights(bool val) { _moreThanTwoWeights = val; }

bool DGPW::GetHasMoreThanTwoWeights() { return _moreThanTwoWeights; }

void DGPW::SetTopWeight(uint64_t val) {
  assert(val > 0);
  _topWeight = val;
}

void DGPW::SetMinWeight(uint64_t val) {
  // assert( val < (uint64_t) - 1 );
  _minWeight = val;
}

uint64_t DGPW::GetMinWeight() { return _minWeight; }

void DGPW::SetMaxWeight(uint64_t val) {
  // assert( val > 0 );
  _maxWeight = val;
}

uint64_t DGPW::GetMaxWeight() { return _maxWeight; }

uint64_t DGPW::GetSATWeight() { return _satWeight; }

std::vector<std::pair<uint64_t, uint32_t>> DGPW::GetTareVector(
    uint64_t weightDiff) {
  if (_mainCascade == nullptr) return {};

  return _mainCascade->GetTareVector(weightDiff);
}

std::vector<std::pair<uint64_t, uint32_t>> DGPW::GetWatchdogs(
    uint64_t weightDiff) {
  if (_mainCascade == nullptr) {
    return {};
  }

  return _mainCascade->GetWatchdogs(weightDiff);
}

std::vector<unsigned> DGPW::GetLastAssumptions() {
  if (_mainCascade == nullptr) {
    return {};
  }
  return _mainCascade->GetLastAssumptions();
}

void DGPW::SetFixedAssumptions(std::vector<unsigned> fixedAssumptions) {
  _fixedAssumptions = fixedAssumptions;
  if (_dgpwSetting->verbosity > 0) std::cout << "c fixed added assumptions: ";
  for (auto fa : fixedAssumptions) {
    std::cout << fa << " ";
  }
  std::cout << std::endl;
}

void DGPW::RemoveFixedAssumptions() { _fixedAssumptions.clear(); }

void DGPW::SetHasHardClauses(bool val) { _hasHardClauses = val; }

void DGPW::SetSatWeight(uint64_t val) {
  _satWeight = val;
  //  _minWeight = val;
}

bool DGPW::GetHasHardClauses() { return _hasHardClauses; }

void DGPW::SetGreatestCommonDivisor(uint64_t val) {
  _greatestCommonDivisor = static_cast<int64_t>(val);
}
}  // namespace DGPW
