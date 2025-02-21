/********************************************************************************************
Copyright (c) 2017-2020, Tobias Paxian

dPermission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************************************/

#include <assert.h>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <algorithm>
#include <iterator>
#include <random>
#include <vector>
#include "Greedyprepro.h"
#include "Pacose.h"
#include "SATSolverProxy.h"
#include "Settings.h"
#include "Softclause.h"
//#include "dgpw.h"  // for SATISFIABLE UNSATISFIABLE
#include <chrono>

#include "Greedyprepro.h"
#include "timemeasurement.h"
#include "timevariables.h"

#include <iostream>

#define UNKNOWN 0
#define SATISFIABLE 10
#define UNSAT 20

GreedyPrepro::GreedyPrepro(std::vector<SoftClause *> &softClauses,
                           Settings *settings, SATSolverProxy *solverProxy,
                           Pacose *pacose, int fixSCs)
    : _solver(solverProxy),
      _pacose(pacose),
      _settings(settings),
      _softClauses(softClauses),
      _satWeight(0),
      _unsatisfiableSCWeight(0),
      _sumOfSoftWeights(0),
      _unknownSolverCalls(0),
      _fixSoftClauses(fixSCs),
      _addedVariables(0),
      _addedClauses(0),
      _solverCalls(0),
      _unsatisfiableSCs(0),
      _satisfiableSCs(0),
      _previousPosition(1 - (_settings->greedyPPSATPercentage * 0.01)),
      _previousLowerBound(0.0),
      _noClauses(-1),
      _allWeightsSat(false),
      _pos1Unsat(false) {
  _timeVariables = new DGPW::TimeVariables();
  //  _solver = SATSolverProxy::InitSATSolver(SATSolverType::GLUCOSE421);
}

GreedyPrepro::~GreedyPrepro() {
  _solver->~SATSolverProxy();
  _timeVariables->~TimeVariables();
}

uint32_t GreedyPrepro::NewVariable(void) {
  _addedVariables++;
  return static_cast<uint32_t>(_solver->NewVariable());
}

// TODO: optimize -> change interface to Glucose data structures
bool GreedyPrepro::AddClause(std::vector<uint32_t> &clause, uint32_t lbd) {
  _addedClauses++;
  return _solver->AddClause(clause);
}

uint32_t GreedyPrepro::Solve(void) {
  _solver->ClearAssumption();

  _solverCalls++;
  return _solver->Solve();
}

uint32_t GreedyPrepro::Solve(std::vector<uint32_t> &assumptions) {
  _solver->ClearAssumption();
  _solver->AddAssumptions(assumptions);
  _solverCalls++;
  return _solver->Solve();
}

uint32_t GreedyPrepro::SolveLimited(std::vector<uint32_t> &assumptions) {
  _solver->ClearAssumption();
  _solver->AddAssumptions(assumptions);
  _solverCalls++;
  return _solver->SolveLimited();
}

uint32_t GreedyPrepro::Model(uint32_t var) const {
  unsigned value = _solver->GetModel(static_cast<int>(var));

  return value;
}

void GreedyPrepro::AddClausesToSATSolver(
    std::vector<std::vector<int>> &_clauseDB, unsigned nVars) {
  // print all hardClauses and the encoding if yet sent
  _solver->NewVariables(nVars + 1);
  //  std::cout << "NoVars: " << nVars << std::endl;
  //  unsigned var = _solver->NewVariable();
  //  std::cout << var << std::endl;
  if (_settings->verbosity > 0)
    std::cout << "c Adding now " << _clauseDB.size() << " hard clauses!"
              << std::endl;
  for (auto clause : _clauseDB) {
    _solver->ResetClause();
    _solver->NewClause();

    for (int lit : clause) {
      //      std::cout << lit << " ";
      //      unsigned int uLit;
      //      if (lit < 0) {
      //        uLit = static_cast<unsigned int>(((-lit) << 1) ^ 1);
      //      } else {
      //        uLit = static_cast<unsigned int>(lit << 1);
      //      }
      //      uLit = (static_cast<unsigned>(abs(lit)) << 1) ^ (lit < 0);
      //      std::cout << lit << ", "
      //                << ((static_cast<unsigned>(abs(lit)) << 1) ^ (lit < 0))
      //                << "  ";
      //      std::cout << lit << " ";

      if (abs(lit) >= nVars) {
        _solver->NewVariables(abs(lit) - nVars + 1);
        nVars = abs(lit) + 1;
        //        std::cout << "NoVars: " << nVars << std::endl;
      }
      //      unsigned var = (static_cast<unsigned>(abs(lit)) << 1) ^ (lit < 0);
      _solver->AddLiteral(&lit);
    }
    _solver->CommitClause();

    //    std::cout << std::endl;
  }

  if (_settings->verbosity > 0)
    std::cout << "c Adding now " << _softClauses.size() << " soft clauses !"
              << std::endl;

  // print all soft clauses
  for (auto softClause : _softClauses) {
    _solver->ResetClause();
    _solver->NewClause();
    if (softClause->relaxationLit == 0) {
      softClause->relaxationLit = _solver->NewVariable() << 1;
      nVars++;

      //      std::cout << "RL: " << softClause->relaxationLit << std::endl;
      //      softClause->clause.push_back(softClause->relaxationLit);
    }
    //    AddClause(softClause->clause);

    for (auto lit : softClause->clause) {
      //      std::cout << lit << ", " << (lit >> 1) << " ";
      //      std::cout << (lit) << " ";
      int literal = lit >> 1;
      if (!(lit ^ 1)) {
        literal = -literal;
      }
      //      std::cout << literal << ", ";
      _solver->AddLiteral(&literal);
      if ((lit >> 1) >= nVars) {
        nVars = (lit >> 1);
        _solver->NewVariables((lit >> 1) - nVars);
        //        std::cout << "NoVars: " << nVars << std::endl;
      }
    }
    _solver->AddLiteral(&softClause->relaxationLit);
    _solver->CommitClause();
    //    std::cout << std::endl;
  }

  //  std::cout << "Solve: " << _solver->Solve() << std::endl;
  //  std::cout << "Solve: " << Solve(assumptions) << std::endl;
  //  exit(1);
}

uint64_t GreedyPrepro::StartPrepro() {
  if (_settings->verbosity > 2) std::cout << __PRETTY_FUNCTION__ << std::endl;

  //      static_cast<unsigned>(_settings->greedyPPTimeLimit);
  unsigned maxRounds =
      static_cast<unsigned>(ceil(pow(_softClauses.size(), 0.3)) * 100);
  _timeLimit =
      _settings->greedyPPTimeoutFactor * pow(_softClauses.size(), 0.35);
  _preproPropagationLimit = static_cast<unsigned long>(
      _timeLimit * _settings->greedyPPPropagationPerSecond);

  //  _settings->greedyPPPropagationPerSecond
  //  _settings->greedyPPTimeoutFactor
  //  _settings->greedyPPSSCSwitchFactor

  //  _timeLimit = 12 * pow(_softClauses.size(), 0.35);
  //  _preproPropagationLimit = static_cast<unsigned long>(_timeLimit *
  //  2500000);

  GreedyMaxInitSATWeightV2(1, maxRounds);

  return _opti;
}

unsigned GreedyPrepro::GreedyMaxInitSATWeightV2(int greedyPrepro,
                                                unsigned maxRounds) {
  uint64_t softWeights = 0;
  std::for_each(_softClauses.begin(), _softClauses.end(),
                [&](SoftClause *SC) { softWeights += SC->weight; });
  _sumOfSoftWeights = softWeights;

  if (_settings->verbosity > 0)
    std::cout << "c _sumOfSoftWeights : " << _sumOfSoftWeights << std::endl;

  if (_settings->verbosity > 2) std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_settings->verbosity > 0) {
    std::cout << "c maxTime: " << _timeLimit << std::endl;
    std::cout << "c maxRounds: " << maxRounds << std::endl;
  }

  //  DGPW::TimeMeasurement timeSolvedFirst(&_timeVariables->solvedFirst, true);
  _timeSolvedFirst =
      new DGPW::TimeMeasurement(&_timeVariables->solvedFirst, true);
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

  //  uint64_t unsatisfiableSCWeight = 0;

  unsigned lastNewMaxRound = 0;
  bool isMaxSat = false;

  bool addedUnitClausesWithoutGettingNewMax = false;
  unsigned onlyOneSCSatByChance = 0;

  double lastNewMaxTime = 0;
  double lastSCDeletedTime = 0;
  _opti = static_cast<unsigned long>(-1);

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
    if (_settings->verbosity > 0)
      std::cout << "c it is unweighted MaxSAT!!!" << std::endl;
  }

  // first solve, with external assumptions, if given
  // could be the last satisfying model of another DGPW
  // for DivideDGPW mode!
  currentresult = Solve();
  if (greedyPrepro == 0) {
    return currentresult;
  }

  if (_settings->verbosity > 0) {
    std::cout << "c PrePro _timeLimit: " << _timeLimit << std::endl;
    std::cout << "c PrePro propagation limit: " << _preproPropagationLimit
              << std::endl;
  }

  if (_settings->verbosity > 1)
    std::cout << "c start solving model!" << std::endl;
  for (unsigned rounds = 0; rounds < maxRounds; ++rounds) {
    // informations...
    if (_settings->verbosity > 0) {
      std::cout << "c time of greedy prepro: "
                << _timeSolvedFirst->CurrentTimeDiff() << std::endl;
      if (_settings->verbosity > 4)
        std::cout << "c local Max: " << localMax
                  << "  rounds / 2: " << rounds / 2
                  << "  rounds - lastNewMaxRound: " << rounds - lastNewMaxRound
                  << "  ceil(maxRounds/40): " << ceil(maxRounds / 40)
                  << std::endl;
      if (_settings->verbosity > 1)
        std::cout << "c rounds: " << rounds << std::endl;
    }

    //        std::vector<unsigned int> indices;
    //        for (unsigned int i = 0; i < UNSATSCs.size(); i++) {
    //            indices.push_back(i);
    //        }

    assert(currentresult == SATISFIABLE);

    // get Info about the actual SC Model, the return values are sorted
    // highest weight first!
    tie(nextAssumptions, UNSATSCs, actualSATWeight) =
        SatisfiedSCsInfo(&sortedSCIndices);

    if (rounds == 0) {
      //      std::cout << "firstSATWeight: " << actualSATWeight << std::endl;
      //      std::cout << "_SATWeight: " << _satWeight << std::endl;
      firstSATWeight = actualSATWeight;
    }

    //        std::cout << "actualSATWeight: " << actualSATWeight << "
    //        _satWeight: " << _satWeight
    //                  << std::endl;

    if (actualSATWeight > _satWeight) {
      _satWeight = actualSATWeight;
      //      std::cout << "ASW: " << actualSATWeight << std::endl;

      addedUnitClausesWithoutGettingNewMax = false;
      _pacose->CalculateSATWeight();
      //      highestSATAssignment = GetLastSatisfiableAssignment();
      if (_settings->verbosity > 0) {
        std::cout << "c                                        NEW MAX FOUND: "
                  << _satWeight << std::endl;
        std::cout << "c local o " << _opti << std::endl;
      }
      lastNewMaxTime = _timeSolvedFirst->CurrentTimeDiff();
      lastNewMaxRound = rounds;
      if (_satWeight == _sumOfSoftWeights) {
        if (_settings->verbosity > 0)
          std::cout << "c All weights could be satisfied!!!" << std::endl;

        break;
      }
    }

    //    if ((_timeLimit < _timeSolvedFirst->CurrentTimeDiff() - lastNewMaxTime
    //    &&
    //         _timeLimit < _timeSolvedFirst->CurrentTimeDiff() -
    //         lastSCDeletedTime) ||
    if (_timeSolvedFirst->CurrentTimeDiff() > _timeLimit) {
      if (_settings->verbosity > 0)
        std::cout << "c timeout in greedy preprocessor!" << std::endl;
      break;
    } else if (_fixSoftClauses != 0 &&
               (_timeSolvedFirst->CurrentTimeDiff() >
                (_timeLimit * _settings->greedyPPSSCSwitchFactor * 0.01))) {
      // after 2/3 of the time, change strategy to only satisfy all SC at least
      // once!
      if (_settings->verbosity > 0) {
        std::cout << "_timeSolvedFirst->CurrentTimeDiff() > (_timeLimit * "
                     "_settings->greedyPPSSCSwitchFactor * 0.01))"
                  << _timeSolvedFirst->CurrentTimeDiff() << " < " << _timeLimit
                  << " * " << _settings->greedyPPSSCSwitchFactor << " * 0.01 "
                  << std::endl;
      }
      _fixSoftClauses = 0;
    }

    //    currentresult = TestNextSoftClauses(nextAssumptions, UNSATSCs,
    //    isMaxSat,
    //                                        tryToSolveMoreThanOneWeight);

    //    std::cout << std::setw(100) << "previousMax: " << previousMax <<
    //    std::endl; std::cout << std::setw(100) << "localMax: " << localMax <<
    //    std::endl;
    if (previousMax != localMax) {
      _previousPosition = 1 - (_settings->greedyPPSATPercentage * 0.01);
      _previousLowerBound = 0.0;
    }
    if (_fixSoftClauses == 0 || (_fixSoftClauses == 1 && localMax > 0)) {
      nextAssumptions.clear();
      neverSATSCs = BuildOrderedIntersection(&UNSATSCs, &neverSATSCs);
      currentresult = BinarySearchSatisfySCs(nextAssumptions, &neverSATSCs);
    } else {
      currentresult = BinarySearchSatisfySCs(nextAssumptions, &UNSATSCs);
    }
    previousMax = localMax;

    if (currentresult != SATISFIABLE) {
      localMax++;
      if (_settings->verbosity > 1 && currentresult == UNSAT)
        std::cout << "c                                              The "
                  << localMax << " th local Max is found: " << actualSATWeight
                  << std::endl;

      neverSATSCs = BuildOrderedIntersection(&UNSATSCs, &neverSATSCs);
      if (neverSATSCs.empty()) {
        if (_settings->verbosity > 1)
          std::cout << "c All clauses were at least once satisfied!"
                    << std::endl;
        break;
      }

      if (_timeSolvedFirst->CurrentTimeDiff() > _timeLimit) {
        if (_settings->verbosity > 0)
          std::cout << "c timeout in greedy preprocessor!" << std::endl;
        break;
      }
      if (_settings->verbosity > 1)
        std::cout << "c neverSATSCs.size(): " << neverSATSCs.size()
                  << std::endl;

      //      if (_noClauses != 1 || localMax == 1) {
      if ((currentresult == UNKNOWN && _fixSoftClauses == 0) ||
          (_fixSoftClauses == 0 && _noClauses != 1) ||
          (_fixSoftClauses == 1 && (localMax == 1 || _noClauses != 1))) {
        nextAssumptions.clear();
        currentresult = Solve(nextAssumptions);
        continue;
      }

      for (unsigned iter = 0; iter < neverSATSCs.size(); ++iter) {
        if (_fixSoftClauses == 2) {
          // no additional SCs can be satisfied!
          // all remaining SCs are unsat!

          nextAssumptions.clear();
          nextAssumptions.push_back(neverSATSCs[iter]->relaxationLit ^ 1);
          //                  _solver->DeactivateLimits();
          currentresult = Solve(nextAssumptions);
          if (currentresult == SATISFIABLE) {
            if (_settings->verbosity > 0) {
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
          _unsatisfiableSCWeight += neverSATSCs[iter]->weight;

          if (_settings->verbosity > 0) {
            std::cout << "c Softclause with weight "
                      << neverSATSCs[iter]->weight
                      << " couldn't be satisfied! - Remove this softclause! - "
                         "new SC.size: "
                      << _softClauses.size() - 1 << "  neverSATSCs.size() "
                      << neverSATSCs.size() << "  UNSATSCs.size() "
                      << UNSATSCs.size()
                      << "   sortedSize: " << sortedSCIndices.size()
                      << std::endl;
          }
          addedUnitClausesWithoutGettingNewMax = true;
          lastSCDeletedTime = _timeSolvedFirst->CurrentTimeDiff();
          std::vector<uint32_t> unitclause;
          unitclause.push_back(neverSATSCs[iter]->relaxationLit);
          AddClause(unitclause);

          //           find SC in sortedSCIndices
          //                    for (unsigned j = 0; j < sortedSCIndices.size();
          //                    j++) {
          //                      if (neverSATSCs[iter]->relaxationLit ==
          //                          _softClauses[sortedSCIndices[j]]->relaxationLit)
          //                          {
          //                        std::cout << sortedSCIndices[j] <<
          //                        std::endl; std::cout <<
          //                        *(sortedSCIndices.begin() + j) << std::endl;
          //                        // remove iter from list and break
          //                        _sumOfSoftWeights -=
          //                        _softClauses[sortedSCIndices[j]]->weight;
          //                        if (sortedSCIndices[j] !=
          //                        _softClauses.size())
          //                          _softClauses.erase(_softClauses.begin() +
          //                          sortedSCIndices[j]);
          //                        else
          //                          _softClauses.pop_back();
          //                        sortedSCIndices.erase(sortedSCIndices.begin()
          //                        + j);
          //                        // the indices may be too big - shrink them
          //                        if bigger than j
          //                            for (unsigned i = 0; i <
          //                        sortedSCIndices.size(); i++) {
          //                          if (sortedSCIndices[i] > j)
          //                          sortedSCIndices[i]--;
          //                        }
          //                        break;
          //                      }
          //                    }

          for (unsigned j = 0; j < _softClauses.size(); ++j) {
            if (_softClauses[j]->relaxationLit ==
                neverSATSCs[iter]->relaxationLit) {
              _sumOfSoftWeights -= _softClauses[j]->weight;
              _softClauses.erase(_softClauses.begin() + j);
              break;
            }
          }

          neverSATSCs.erase(neverSATSCs.begin() + iter);
          _pos1Unsat = false;

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
        if (_timeSolvedFirst->CurrentTimeDiff() > _timeLimit) {
          if (_settings->verbosity > 0)
            std::cout << "c timeout in greedy preprocessor!" << std::endl;
          break;
        }
      }

      if (currentresult != SATISFIABLE) {
        if (_settings->verbosity > 0)
          std::cout << "c no more softclauses can be satisfieda!" << std::endl;
        break;
      }
    }

    //    if (neverSATSCs.empty() && localMax > 0) {

    if (neverSATSCs.empty() &&
        (_fixSoftClauses == 0 || (_fixSoftClauses != 1 && localMax > 0))) {
      if (_settings->verbosity > 0)
        std::cout << "c no more softclauses are satisfiable!" << std::endl;

      break;
    }
  }

  //      _solver->DeactivateLimits();
  // TODO!! maybe faster without this line!
  //  if (addedUnitClausesWithoutGettingNewMax) {
  //    currentresult = Solve();
  //    if (currentresult != 10) {
  //      std::cout << "c Strange result UNSAT after greedy Prepro, shouldn't
  //      be!"
  //                << std::endl;
  //      assert(false);
  //    }
  //    _pacose->CalculateSATWeight();
  //    tie(nextAssumptions, UNSATSCs, actualSATWeight) =
  //        SatisfiedSCsInfo(&sortedSCIndices);
  //  }

  while (!_softClauses.empty() &&
         _opti < _softClauses[sortedSCIndices.back()]->weight) {
    //    std::cout << "c XXA SC can be removed because it is SAT anyway!!!!"
    //              << std::endl;
    //    std::cout << "c o=" << _sumOfSoftWeights - actualSATWeight
    //              << "  _softClauses[sortedSCIndices.back()]->weight: "
    //              << _softClauses[sortedSCIndices.back()]->weight
    //              << "  _sumOfSoftWeights: " << _sumOfSoftWeights
    //              << "  actualSATWeight: " << actualSATWeight << std::endl;

    _satisfiableSCs++;

    std::vector<uint32_t> unitclause;
    unitclause.push_back(_softClauses[sortedSCIndices.back()]->relaxationLit ^
                         1);
    AddClause(unitclause);

    for (unsigned i = 0; i < sortedSCIndices.size(); i++) {
      if (sortedSCIndices[i] > sortedSCIndices.back()) sortedSCIndices[i]--;
    }
    //    _sumOfSoftWeights -= _softClauses[sortedSCIndices.back()]->weight;
    //    _satWeight -= _softClauses[sortedSCIndices.back()]->weight;

    if (_settings->verbosity > 0) {
      std::cout << "c Softclause with weight "
                << _softClauses[sortedSCIndices.back()]->weight
                << " is always satisfied! - Remove this softclause!"
                << std::endl;
    }
    _softClauses.erase(_softClauses.begin() + sortedSCIndices.back());
    sortedSCIndices.pop_back();

    if (_softClauses.empty() ||
        _opti >= _softClauses[sortedSCIndices.back()]->weight) {
      _solver->ClearAssumption();
      unsigned rv = _solver->Solve();
      if (rv != SATISFIABLE) {
        std::cout << "ERROR: SHOULD BE SATISFIABLE IN GREEDY PPA!" << std::endl;
        exit(1);
      }
      tie(nextAssumptions, UNSATSCs, actualSATWeight) =
          SatisfiedSCsInfo(&sortedSCIndices);
    }
  }

  if (_settings->verbosity > 0) {
    if (_satisfiableSCs > 0) {
      std::cout << "c in any case SAT SCs....: " << _satisfiableSCs
                << std::endl;
    }

    if (_unsatisfiableSCs > 0) {
      std::cout << "c unsatisfiable SCs......: " << _unsatisfiableSCs
                << std::endl;
      std::cout << "c remaining SCs..........: " << _softClauses.size()
                << std::endl;
      std::cout << "c weight of unsat SCs....: " << _unsatisfiableSCWeight
                << std::endl;
    }

    if (onlyOneSCSatByChance > 0)
      std::cout << "c only 1 SC SAT by chance: " << onlyOneSCSatByChance
                << std::endl;

    //  std::cout << "c max SatWeight Prepro...: " << actualSATWeight <<
    //  std::endl;

    std::cout << "c max SatWeight Prepro...: " << _satWeight << std::endl;

    std::cout << "c sumOfSoftWeights.......: " << _sumOfSoftWeights
              << std::endl;

    std::cout << "c localPrePro o..........: " << _opti << std::endl;
    //  std::cout << "c add to 1st solved......: " << _satWeight -
    //  firstSATWeight
    //            << std::endl;
    std::cout << "c time greedyPrepro......: "
              << _timeSolvedFirst->CurrentTimeDiff() << std::endl;
  }

  return SATISFIABLE;
}

std::tuple<std::vector<unsigned>, std::vector<SoftClause *>, uint64_t>
GreedyPrepro::SatisfiedSCsInfo(std::vector<unsigned> *sortedSCIndices) {
  if (_settings->verbosity > 2) std::cout << __PRETTY_FUNCTION__ << std::endl;

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
  if (_settings->verbosity > 0)
    std::cout << "c local o: " << unsatWeight << std::endl;
  //    _sumOfSoftWeights << std::endl;

  _sumOfSoftWeights = unsatWeight + actualWeight;

  if (_opti > unsatWeight || _opti == static_cast<unsigned long>(-1))
    _opti = unsatWeight;

  //  if (actualWeight > _satWeight) {
  //    _satWeight = actualWeight;
  //  }
  //    std::cout << "DONE!" << std::endl;
  return std::make_tuple(nextAssumptions, UNSATSCs, actualWeight);
}

std::vector<SoftClause *> GreedyPrepro::BuildOrderedIntersection(
    std::vector<SoftClause *> *UNSATSCs,
    std::vector<SoftClause *> *neverSATSCs) {
  if (_settings->verbosity > 2) std::cout << __PRETTY_FUNCTION__ << std::endl;

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

unsigned GreedyPrepro::BinarySearchSatisfySCs(
    std::vector<unsigned> &nextAssumptions,
    std::vector<SoftClause *> *unsatSCs) {
  if (_settings->verbosity > 2) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    std::cout << "c greedy strategy: kind of binary search to force optimal "
                 "number of SC to be true!"
              << std::endl;
  }
  double noVars(-1);

  if (unsatSCs->empty()) {
    if (_settings->verbosity > 1) {
      std::cout << "No more SCs to satisfy!" << std::endl;
    }
    return SATISFIABLE;
  }

  //  unsatSCs = &_softClauses;

  if (_settings->verbosity > 0) {
    std::cout << std::setw(30) << "_previousPosition: " << std::setw(10)
              << _previousPosition << std::endl;
    std::cout << std::setw(30) << "_previousLowerBound: " << std::setw(10)
              << _previousLowerBound << std::endl;
  }
  assert(_previousPosition <= 1);
  assert(_previousPosition > 0);

  unsigned noUnsatSCs = unsatSCs->size();
  //  noUnsatSCs = 25;
  //  _previousPosition = 0.18;

  double sqrt = std::sqrt(noUnsatSCs);
  unsigned numberOfPositions = (floor(sqrt * 2)) - 1;
  if (noUnsatSCs == 2) {
    numberOfPositions++;
  }

  unsigned newPosition = 0;

  if (_settings->verbosity > 5) {
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
      noVars = noUnsatSCs / _noClauses;
      std::cout << "Position, noClauses, noVars: " << std::setw(5) << i
                << std::setw(8) << _noClauses << std::setw(8)
                << std::round(noVars * 100) / 100 << std::endl;
    }
  }

  newPosition = ceil(_previousPosition * numberOfPositions);
  if (_pos1Unsat && newPosition == 1 && numberOfPositions > 1) newPosition++;

  if (newPosition < sqrt) {
    //    noClauses = ceil(noUnsatSCs / newPosition);
    assert(newPosition != 0);
    //    if (newPosition == 0) newPosition = 1;
    _noClauses = noUnsatSCs / newPosition;
  } else {
    _noClauses = numberOfPositions - newPosition + 1;
  }
  noVars = (noUnsatSCs / _noClauses);
  if (_settings->verbosity > 2) {
    std::cout << std::endl
              << "Position, noClauses, noVars: " << std::setw(5) << newPosition
              << std::setw(8) << _noClauses << std::setw(8)
              << std::round(noVars * 100) / 100 << std::endl;
  }

  unsigned upperBound = noUnsatSCs % static_cast<int>(floor(noVars));
  if (_settings->verbosity > 3) {
    std::cout << "Number of clauses  X  Clause Size: " << std::setw(10)
              << upperBound << " X " << ceil(noVars) << std::endl;
    std::cout << "Number of clauses  X  Clause Size: " << std::setw(10)
              << _noClauses - upperBound << "  X  " << floor(noVars)
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

  // unsigned seed =
  // std::chrono::system_clock::now().time_since_epoch().count();
  unsigned seed = 1234567890;
  std::shuffle(indices.begin(), indices.end(),
               std::default_random_engine(seed));

  int varsPerClause = ceil(noVars);
  //  std::cout << "unsatSCs->size(): " << unsatSCs->size()
  //            << " clauses.size(): " << clauses.size() << std::endl;

  int i = 0;
  for (unsigned j = 0; j < clauses.size(); j++) {
    if (j == upperBound) varsPerClause = floor(noVars);
    for (int k = 0; k < varsPerClause; k++) {
      //      std::cout << "i: " << i << " varsPerClause: " << varsPerClause
      //                << " j: " << j << " UB: " << upperBound << "; ";
      if (_settings->verbosity > 4)
        std::cout << ((*unsatSCs)[indices[i]]->relaxationLit ^ 1) << " ";
      //                << std::endl;
      clauses[j].push_back((*unsatSCs)[indices[i]]->relaxationLit ^ 1);
      i++;
    }
    AddClause(clauses[j]);
    if (_settings->verbosity > 4) std::cout << std::endl;
  }

  assert(i == unsatSCs->size());
  //  unsigned nextAssumptions nextAssumptions.push_back(clauses[index].back() ^
  //  1);
  //    AddClause(clauses[index]);
  //  }

  //  //      _solver->SetMaxPropagations(_preproPropagationLimit);

  nextAssumptions.push_back(relaxlit ^ 1);

  //  unsigned currentresult = Solve(nextAssumptions);
  double solvingTime = _timeSolvedFirst->CurrentTimeDiff();
  _solver->SetPropagationBudget(_preproPropagationLimit);
  unsigned currentresult = SolveLimited(nextAssumptions);
  solvingTime = _timeSolvedFirst->CurrentTimeDiff() - solvingTime;
  if (_settings->verbosity > 0)
    std::cout << "c Solver Call Time: " << solvingTime << std::endl;
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

  bool firstTimeLimitReached =
      (_timeSolvedFirst->CurrentTimeDiff() >
       (_timeLimit * _settings->greedyPPSSCSwitchFactor * 0.01));
  //  bool firstTimeLimitReached =
  //      (_timeSolvedFirst->CurrentTimeDiff() > (_timeLimit / 3) * 2);
  bool timeLimitReached = (_timeSolvedFirst->CurrentTimeDiff() > _timeLimit);

  if (currentresult == SATISFIABLE) {
    if (_settings->verbosity > 0) {
      std::cout << std::setw(100) << "SAT" << std::endl;
      if (noVars == 1) {
        std::cout << std::setw(100) << "c ALL REMAINING SCs COULD BE SATISFIED!"
                  << std::endl;
      }
    }
    _previousPosition -= (_previousPosition - _previousLowerBound) /
                         (1 / (_settings->greedyPPSATPercentage * 0.01));
    if (_previousPosition == 0) _previousPosition = 0.000001;

    return currentresult;
  } else if (currentresult == UNSAT) {
    if (_settings->verbosity > 0)
      std::cout << std::setw(100) << "UNSAT" << std::endl;
    _previousLowerBound = _previousPosition;
    _previousPosition += (1 - _previousLowerBound) /
                         (1 / (_settings->greedyPPUNSATPercentage * 0.01));
    if (_previousPosition > 1) _previousPosition = 1;

    // no more soft clauses can be satisfied
    if (_noClauses == 1) {
      if (_settings->verbosity > 0)
        std::cout
            << std::setw(100)
            << "c NO ADDITIONAL SOFTCLAUSES FROM THIS SET ARE SATISFIABLE!"
            << std::endl;

      return currentresult;
    }

    if (newPosition == 1) {
      // we do not try again to satisfy all remaining softclauses at once
      _pos1Unsat = true;
    }

    if (timeLimitReached || (firstTimeLimitReached && _fixSoftClauses == 2)) {
      return UNKNOWN;
    }
    return BinarySearchSatisfySCs(nextAssumptions, unsatSCs);
  } else {
    _unknownSolverCalls++;
    if (_settings->verbosity > 0)
      std::cout << std::setw(100)
                << "PROPAGATION LIMIT REACHED FOR SINGLE SOLVER CALL"
                << std::endl;
    _previousLowerBound = _previousPosition;
    _previousPosition += (1 - _previousLowerBound) /
                         (1 / (_settings->greedyPPUNSATPercentage * 0.01));
    _previousLowerBound = _previousPosition;
    _previousPosition += (1 - _previousLowerBound) /
                         (1 / (_settings->greedyPPUNSATPercentage * 0.01));
    if (_previousPosition > 1) _previousPosition = 1;

    // no more soft clauses can be satisfied
    if (_noClauses == 1 && _fixSoftClauses == 0) {
      // the last remaining SC could not be satisfied, reduce timelimit to kill
      // prepro!
      _timeLimit = 0;
      return currentresult;
    } else if (_noClauses == 1) {
      _fixSoftClauses = 0;
      return currentresult;
    }

    if (_unknownSolverCalls < 2 &&
        ((!firstTimeLimitReached && _fixSoftClauses == 2) ||
         !timeLimitReached)) {
      return BinarySearchSatisfySCs(nextAssumptions, unsatSCs);
    } else if (_unknownSolverCalls == 2 || firstTimeLimitReached) {
      if (_fixSoftClauses != 0 || firstTimeLimitReached) {
        _fixSoftClauses = 0;
        return currentresult;
      } else {
        return BinarySearchSatisfySCs(nextAssumptions, unsatSCs);
      }
    } else if (_unknownSolverCalls > 3) {
      _timeLimit = 0;
      return currentresult;
    }
  }
  return 0;
}
