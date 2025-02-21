/********************************************************************************************
Copyright (c) 2014-2016, Sven Reimer,
Copyright (c) 2017-2020, Tobias Paxian

dPermission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights to
    use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is furnished
to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
********************************************************************************************/
#include <unistd.h>
#include <algorithm>  //std::sort
#include <cmath>
#include <iomanip>
#include <numeric>

// Include dgpw related headers.
#include "bucket.h"
#include "cascade.h"
#include "dgpw.h"
#include "multiplecascade.h"
#include "softclausenodes.h"
#include "timemeasurement.h"
#include "timevariables.h"
#include "totalizerencodetree.h"

namespace DGPW {

// Constructor
Cascade::Cascade(DGPW *dgpw, MultipleCascade *multipleCascade, bool onlyByTares)
    : _dgpw(dgpw),
      //_control(dgpw->_control),
      _setting(dgpw->_dgpwSetting),
      _multipleCascade(multipleCascade),
      _base(_setting->base),
      _onlyByTares(onlyByTares),
      _maxPos(-1),
      _satWeight(0),
      _tareWeight(0),
      _weightToSubstract(0),
      _sumOfSoftWeights(0),
      _softClauseTreeCreated(false),
      _collectedCascadeAssumptions(),
      _structure(),
      _numberOfBuckets(0),
      _totalBucketEntriesperWeight(0),
      _totalBucketOccurrences(0),
      _totalBucketEntries(),
      _maxSorterDepth(0),
      _highestBucketMultiplicator(0),
      _upperWeightBoundAllLowerCascades(0),
      _howManyPositionsCuttedAtBottom(0),
      _softClauses(),
      _softClauseTree(),
      _processingSoftClauseTree(),
      _processingPercentOffTree(),
      _howOftenReinsertedFromProcessingPercentOffTree(0),
      _tareAssumptions(),
      _fixedTareAssumption() {
  assert(dgpw != nullptr);
}

void Cascade::Fill(std::vector<SoftClause *> *softClauses,
                   PartitionStrategy partitionStrategy,
                   EncodeStrategy encodeStrategy) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_onlyByTares) encodeStrategy = ENCODEONLYIFNEEDED;

  TimeMeasurement timeFillingBuckets(&_dgpw->_timeVariables->fillingBuckets,
                                     true);

  CountSumOfSoftWeights(softClauses);

  FillStructure(partitionStrategy, encodeStrategy);

  if (_setting->verbosity > 3)
    std::cout << std::endl
              << "Buckets are filled - structure is dumped" << std::endl
              << std::endl;
}

void Cascade::CountSumOfSoftWeights(std::vector<SoftClause *> *softClauses) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  _softClauses = *softClauses;
  _sumOfSoftWeights = 0;

  // sum over all SC weights of CASCADE!
  std::for_each(_softClauses.begin(), _softClauses.end(),
                [&](SoftClause *SC) { _sumOfSoftWeights += SC->weight; });
  _satWeight = CountSatisfiedSoftClauses(nullptr, _dgpw->_lastModel);

  if (_setting->verbosity > 4) {
    std::cout << "c SAT Weight of Cascade..: " << _satWeight << std::endl;
    std::cout << "c SAT Weight of dgpw.....: " << _dgpw->_satWeight
              << std::endl;
    std::cout << "c softWeights of Cascade.: " << _sumOfSoftWeights
              << std::endl;
    std::cout << "c softWeights of dgpw....: " << _dgpw->_sumOfSoftWeights
              << std::endl;
  }
}

uint32_t Cascade::CutMaxPos(bool solve) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  _structure.back()->_isLastBucket = true;
  return _structure.back()->CutMaxPos(solve);
}

uint32_t Cascade::CutMinPos(bool solve) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  _structure.back()->_isLastBucket = true;
  return _structure.back()->CutMinPos(solve);
}

void Cascade::FillStructure(PartitionStrategy partitionStrategy,
                            EncodeStrategy encodeStrategy) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_onlyByTares) encodeStrategy = ENCODEONLYIFNEEDED;

  PartitionSoftClauses(partitionStrategy);
  FillBuckets();

  AddTaresToBuckets();

  if (encodeStrategy == ENCODEONLYIFNEEDED) {
    UnionBucketsIntoLast();
    if (_onlyByTares) AddAsManyBucketsAsPossible();
    DumpBucketStructure(true, 3);

  } else {
    DumpBucketStructure(false, 3);
  }
}

void Cascade::AddAsManyBucketsAsPossible() {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  //    std::cout << "size: " << _structure.back()->size() << std::endl;
  //    std::cout << "_structure.back()->_tares[0]: " <<
  //    _structure.back()->_tares[0] << std::endl;

  if (_structure.back()->_tares.empty()) AddTare(_structure.size() - 1);

  if (_setting->interimResult == CUTATTOP) {
    CutMaxPos();
  }
  _structure.back()->_encodeTreeGenerated = false;

  /** ATTENTION NOT ALWAYS GIVEN
   *
   * If e.g. highest number is close to uint64_t...
   */
  if (AddNewBucketsTillMultiplicatorMatches(static_cast<uint64_t>(-1), true)) {
    std::cout << "EXCEPTION - cannot be solved only by tares because of "
                 "multiplicator limit of "
                 "64 Bit!";
    _onlyByTares = false;
    return;
  }
  DumpBucketStructure(true);

  _estimatedWeightBoundaries[0] = -_highestBucketMultiplicator + 1;
  _estimatedWeightBoundaries[1] = _highestBucketMultiplicator;
}

bool Cascade::Encode() {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  switch (_setting->encodeStrategy) {
    case ENCODEALL:
      EncodeTopBuckets();
      /*
      if ( _control->ReachedLimits() )
          return false;
              */

      EncodeBottomBuckets();
      CalculateBucketEntries();
      DumpBucketStructure(false, 4);
      break;
    case ENCODEONLYIFNEEDED:
      if (_onlyByTares) return true;
      CreateTotalizerEncodeTree();
      CalculateBucketEntries();
      DumpBucketStructure(true, 4);
      break;
  }
  /*
  if ( _control->ReachedLimits() )
      return false;
      */

  if (_setting->verbosity < 4) return true;

  std::cout << std::endl
            << "Buckets are encoded - structure is dumped" << std::endl
            << std::endl;
  return true;
}

uint32_t Cascade::Solve(bool onlyWithAssumptions, bool solveTares) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  uint32_t currentresult(1);

  if (_dgpw->_satWeight == _dgpw->_sumOfSoftWeights) {
    _maxPos = _structure.back()->size() - 1;
    return SATISFIABLE;
  }

  //    TimeMeasurement
  //    TimeSolvingLastBucket(&_dgpw->_timeVariables->solvingLastBucket);
  if (!_onlyByTares) {
    _maxPos = static_cast<int>(_structure.back()->SolveBucketReturnMaxPosition(
        onlyWithAssumptions, false));
    if (_maxPos != -1 && onlyWithAssumptions)
      _fixedTareAssumption.push_back(
          (_structure.back()->_sorter->GetOrEncodeOutput(_maxPos) << 1) ^ 1);
    //        std::cout << "literal of mP == " << litmP << std::endl <<
    //        std::endl;
    //        _dgpw->AddUnit((_structure.back()->_sorter->GetOrEncodeOutput(_maxPos)
    //        << 1) ^ 1); exit(0);

    if (_dgpw->_resultUnknown) {
      return UNKNOWN;
    }
    if (_setting->encodeStrategy == ENCODEONLYIFNEEDED &&
        _setting->createGraphFile != "")
      _structure.back()->_sorter->_outputTree->DumpOutputTree(
          _setting->createGraphFile + "_withOutputs.tgf", true);

    //    if (_setting->encodeStrategy != ENCODEALL) {
    //      std::cout << "c #clauses of coding.....: " << _dgpw->_addedClauses
    //                << std::endl;
    //      std::cout << "c #variables of coding...: "
    //                << _dgpw->Variables() - _dgpw->_variablesBefore <<
    //                std::endl;
    //    }
  }
  if (_dgpw->_satWeight == _dgpw->_sumOfSoftWeights) {
    //        std::cout << "_dgpw->_satWeight == _dgpw->_sumOfSoftWeights" <<
    //        std::endl;
    return SATISFIABLE;
  }
  //        std::cout << "_dgpw->_satWeight != _dgpw->_sumOfSoftWeights" <<
  //        std::endl;

  //    _dgpw->_glucose->printIncrementalStats(1);
  //      TimeSolvingLastBucket.~TimeMeasurement();

  if (solveTares) {
    //        std::cout << "c CURRENTRESULT: SOLVE TARES! " << currentresult <<
    //        std::endl;
    currentresult = SolveTares(onlyWithAssumptions);
    //        std::cout << "c CURRENTRESULT: AFTER SOLVE TARES! " <<
    //        currentresult << std::endl;
  }
  //    currentresult = SolveAllTares();

  return currentresult;
}

void Cascade::CreateSoftClauseTree(std::vector<SoftClause *> *softClauses,
                                   bool split) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  _processingSoftClauseTree.clear();
  _softClauseTree.clear();

  // sort without changing order of elements
  //    std::stable_sort(softClauses->begin(), softClauses->end(),
  //    SoftClause::bigger ); _numberOfBuckets =
  //    static_cast<uint32_t>(floor(log2(_softClauses[0]->weight)/log2(_base)));

  //    if (_setting->featureTest)
  //    {
  //        //sort with changing order of elements (better worst case runtime)
  //        std::sort(softClauses->begin(), softClauses->end(),
  //        SoftClause::bigger ); _numberOfBuckets =
  //        static_cast<uint32_t>(floor(log2(_softClauses[0]->weight)/log2(_base)));
  //    } else
  //    {
  uint64_t maxValue = 0;
  for (auto softclause : *softClauses) {
    maxValue = softclause->weight > maxValue ? softclause->weight : maxValue;
  }
  _numberOfBuckets = static_cast<uint32_t>(floor(log2(maxValue) / log2(_base)));

  // SoftClauseNode Structure is created. Extract function!
  for (uint32_t i = 0; i != softClauses->size(); ++i) {
    SoftClauseNodes *softClauseNode =
        new SoftClauseNodes((*softClauses)[i], _base);
    // softclause is only in one bucket
    if (split && softClauseNode->inHowManyBuckets != 1) {
      _processingSoftClauseTree.push_back(softClauseNode);
    } else {
      _softClauseTree.push_back(softClauseNode);
    }
  }
  if (split) _softClauseTreeCreated = true;
}

void Cascade::PartitionSoftClauseTree(
    std::vector<SoftClauseNodes *> *tmpSoftClausesTree) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  _processingSoftClauseTree.clear();
  _softClauseTree.clear();

  // SoftClauseNode Structure is created. Extract function!
  for (uint32_t i = 0; i != tmpSoftClausesTree->size(); ++i) {
    // SoftClauseNodes* softClauseNode = new SoftClauseNodes(_softClauses[i],
    // _base);
    // softclause is only in one bucket

    if ((*tmpSoftClausesTree)[i]->inHowManyBuckets != 1) {
      _processingSoftClauseTree.push_back((*tmpSoftClausesTree)[i]);
    } else {
      _softClauseTree.push_back((*tmpSoftClausesTree)[i]);
    }
  }
  _softClauseTreeCreated = true;
}

std::vector<std::vector<uint32_t>> Cascade::GetTareVectors() {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  std::vector<std::vector<uint32_t>> tares;

  uint32_t addNoTareToLastBucket = (_setting->cascadeDivider > 0) ? 0 : 1;
  if (_setting->verbosity > 3)
    std::cout << "addNoTareToLastBucket: " << addNoTareToLastBucket
              << std::endl;
  for (uint32_t ind = 0; ind < _structure.size() - addNoTareToLastBucket;
       ind++) {
    tares.push_back(_structure[ind]->_tares);
  }
  if (_setting->verbosity > 3) std::cout << "returning tares" << std::endl;
  return tares;
}

std::vector<std::pair<uint64_t, uint32_t>> Cascade::GetTareVector(
    uint64_t weightDiff) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  int64_t wd = static_cast<int64_t>(weightDiff);
  std::vector<std::pair<uint64_t, uint32_t>> tares = {};
  if (_softClauses.size() == 0) {
    return tares;
  }

  uint32_t addNoTareToLastBucket = (_setting->cascadeDivider > 0) ? 0 : 1;
  if (_setting->verbosity > 3)
    std::cout << "addNoTareToLastBucket: " << addNoTareToLastBucket
              << std::endl;
  for (uint32_t ind = 0; ind < _structure.size() - addNoTareToLastBucket;
       ind++) {
    if (wd > 0) {
      tares.push_back(std::pair<uint64_t, uint32_t>(
          _structure[ind]->_multiplicator * _dgpw->_greatestCommonDivisor,
          _structure[ind]->_tares[0] << 1));
      // if I only solve the watchdogs I have to return all tares!
      // else test if tare is used to minimize weight, then we can lift this
      // weight!
      if (!(_dgpw->_dgpwSetting->divideDGPW == DIVIDEALLSOLVEONLYWATCHDOGS) &&
          _dgpw->Model(_structure[ind]->_tares[0]) ==
              ((_structure[ind]->_tares[0] << 1) ^ 0)) {
        wd -= static_cast<int64_t>(_structure[ind]->_multiplicator) *
              _dgpw->_greatestCommonDivisor;
      }
      if (_setting->verbosity > 4)
        std::cout << "WD: " << wd
                  << "  multi: " << _structure[ind]->_multiplicator
                  << "  tare: " << _structure[ind]->_tares[0] * 2 << std::endl;
    } else {
      // tare doesn't need to be added anymore, can be set to the model value!
      _dgpw->AddUnit(_dgpw->Model(_structure[ind]->_tares[0]));
    }
  }

  if (_setting->verbosity > 3) std::cout << "c returning tares" << std::endl;
  return tares;
}

std::vector<std::pair<uint64_t, uint32_t>> Cascade::GetWatchdogs(
    uint64_t weightDiff) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  std::vector<std::pair<uint64_t, uint32_t>> watchdogs = {};
  if (_softClauses.size() == 0) {
    return watchdogs;
  }

  if (_setting->verbosity > 1) {
    std::cout
        << "(_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) - weightDiff: "
        << (_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) - weightDiff
        << std::endl;
    std::cout << "_dgpw->_greatestCommonDivisor * (_dgpw->_satWeight - "
                 "_estimatedWeightBoundaries[0]): "
              << _dgpw->_greatestCommonDivisor *
                     (_dgpw->_satWeight - _estimatedWeightBoundaries[0])
              << std::endl;

    std::cout << "((_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) - "
                 "weightDiff) / "
                 "(_structure.back()->_multiplicator * "
                 "_dgpw->_greatestCommonDivisor): "
              << ((_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) -
                  weightDiff) /
                     (_structure.back()->_multiplicator *
                      _dgpw->_greatestCommonDivisor)
              << std::endl;
    std::cout << "floor(((_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) - "
                 "weightDiff) / "
                 "(_structure.back()->_multiplicator * "
                 "_dgpw->_greatestCommonDivisor)): "
              << floor(((_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) -
                        weightDiff) /
                       (_structure.back()->_multiplicator *
                        _dgpw->_greatestCommonDivisor))
              << std::endl;
  }

  unsigned firstWatchdog = 0;
  unsigned lastWatchdog =
      static_cast<unsigned>(_maxPos) < _structure.back()->size()
          ? static_cast<unsigned>(_maxPos)
          : _structure.back()->size() - 1;

  if (static_cast<unsigned long>(_dgpw->_greatestCommonDivisor) *
          _dgpw->_satWeight >
      weightDiff) {
    firstWatchdog = static_cast<unsigned>(floor(
        ((_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) - weightDiff) /
        (_structure.back()->_multiplicator * _dgpw->_greatestCommonDivisor)));
  } /*else {
      firstWatchdog = 0;
  }*/

  if (_setting->verbosity > 0) {
    std::cout << "firstWatchdog: " << firstWatchdog << std::endl;
    std::cout << "lastWatchdog: " << lastWatchdog << std::endl;
    std::cout << "(_dgpw->_greatestCommonDivisor * _dgpw->_satWeight): "
              << static_cast<long int>(_dgpw->_greatestCommonDivisor) *
                     static_cast<long int>(_dgpw->_satWeight)
              << std::endl;
    std::cout << std::endl;
    std::cout << "weightDiff: " << weightDiff << std::endl;
    std::cout << "_structure.back()->_multiplicator: "
              << _structure.back()->_multiplicator << std::endl;
    std::cout << "_structure.back()->_localMaxPos: "
              << _structure.back()->_localMaxPos << std::endl;
    std::cout << "weight boundaries: ( " << _estimatedWeightBoundaries[0]
              << " / " << _estimatedWeightBoundaries[1] << " )" << std::endl;
    std::cout << "_maxPos: " << _maxPos << std::endl;
    std::cout << "dgpwSatWeight: " << _dgpw->_satWeight << std::endl;
    std::cout << "_GGT: " << _dgpw->_greatestCommonDivisor << std::endl
              << std::endl;
  }
  assert(firstWatchdog <= lastWatchdog);

  if (_maxPos >= 0 || (firstWatchdog) > 0) {
    std::cout << "c watchdog of position " << firstWatchdog
              << " added as unit clause! " << std::endl;
    _dgpw->AddUnit(
        (_structure.back()->_sorter->GetOrEncodeOutput(firstWatchdog) << 1) ^
        1);
  }

  firstWatchdog++;
  if (firstWatchdog <= lastWatchdog) {
    for (unsigned i = firstWatchdog; i <= lastWatchdog; ++i) {
      watchdogs.push_back(std::pair<uint64_t, uint32_t>(
          _structure.back()->_multiplicator * _dgpw->_greatestCommonDivisor,
          (_structure.back()->_sorter->GetOrEncodeOutput(i) << 1) ^ 1));
    }
    //    std::cout << "watchdogs.back(): " << watchdogs.back().second <<
    //    std::endl;
  }

  if (_setting->verbosity > 0) {
  }
  return watchdogs;
}

std::vector<unsigned> Cascade::GetLastAssumptions() {
  std::cout << "_fixedTareAssumption: " << _fixedTareAssumption.size()
            << std::endl;

  if (_fixedTareAssumption.empty()) {
    _fixedTareAssumption.push_back(
        (_structure.back()->_sorter->GetOrEncodeOutput(_maxPos) << 1) ^ 1);
  }
  std::cout << "_fixedTareAssumption: " << _fixedTareAssumption.size() << ":"
            << std::endl;
  for (auto i : _fixedTareAssumption) {
    std::cout << i << " ";
  }
  std::cout << std::endl;
  return _fixedTareAssumption;
}

void Cascade::PartitionSoftClauses(PartitionStrategy partitionStrategy) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (!_softClauseTreeCreated) CreateSoftClauseTree(&_softClauses, true);

  //    std::cout << "start sorting!" << std::endl;
  if (partitionStrategy == GROUPBYWEIGHT &&
      _setting->atLeastnEqualWeights > 1) {
    // test if there are at least 8 equal elements - otherwise don't sort!
    // get indices of sorted Bucket entries
    // generate Indice Vector to sort that vector!
    std::vector<unsigned> sortedSCIndices(_softClauses.size());
    std::size_t n(0);
    std::generate(std::begin(sortedSCIndices),
                  std::begin(sortedSCIndices) +
                      static_cast<unsigned>(_softClauses.size()),
                  [&] { return n++; });

    // stable sort - not changing order of SC's important for some of the
    // instances! especially spot5!
    std::stable_sort(
        std::begin(sortedSCIndices), std::end(sortedSCIndices),
        [&](std::size_t i1, std::size_t i2) {
          return (_softClauses[i2]->weight > _softClauses[i1]->weight);
        });

    //        std::cout << "before for!" << std::endl;
    //        std::cout << "softclauses.size(): " << _softClauses.size() <<
    //        std::endl;
    partitionStrategy = NOPARTITION;
    for (unsigned i = _setting->atLeastnEqualWeights; i < _softClauses.size();
         ++i) {
      //            std::cout << "i: " << i << "  weight: " <<
      //            _softClauses[sortedSCIndices[i]]->weight << "   i-n: " <<
      //            i-_setting->atLeastnEqualWeights << "  weight: " <<
      //            _softClauses[sortedSCIndices[i-_setting->atLeastnEqualWeights]]->weight
      //            << "  ((_softClauses[sortedSCIndices[i]]->weight &
      //            (_softClauses[sortedSCIndices[i]]->weight - 1)): " <<
      //            (_softClauses[sortedSCIndices[i]]->weight &
      //            (_softClauses[sortedSCIndices[i]]->weight - 1)) <<
      //            std::endl;
      if (_softClauses[sortedSCIndices[i]]->weight ==
              _softClauses[sortedSCIndices[i - _setting->atLeastnEqualWeights]]
                  ->weight &&
          // weight is not power of two!
          !((_softClauses[sortedSCIndices[i]]->weight &
             (_softClauses[sortedSCIndices[i]]->weight - 1)) == 0)) {
        partitionStrategy = GROUPBYWEIGHT;
        break;
      }
    }
    //        std::cout << "PartitionStrategy changed to: " << partitionStrategy
    //        << std::endl;
  }

  switch (partitionStrategy) {
    case NOPARTITION:
      _softClauseTree.insert(_softClauseTree.end(),
                             _processingSoftClauseTree.begin(),
                             _processingSoftClauseTree.end());
      _processingSoftClauseTree.clear();
      break;
    // both cases have the same grouping algorithm, but differ in the way they
    // connect the weights.
    case GROUPBYWEIGHTADDATLAST:
    case GROUPBYWEIGHT:
      GroupByWeight();
      _softClauseTree.insert(_softClauseTree.end(),
                             _processingSoftClauseTree.begin(),
                             _processingSoftClauseTree.end());
      _processingSoftClauseTree.clear();
      break;
    case GROUPBYBIGGESTREPEATINGENTRY:
      //        _setting->verbosity = 7;
      std::cout << "c GROUPBYBIGGESTREPEATINGENTRY" << std::endl;
      GroupByWeight();

      // actualize values as sum of both trees.
      CalculateTotalBucketEntries(&_processingSoftClauseTree, false);
      CalculateTotalBucketEntries(&_softClauseTree, true);

      DumpSCNodeStructure(&_processingSoftClauseTree, 5);

      GroupByBiggestRepeatingEntry();

      assert(_processingSoftClauseTree.empty());
      //        _setting->verbosity = 0;
      break;
    case GROUPBYCOLUMNS:

      std::cout << "c GROUPBYCOLUMNS" << std::endl;
      //        _setting->verbosity = 7;
      if (_dgpw->_dgpwSetting->featureTest == 1) {
        GroupByWeight();
      }

      // actualize values as sum of both trees.
      CalculateTotalBucketEntries(&_processingSoftClauseTree, false);
      CalculateTotalBucketEntries(&_softClauseTree, true);

      //        DumpSCNodeStructure( &_processingSoftClauseTree, 2 );

      GroupByColumns();

      assert(_processingSoftClauseTree.empty());
      //        _setting->verbosity = 0;
      break;
  }

  CalculateTotalBucketEntries(&_softClauseTree, false);

  DumpSCNodeStructure(&_softClauseTree, 5);
}

void Cascade::GroupByBiggestRepeatingEntry() {
  if (_processingSoftClauseTree.empty()) {
    return;
  }

  // std::cout << std::endl << __func__ << std::endl;

  // get indices of sorted Bucket entries
  // generate Indice Vector to sort that vector!
  std::vector<uint16_t> sortedBucketIndices(
      _processingSoftClauseTree[0]->highestBucket + 1);
  std::size_t n(0);
  std::generate(std::begin(sortedBucketIndices),
                std::begin(sortedBucketIndices) +
                    _processingSoftClauseTree[0]->highestBucket + 1,
                [&] { return n++; });

  // sort Buckets by total entries.
  std::sort(std::begin(sortedBucketIndices), std::end(sortedBucketIndices),
            [&](std::size_t i1, std::size_t i2) {
              return (_totalBucketEntries[i1] > _totalBucketEntries[i2]);
            });

  uint32_t ind(0);
  while (true) {
    if (_processingSoftClauseTree.size() == 1) {
      _softClauseTree.push_back(_processingSoftClauseTree.back());
      _processingSoftClauseTree.clear();
      std::cout << "c grouping iterations....: " << ind << std::endl;
      break;
    }

    uint16_t maxNodeIndex = GetMaxNodeIndex();

    assert(maxNodeIndex < _processingSoftClauseTree.size());

    // form a new vector with all indices of buckets containing a max Bucket
    // entry
    std::vector<uint16_t> tmpBucketIndices;
    for (auto v : sortedBucketIndices) {
      if (v <= _processingSoftClauseTree[maxNodeIndex]->highestBucket &&
          _processingSoftClauseTree[maxNodeIndex]->occursHowOftenInBucket[v] >
              0)
        tmpBucketIndices.push_back(v);
    }

    if ((ind % 500 == 0 && _setting->verbosity > 2) ||
        (ind % 250 == 0 && _setting->verbosity > 3) || _setting->verbosity > 4)
      DumpMaxNodeOverlappingsAndHeuristicValues(maxNodeIndex,
                                                &tmpBucketIndices);

    // TOBI: Is there another node with exactly the same overlapping? -> then
    // merge more nodes!
    // Maybe there is even a bigger Set in the same subset to merge first with
    // --> see bwt3cc
    int32_t nodeIndexToMergeWith =
        CalculateNodeIndexToMergeWith(maxNodeIndex, &tmpBucketIndices);

    MergeNodes(&tmpBucketIndices, nodeIndexToMergeWith, maxNodeIndex);

    if (_setting->equalWeight > 0 && (ind % _setting->equalWeight == 0))
      GroupByWeight();

    if ((ind % 1000 == 0 && _setting->verbosity > 3) ||
        _setting->verbosity > 4) {
      std::cout << "ProcessingSoftClauseTree" << std::endl;
      CalculateTotalBucketEntries(&_processingSoftClauseTree, false);
      DumpSCNodeStructure(&_processingSoftClauseTree, 6);

      std::cout << "FinalSoftClauseTree" << std::endl;
      CalculateTotalBucketEntries(&_softClauseTree, false);
      DumpSCNodeStructure(&_softClauseTree, 6);

      std::cout << "_processingSoftClauseTree.size(): "
                << _processingSoftClauseTree.size() << std::endl;
    }

    ind++;
  }
  for (std::map<uint32_t, SoftClauseNodes *>::iterator it =
           _processingPercentOffTree.begin();
       it != _processingPercentOffTree.end(); ++it) {
    _softClauseTree.push_back(it->second);
  }
  DumpSCNodeStructure(&_processingSoftClauseTree, 6);
  std::cout << "c percentOff reinsertions: "
            << _howOftenReinsertedFromProcessingPercentOffTree << std::endl;
}

Cascade::BucketOverlaps Cascade::OverlappingBuckets(
    std::vector<unsigned> *bucketIndices, std::vector<unsigned> *pSCTIndices,
    bool testOnlyLastBucket) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;
  BucketOverlaps overlaps;
  if (pSCTIndices == nullptr) {
    pSCTIndices = new std::vector<unsigned>;
    for (unsigned i = 0; i < _processingSoftClauseTree.size(); i++) {
      pSCTIndices->push_back(i);
    }
  }
  //    std::cout << "pSCTIndices: " << std::endl;
  //    std::cout << std::setw(4) << "s" << std::setw(8) << "s" << std::setw(15)
  //    << "Weight" << std::setw(3) << "|"; for (uint32_t i = 0; i <=
  //    _numberOfBuckets; ++i)
  //    {
  //        std::cout << std::setw(4) << i;
  //    }
  //    std::cout << std::endl;
  for (auto j : (*pSCTIndices)) {
    if (testOnlyLastBucket) {
      unsigned k = (*bucketIndices).back();
      if ((_processingSoftClauseTree[j]->weight & 1UL << k)) {
        assert(_processingSoftClauseTree[j]->occursHowOftenInBucket[k] > 0);
        assert(_processingSoftClauseTree[j]->weight >
               static_cast<unsigned long>(pow(2, k)));
        overlaps.noOverlaps++;
        //                std::cout << "weight: " <<
        //                _processingSoftClauseTree[j]->weight << "  PSCTindex:
        //                " << j << " position: " << k << "  pow^position: " <<
        //                std::fixed << static_cast<unsigned long>(pow(2, k)) <<
        //                std::endl;
        //                _processingSoftClauseTree[j]->dumpStructure(true, j);

        overlaps.overlappingSCTIndices.push_back(j);
      }
      continue;
    }
    for (auto i : (*bucketIndices)) {
      if (!(_processingSoftClauseTree[j]->weight & (1UL << i))) {
        break;
      } else if (i == (*bucketIndices).back()) {
        overlaps.noOverlaps++;
        //                std::cout << j << " ";
        overlaps.overlappingSCTIndices.push_back(j);
      }
    }
  }
  //    std::cout << std::endl;

  //    if (overlaps.noOverlaps > 0) {
  //        std::cout << "Compare Buckets: ";
  //        for (auto i : (*bucketIndices)) {
  //            std::cout << i << "  ";
  //        }
  //        std::cout << "OverlSCTInd: ";
  //        for (auto i : overlaps.overlappingSCTIndices) {
  //            std::cout << i << "  ";
  //        }

  //        std::cout << "   - noOverlaps: " << overlaps.noOverlaps <<
  //        std::endl;
  //    }

  overlaps.overlappingBucketIndices = (*bucketIndices);
  return overlaps;
}

unsigned Cascade::CalcExactTernaryClauseCosts(unsigned size) {
  if (size == 0) return 0;
  //    std::cout << "size: " << size << std::endl;
  unsigned a = size / 2;
  unsigned b = size - a;
  unsigned c = a * b;
  //    std::cout << "a: " << a << std::endl;
  //    std::cout << "b: " << b << std::endl;
  //    std::cout << "c: " << c << std::endl;
  unsigned d = 0;
  if (c == 0) {
    d = c;
    //        return c;
    //    } else if (a == b) {
    //        return c + 2 * CalcExactTernaryClauseCosts(a);
  } else if (a == 0) {
    d = c + CalcExactTernaryClauseCosts(b);
    //        return c + CalcExactTernaryClauseCosts(b);
  } else if (b == 0) {
    d = c + CalcExactTernaryClauseCosts(a);
    //        return c + CalcExactTernaryClauseCosts(a);
  } else {
    d = c + CalcExactTernaryClauseCosts(a) + CalcExactTernaryClauseCosts(b);
  }
  //    std::cout << "d: " << d << std::endl;
  return d;
  //    return c + CalcExactTernaryClauseCosts(a) +
  //    CalcExactTernaryClauseCosts(b);
}

void Cascade::GroupByColumns() {
  struct GreaterThan {
    bool operator()(const unsigned &left, const unsigned &right) const {
      return (left > right);
    }
  };

  std::vector<std::multimap<unsigned, BucketOverlaps, GreaterThan>>
      sortedOverlapBuckets;
  std::multimap<unsigned, BucketOverlaps, GreaterThan> sortedOverlaps;
  unsigned ind(0);

  if (_processingSoftClauseTree.empty()) {
    return;
  }

  // std::cout << std::endl << __func__ << std::endl;

  // get indices of sorted Bucket entries
  // generate Indice Vector to sort that vector!
  //    std::cout << "_numberOfBuckets: " << _numberOfBuckets << std::endl;
  std::vector<unsigned> sortedBucketIndices(_numberOfBuckets + 1);
  std::size_t n(0);
  std::generate(std::begin(sortedBucketIndices),
                std::begin(sortedBucketIndices) + _numberOfBuckets + 1,
                [&] { return n++; });
  //    std::cout << "sortedBucketIndices.size(): " <<
  //    sortedBucketIndices.size() << std::endl;

  // sort Buckets by total entries.
  std::sort(std::begin(sortedBucketIndices), std::end(sortedBucketIndices),
            [&](std::size_t i1, std::size_t i2) {
              return (_totalBucketEntries[i1] > _totalBucketEntries[i2]);
            });

  std::cout << std::endl;
  //    for (unsigned j = 0; j < sortedBucketIndices.size(); j++) {
  //        std::cout << "sortedBucketIndices[" << j << "]: " <<
  //        sortedBucketIndices[j] << std::endl;
  //    }
  //    std::cout << std::endl;

  unsigned totalSavedTernaryClauses = 0;
  unsigned totalSavedBinaryClauses = 0;
  while (true) {
    ind++;
    // calculate in which buckets the first bucket is in!
    std::vector<unsigned> tmpBucketIndices;
    for (unsigned i = 0; i <= _numberOfBuckets; i++) {
      tmpBucketIndices.push_back(i);
      BucketOverlaps bOverlaps = OverlappingBuckets(&tmpBucketIndices);
      sortedOverlaps.insert(
          std::pair<unsigned, BucketOverlaps>(bOverlaps.noOverlaps, bOverlaps));

      tmpBucketIndices.pop_back();
    }

    sortedOverlapBuckets.clear();
    sortedOverlapBuckets.push_back(sortedOverlaps);

    while (true) {
      //        sleep(1);

      sortedOverlaps.clear();
      unsigned counter = 0;
      for (auto overlapMM : sortedOverlapBuckets.back()) {
        counter++;
        if (counter > 40) {
          break;
        }
        for (unsigned bucketIndex = 0; bucketIndex <= _numberOfBuckets;
             bucketIndex++) {
          if (bucketIndex <= overlapMM.second.overlappingBucketIndices.back() ||
              find(overlapMM.second.overlappingBucketIndices.begin(),
                   overlapMM.second.overlappingBucketIndices.end(),
                   bucketIndex) !=
                  overlapMM.second.overlappingBucketIndices.end()) {
            continue;
          }
          counter++;
          tmpBucketIndices = overlapMM.second.overlappingBucketIndices;
          tmpBucketIndices.push_back(bucketIndex);
          BucketOverlaps bOverlaps = OverlappingBuckets(
              &tmpBucketIndices, &overlapMM.second.overlappingSCTIndices, true);

          unsigned minOverlaps = 3;
          if (_dgpw->_dgpwSetting->featureTest == 3) {
            minOverlaps = 2;
          } else if (_dgpw->_dgpwSetting->featureTest == 4) {
            minOverlaps = 4;
          }

          if (bOverlaps.noOverlaps < minOverlaps) {
            continue;
          }
          sortedOverlaps.insert(std::pair<unsigned, BucketOverlaps>(
              bOverlaps.noOverlaps, bOverlaps));
        }
      }

      if (sortedOverlaps.empty()) {
        break;
      }
      sortedOverlapBuckets.push_back(sortedOverlaps);
      //        for (std::multimap<unsigned, BucketOverlaps>::iterator it =
      //        sortedOverlaps.begin(); it != sortedOverlaps.end(); ++it) {
      //            std::cout << "  [" << it->first << ", (";
      //            for (auto j : it->second.overlappingBucketIndices) {
      //                std::cout << j << " ";
      //            }
      //            std::cout << "), (";
      //            for (auto j : it->second.overlappingSCTIndices) {
      //                std::cout << j << " ";
      //            }
      //            std::cout << ")]" << std::endl;
      //        }
      //        if (counter == 8) {
      //            break;
      //        }
    }
    //    std::cout << std::endl << "biggest Values: " << std::endl;
    unsigned maxTernarySavings = 0;
    unsigned maxBinarySavings = 0;
    BucketOverlaps chosenBucketOverlaps;
    for (auto overlapMMs : sortedOverlapBuckets) {
      //        std::cout << "Size: " << overlapMMs.size() << std::endl;
      //        std::cout << "  ["  << std::setw(3) << overlapMMs.begin()->first
      //        << ", "; std::cout <<
      //        overlapMMs.begin()->second.overlappingBucketIndices.size() << ":
      //        (" << std::setw(3); for (auto j :
      //        overlapMMs.begin()->second.overlappingBucketIndices) {
      //            std::cout << j << " " << std::setw(3);
      //        }
      //        std::cout << "), " << "(" << std::setw(3);
      //        for (auto j : overlapMMs.begin()->second.overlappingSCTIndices)
      //        {
      //            std::cout << j << " " << std::setw(3);
      //        }
      //        std::cout << ")]" << std::endl;
      unsigned savedTernaryClauses = 0;
      if (_dgpw->_dgpwSetting->featureTest == 2) {
        savedTernaryClauses =
            ((static_cast<unsigned>(pow(overlapMMs.begin()->first, 2)) -
              overlapMMs.begin()->first) /
             2) *
            (2 *
             (static_cast<unsigned>(
                 overlapMMs.begin()->second.overlappingBucketIndices.size() -
                 1)));
      } else if (_dgpw->_dgpwSetting->featureTest == 5) {
        savedTernaryClauses =
            ((static_cast<unsigned>(pow(overlapMMs.begin()->first, 2)) -
              overlapMMs.begin()->first) /
             2) *
            static_cast<unsigned>(pow(
                overlapMMs.begin()->second.overlappingBucketIndices.size() - 1,
                2));
      } else if (_dgpw->_dgpwSetting->featureTest == 6) {
        savedTernaryClauses =
            ((static_cast<unsigned>(pow(overlapMMs.begin()->first, 2)) -
              overlapMMs.begin()->first) /
             2) *
            static_cast<unsigned>(pow(
                overlapMMs.begin()->second.overlappingBucketIndices.size() - 1,
                3));
      } else {
        savedTernaryClauses =
            (static_cast<unsigned>(pow(overlapMMs.begin()->first, 2)) -
             overlapMMs.begin()->first) /
            2 *
            static_cast<unsigned>(
                overlapMMs.begin()->second.overlappingBucketIndices.size() - 1);
      }
      unsigned savedBinaryClauses =
          (static_cast<unsigned>(log2(overlapMMs.begin()->first)) *
           overlapMMs.begin()->first *
           static_cast<unsigned>(
               overlapMMs.begin()->second.overlappingBucketIndices.size() - 1));
      //        unsigned savedVariables = savedBinaryClauses;
      if (savedTernaryClauses > maxTernarySavings) {
        //            std::cout << "This bucket is chosen" << std::endl;
        //            for (auto j :
        //            overlapMMs.begin()->second.overlappingBucketIndices) {
        ////                std::cout << j << " " << std::setw(3);
        //            }
        maxTernarySavings = savedTernaryClauses;
        maxBinarySavings = savedBinaryClauses;
        chosenBucketOverlaps = overlapMMs.begin()->second;
      }
      //        std::cout << "Saved ternary clauses: " << std::setw(5) <<
      //        savedTernaryClauses; std::cout << "  Saved binary clauses: " <<
      //        std::setw(5) << savedBinaryClauses << std::endl; std::cout << "
      //        Saved variables: " << std::setw(5) << savedVariables <<
      //        std::endl << std::endl; savedTernaryClauses = 0;
      //        savedTernaryClauses =
      //        CalcExactTernaryClauseCosts(overlapMMs.begin()->first) *
      //        static_cast<unsigned>(overlapMMs.begin()->second.overlappingBucketIndices.size()
      //        - 1); savedBinaryClauses =
      //        (static_cast<unsigned>(log2(overlapMMs.begin()->first)) *
      //        overlapMMs.begin()->first *
      //        static_cast<unsigned>(overlapMMs.begin()->second.overlappingBucketIndices.size()
      //        - 1)); savedVariables = savedBinaryClauses; std::cout << "Saved
      //        ternary clauses: " << std::setw(5) << savedTernaryClauses;
      //        std::cout << "   Saved binary clauses: " << std::setw(5) <<
      //        savedBinaryClauses; std::cout << "   Saved variables: " <<
      //        std::setw(5) << savedVariables << std::endl << std::endl;
    }
    totalSavedTernaryClauses += maxTernarySavings;
    totalSavedBinaryClauses += maxBinarySavings;
    //    std::cout << "maxSavings: " << maxTernarySavings << std::endl;
    if (maxTernarySavings == 0) {
      break;
    }

    if (_processingSoftClauseTree.size() == 1) {
      _softClauseTree.push_back(_processingSoftClauseTree.back());
      _processingSoftClauseTree.clear();
      break;
    }

    //    std::cout << "Chosen Bucket overlaps!!!" << std::endl;
    //    std::cout << "  ["  << std::setw(3) << chosenBucketOverlaps.noOverlaps
    //    << ", "; std::cout <<
    //    chosenBucketOverlaps.overlappingBucketIndices.size() << ": (" <<
    //    std::setw(3); for (auto j :
    //    chosenBucketOverlaps.overlappingBucketIndices) {
    //        std::cout << j << " " << std::setw(3);
    //    }
    //    std::cout << "), " << "(" << std::setw(3);
    //    for (auto j : chosenBucketOverlaps.overlappingSCTIndices) {
    //        std::cout << j << " " << std::setw(3);
    //    }
    //    std::cout << ")]" << std::endl;

    MergeMultipleNodes(&chosenBucketOverlaps);

    //    if ((ind % 1000 == 0 && _setting->verbosity > 3)  ||
    //    _setting->verbosity > 4)
    //    {
    //        std::cout << "ProcessingSoftClauseTree" << std::endl;
    //        CalculateTotalBucketEntries( &_processingSoftClauseTree, false);
    //        DumpSCNodeStructure(&_processingSoftClauseTree, 3);

    //        std::cout << "FinalSoftClauseTree" << std::endl;
    //        CalculateTotalBucketEntries( &_softClauseTree, false);
    //        DumpSCNodeStructure(&_softClauseTree, 3);
    //    }

    //    if (ind == 7) {
    //        break;
    //    }
  }

  std::cout << "c adder caching on all rows of the soft clause tree in " << ind
            << " rounds!" << std::endl;

  // append all elements of the _processingSoftClauseTree to _softClauseTree
  _softClauseTree.insert(std::end(_softClauseTree),
                         std::begin(_processingSoftClauseTree),
                         std::end(_processingSoftClauseTree));
  _processingSoftClauseTree.clear();

  if ((ind % 1000 == 0 && _setting->verbosity > 3) || _setting->verbosity > 4) {
    std::cout << "FinalSoftClauseTree" << std::endl;
    DumpSCNodeStructure(&_softClauseTree, 6);
  }

  std::cout << "c TotalSavedTernaryClause: " << totalSavedTernaryClauses
            << std::endl;
  std::cout << "c TotalSavedBinaryClauses: " << totalSavedBinaryClauses
            << std::endl;
  std::cout << "c grouping iterations....: " << ind << std::endl;
  //    exit(0);
}

void Cascade::MergeMultipleNodes(BucketOverlaps *overlaps) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;
  assert(overlaps->overlappingBucketIndices.size() > 1);
  assert(overlaps->overlappingSCTIndices.size() > 1);

  // calc new weight for the two subweights to merge
  uint64_t newOverlappingWeight(0);
  for (auto v : overlaps->overlappingBucketIndices) {
    newOverlappingWeight += pow(_base, v);
    //        std::cout << "v: " << v << "   pow(base, v): " << std::fixed <<
    //        static_cast<unsigned long>(pow(_base,v)) << "
    //        newOverlappingWeight: " << newOverlappingWeight << std::endl;
  }
  //    std::cout << "newOverlappingWeight: " << newOverlappingWeight <<
  //    std::endl;
  Bucket *tmpBucket = new Bucket(_dgpw, this);
  for (auto w : overlaps->overlappingSCTIndices) {
    tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[w]);
  }
  SoftClauseNodes *sCNode =
      new SoftClauseNodes(tmpBucket, newOverlappingWeight, _base);
  //    sCNode->dumpStructure(true);
  // is at least in two buckets - otherwise no merge
  _processingSoftClauseTree.push_back(sCNode);

  std::vector<unsigned> eraseList;

  // erase the old SoftClauseNodes.
  for (auto w : overlaps->overlappingSCTIndices) {
    //        std::cout << "_processingSoftClauseTree[w]->weight: " <<
    //        _processingSoftClauseTree[w]->weight << "  index: " << w <<
    //        std::endl; std::cout << "newOverlappingWeight: " <<
    //        newOverlappingWeight << "  index: " << w << std::endl; std::cout
    //        << "newWeight: " << newWeight << "  index: " << w << std::endl <<
    //        std::endl; std::cout << "pSCT[w]->weight: " << std::setw(20) <<
    //        _processingSoftClauseTree[w]->weight; std::cout << "
    //        newOverlappingWeight: "  << std::setw(20) << newOverlappingWeight
    //        << "  index: "  << std::setw(5) << w << std::endl;
    assert(_processingSoftClauseTree[w]->weight >= newOverlappingWeight);
    uint64_t newWeight =
        _processingSoftClauseTree[w]->weight - newOverlappingWeight;
    //        std::cout << "   newWeight: "  << std::setw(20) << newWeight <<
    //        std::endl;
    if (newWeight == 0) {
      eraseList.push_back(w);
    } else {
      _processingSoftClauseTree[w]->setWeight(newWeight);
      if (_processingSoftClauseTree[w]->inHowManyBuckets == 1) {
        _softClauseTree.push_back(_processingSoftClauseTree[w]);
        eraseList.push_back(w);
      }
    }
  }

  for (unsigned eraseInd = static_cast<unsigned>(eraseList.size());
       eraseInd > 0; eraseInd--) {
    //        std::cout << "EraseIndex: " << eraseList[eraseInd - 1] <<
    //        std::endl; _processingSoftClauseTree[eraseList[eraseInd -
    //        1]]->dumpStructure(true);
    _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() +
                                    eraseList[eraseInd - 1]);
  }

  if (_setting->verbosity < 4) return;

  //    std::cout << std::setw(30) << "new Weight: " << newOverlappingWeight <<
  //    std::endl;
}

void Cascade::GroupByWeight() {
  std::vector<SoftClauseNodes *> tmpSCTree;
  // sort Buckets by total entries.
  // use stable sort only if there are identical weights!
  std::stable_sort(std::begin(_processingSoftClauseTree),
                   std::end(_processingSoftClauseTree),
                   [&](SoftClauseNodes *sCN1, SoftClauseNodes *sCN2) {
                     return (sCN1->weight > sCN2->weight);
                   });

  for (uint32_t i = 0; i < _processingSoftClauseTree.size(); ++i) {
    if (i + 1 == _processingSoftClauseTree.size() ||
        _processingSoftClauseTree[i + 1]->weight !=
            _processingSoftClauseTree[i]->weight) {
      tmpSCTree.push_back(_processingSoftClauseTree[i]);
      continue;
    }
    Bucket *tmpBucket = new Bucket(_dgpw, this);
    while (i + 1 < _processingSoftClauseTree.size() &&
           _processingSoftClauseTree[i + 1]->weight ==
               _processingSoftClauseTree[i]->weight) {
      tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[i]);
      i++;
    }
    tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[i]);
    SoftClauseNodes *sCNode = new SoftClauseNodes(
        tmpBucket, _processingSoftClauseTree[i]->weight, _base);
    tmpSCTree.push_back(sCNode);
  }
  _processingSoftClauseTree = tmpSCTree;
}

uint16_t Cascade::GetMaxNodeIndex() {
  // generate Indice Vector to get highest index from!
  std::vector<uint16_t> maxSCTIndices(_processingSoftClauseTree.size());
  std::size_t m(0);
  std::generate(std::begin(maxSCTIndices),
                std::begin(maxSCTIndices) + _processingSoftClauseTree.size(),
                [&] { return m++; });

  // Get index of SoftClauseTree with most Bucket entries!
  // If 2 SCN have same # entries, take the one with more occurrences.
  return *std::max_element(
      std::begin(maxSCTIndices), std::end(maxSCTIndices),
      [&](std::size_t i1, std::size_t i2) {
        return ((_processingSoftClauseTree[i1]->inHowManyBuckets <
                 _processingSoftClauseTree[i2]->inHowManyBuckets) ||
                ((_processingSoftClauseTree[i1]->inHowManyBuckets ==
                  _processingSoftClauseTree[i2]->inHowManyBuckets) &&
                 (_processingSoftClauseTree[i1]->GetOccurrences() *
                      _processingSoftClauseTree[i1]->size() <
                  _processingSoftClauseTree[i2]->GetOccurrences() *
                      _processingSoftClauseTree[i2]->size())));
      });

  // First Idea - sort whole tree...
  // sort whole SoftClauseTree by entries and if they are equal by occurences.
  // maybe to look later on only to the top elements. Otherwise it is sufficient
  // to find in O(n) the biggest element.
  // std::sort(_processingSoftClauseTree.begin(),
  // _processingSoftClauseTree.end(),
  // SoftClauseNodes::SortByEntriesOccurrences);
}

void Cascade::DumpMaxNodeOverlappingsAndHeuristicValues(
    uint16_t maxNodeIndex, std::vector<uint16_t> *tmpBucketIndices) {
  if (_setting->verbosity < 1) return;

  std::cout << "_processingSoftClauseTree.size: "
            << _processingSoftClauseTree.size() << std::endl;
  std::cout << "maxBucketEntryIndex: " << maxNodeIndex << "  maxBucketEntries: "
            << _processingSoftClauseTree[maxNodeIndex]->inHowManyBuckets
            << std::endl
            << std::endl;
  for (auto v : *tmpBucketIndices) std::cout << std::setw(4) << v;
  std::cout << std::endl;
  for (uint64_t ind = 0; ind < (*tmpBucketIndices).size(); ind++)
    std::cout << std::setw(4) << "----";
  std::cout << "--------------------------------------------------"
            << std::endl;

  for (auto v : *tmpBucketIndices) {
    if (_processingSoftClauseTree[maxNodeIndex]->occursHowOftenInBucket[v] >
            0 &&
        v < _processingSoftClauseTree[maxNodeIndex]->highestBucket + 1)
      std::cout << std::setw(4)
                << _processingSoftClauseTree[maxNodeIndex]->size();
    else
      std::cout << std::setw(4) << "";
  }
  std::cout << "   |" << std::setw(7) << "0"
            << " |" << std::setw(7) << "1"
            << " |" << std::setw(7) << "2"
            << " |" << std::setw(7) << "3"
            << " |" << std::setw(7) << "4"
            << " |" << std::setw(7) << "5"
            << " |" << std::setw(7) << "6"
            << " |  <-- heuristics" << std::endl;

  for (uint64_t ind = 0; ind < (*tmpBucketIndices).size(); ind++)
    std::cout << std::setw(4) << "----";
  std::cout << "--------------------------------------------------"
            << std::endl;

  int32_t maxCosts(0);
  int32_t estimatedMaxCostIndex(0);

  for (int32_t j = 0; j < (int)_processingSoftClauseTree.size(); j++) {
    if (j == maxNodeIndex) continue;

    std::vector<SoftClauseNodes *> NodesToMerge = {
        _processingSoftClauseTree[maxNodeIndex], _processingSoftClauseTree[j]};
    int32_t usedCosts = GetBenefitOfMergingNodes(NodesToMerge, tmpBucketIndices,
                                                 true, _setting->groupHeuristic,
                                                 _setting->percentOff);

    if (usedCosts == 0) continue;

    if (maxCosts < usedCosts) {
      maxCosts = usedCosts;
      estimatedMaxCostIndex = j;
    }
  }

  for (uint64_t ind = 0; ind < (*tmpBucketIndices).size(); ind++)
    std::cout << std::setw(4) << "----";
  std::cout << "--------------------------------------------------"
            << std::endl;

  for (auto v : *tmpBucketIndices)
    std::cout << std::setw(4) << _totalBucketEntries[v];
  std::cout << std::endl;
  std::cout << std::endl
            << "MaxCosts: " << maxCosts << "  Index: " << estimatedMaxCostIndex
            << std::endl;
}

void Cascade::DumpModelOfTares(uint16_t verbosity) {
  if (_setting->verbosity < verbosity || _dgpw->GetLastResult() != 10) return;

  std::cout << __PRETTY_FUNCTION__ << std::endl;

  // if there wasn't a SAT call making the result bigger after constructing the
  // whole cascade, then the tares are 0.
  //    if (_dgpw->_lastModel[_structure[0]->_tares[0]] == 0)
  //        return;
  // output the variables of the Tare T
  //    _dgpw->Solve();
  std::cout << "Model of Tares: (n...0): ";
  for (int bucketInd = (_structure.size() - 1); bucketInd >= 0; --bucketInd) {
    std::cout << "(";
    for (int tareInd =
             static_cast<int>(_structure[bucketInd]->_tares.size() - 1);
         tareInd >= 0; --tareInd) {
      std::cout << _dgpw->Model(_structure[bucketInd]->_tares[tareInd]);
      if (tareInd != 0) {
        std::cout << ", ";
      }
    }
    if (bucketInd != 0) {
      std::cout << "), ";
    } else {
      std::cout << ")";
    }
  }
  std::cout << std::endl;
}

void Cascade::DumpBucketSolveInformation(uint32_t actualPos, bool _isLastBucket,
                                         uint16_t verbosity) {
  if (_setting->verbosity < verbosity) return;

  if (_isLastBucket) {
    _estimatedWeightBoundaries[0] = _highestBucketMultiplicator * actualPos;
    _estimatedWeightBoundaries[1] =
        _highestBucketMultiplicator * _structure.back()->size();
    std::cout << std::setw(52) << "weight boundaries: ( "
              << _estimatedWeightBoundaries[0] << " / "
              << _estimatedWeightBoundaries[1] << " )" << std::endl;
  }
  std::cout << std::setw(52) << "satisfied weights: ( " << _satWeight << " / "
            << _sumOfSoftWeights << " )" << std::endl;
  std::cout << std::setw(50) << "actualPosition: " << actualPos << std::endl;
  std::cout << "---------------------------------------------------------------"
               "------------------"
               "----------"
            << std::endl;
}

int32_t Cascade::CalculateNodeIndexToMergeWith(
    uint16_t maxNodeIndex, std::vector<uint16_t> *tmpBucketIndices) {
  int32_t maxCosts(0);
  // highest possible number to indicate that there is no clause node to merge
  // with!
  int32_t estimatedMaxCostIndex(-1);

  int32_t sumOfSizesOfPercentageFails(0);
  for (uint32_t j = 0; j < _processingSoftClauseTree.size(); j++) {
    if (j == maxNodeIndex) continue;

    std::vector<SoftClauseNodes *> NodesToMerge = {
        _processingSoftClauseTree[maxNodeIndex], _processingSoftClauseTree[j]};
    int32_t usedCosts = GetBenefitOfMergingNodes(
        NodesToMerge, tmpBucketIndices, false, _setting->groupHeuristic,
        _setting->percentOff);
    if (usedCosts == 0) {
      continue;
    } else if (usedCosts < 0) {
      sumOfSizesOfPercentageFails -= usedCosts;
    }

    if (maxCosts < usedCosts) {
      maxCosts = usedCosts;
      estimatedMaxCostIndex = j;
    }
  }
  // If no merge is possible
  // check if there is the possibility that the node can be merged later on
  if (estimatedMaxCostIndex == -1 && sumOfSizesOfPercentageFails > 0 &&
      _setting->percentOffReinsert) {
    if ((double)_processingSoftClauseTree[maxNodeIndex]->size() *
            (double)((100 - _setting->percentOff) / 100) <=
        (double)sumOfSizesOfPercentageFails) {
      return -2;
    }
  }

  // std::cout << "sumOfSizesOfPercentageFails: " << sumOfSizesOfPercentageFails
  // << std::endl;
  return estimatedMaxCostIndex;
}

int32_t Cascade::GetBenefitOfMergingNodes(
    std::vector<SoftClauseNodes *> NodesToMerge,
    std::vector<uint16_t> *tmpBucketIndices, bool dump, uint32_t heuristic,
    int32_t minPercentOff) {
  int32_t maxSize = (NodesToMerge[0]->size() > NodesToMerge[1]->size())
                        ? NodesToMerge[0]->size()
                        : NodesToMerge[1]->size();
  int32_t minSize = (NodesToMerge[0]->size() > NodesToMerge[1]->size())
                        ? NodesToMerge[1]->size()
                        : NodesToMerge[0]->size();
  //    std::cout << "maxSize: " << maxSize << "  minSize: " << minSize <<
  //    std::endl;
  int calcPercentOff = 100 - (int)(((double)minSize / (double)maxSize) * 100);
  int32_t binaryClauses = NodesToMerge[1]->size() + NodesToMerge[0]->size();
  int32_t ternaryClauses = NodesToMerge[1]->size() * NodesToMerge[0]->size();
  int32_t clauseCosts = 2 * ternaryClauses + binaryClauses;
  int32_t usedCosts(0);
  int32_t howManyBucketsOverlapping(0);
  int32_t sumOfBucketSizes(1);
  int32_t howOftenBucketsOverlapping(0);
  int32_t bucketSizesFactor(1);
  for (auto v : *tmpBucketIndices) {
    if (!(NodesToMerge[1]->occursHowOftenInBucket[v] > 0 &&
          v < NodesToMerge[1]->highestBucket + 1))
      continue;

    bucketSizesFactor +=
        static_cast<uint32_t>(pow(_totalBucketEntries[v], 1.7) * 0.1);

    howManyBucketsOverlapping += NodesToMerge[1]->size();
    sumOfBucketSizes += _totalBucketEntries[v];
    howOftenBucketsOverlapping++;
  }
  if (howOftenBucketsOverlapping <= 1) return 0;

  switch (heuristic) {
    // Standard combination of other heuristics
    case 0:
      usedCosts =
          bucketSizesFactor *
          static_cast<int32_t>(pow((howOftenBucketsOverlapping - 1), 1.7) *
                               pow((clauseCosts), 0.5));
      break;
    // number of reduced (ternary clauses * 2 + binary clauses)
    case 1:
      usedCosts = (howOftenBucketsOverlapping - 1) * clauseCosts;
      //        std::cout << "howOftenBucketsOverlapping: " <<
      //        howOftenBucketsOverlapping << "   clauseCosts: " << clauseCosts
      //        << "   ternaryClauses : " << ternaryClauses << "  binaryClauses:
      //        " << binaryClauses << std::endl;
      break;
    // sum of bucket sizes the SC occurs
    case 2:
      usedCosts = sumOfBucketSizes;
      break;
    // closest size in percentage
    case 3:
      usedCosts = 100 - calcPercentOff;
      break;
    // the one with least occurences in other buckets
    case 4:
      usedCosts =
          (int)_totalBucketEntries.size() -
          ((int)NodesToMerge[1]->inHowManyBuckets - howOftenBucketsOverlapping);
      break;
    // how many merges are possible after this merge
    case 5:
      usedCosts =
          CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 0);
      break;
    // greatest depth of submerges (till depth 3, then adds possible merges of
    // depth 4) can be made more effective - but for the ones with big numbers
    // it is way to much work!
    case 6:
      usedCosts =
          CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 1);
      break;
    //
    case 7:
      // usedCosts = CalculateNumberOfPossibleSubmerges(tmpBucketIndices,
      // NodesToMerge[1]);
      break;
  }

  if (minPercentOff != 100) {
    // at least one difference - should be always possible, or the calculated
    // percentage.
    if (maxSize - minSize != 1 && calcPercentOff > minPercentOff) {
      // if NodesToMerge[0]->size() is bigger, then give back the negative value
      // to calculate if it is theoretically possible to reach that value again
      // if no merge at all is possible.
      if (NodesToMerge[0]->size() > NodesToMerge[1]->size()) {
        usedCosts = -NodesToMerge[1]->size();
      } else {
        usedCosts = 0;
      }
    }
  }

  if (!dump) return usedCosts;

  int32_t heuristic0 =
      bucketSizesFactor *
      static_cast<int32_t>(pow((howOftenBucketsOverlapping - 1), 1.7) *
                           pow((clauseCosts), 0.5));
  int32_t heuristic1 = (howOftenBucketsOverlapping - 1) * clauseCosts;
  //            CalculateHowManySCNHaveTheSameOverlap(tmpBucketIndices,
  //            NodesToMerge);
  int32_t heuristic2 = sumOfBucketSizes;
  int32_t heuristic3 = calcPercentOff;
  int32_t heuristic4 =
      (int)_totalBucketEntries.size() -
      ((int)NodesToMerge[1]->inHowManyBuckets - howOftenBucketsOverlapping);
  ;
  int32_t heuristic5 =
      CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 0);
  int32_t heuristic6 =
      CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 1);
  int32_t heuristic7 =
      CalculateHowManySCNHaveTheSameOverlap(tmpBucketIndices, NodesToMerge);
  //    int32_t heuristic7 =
  //    CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 1);

  for (auto v : *tmpBucketIndices) {
    if (NodesToMerge[1]->occursHowOftenInBucket[v] > 0 &&
        v < NodesToMerge[1]->highestBucket + 1)
      std::cout << std::setw(4) << NodesToMerge[1]->size();
    else
      std::cout << std::setw(4) << "";
  }
  std::cout << "  ";
  std::cout << " |" << std::setw(7) << heuristic0;
  std::cout << " |" << std::setw(7) << heuristic1;
  std::cout << " |" << std::setw(7) << heuristic2;
  std::cout << " |" << std::setw(7) << heuristic3;
  std::cout << " |" << std::setw(7) << heuristic4;
  std::cout << " |" << std::setw(7) << heuristic5;
  std::cout << " |" << std::setw(7) << heuristic6;
  std::cout << " |" << std::setw(7) << heuristic7;

  if (minPercentOff != 100) {
    std::cout << " |" << std::setw(7) << calcPercentOff;

    if (usedCosts < 0)
      std::cout << "  <   " << minPercentOff << " <- too big difference";
    else
      std::cout << "  >=  " << minPercentOff;

    std::cout << " |  usedCosts: " << usedCosts << std::endl;
  } else {
    std::cout << std::endl;
  }

  return usedCosts;
}

int32_t Cascade::CalculateHowManySCNHaveTheSameOverlap(
    std::vector<uint16_t> *tmpBucketIndices,
    std::vector<SoftClauseNodes *> NodesToMerge) {
  //    std::cout << std::endl << __PRETTY_FUNCTION__ << std::endl;

  std::vector<uint16_t> overlappingIndices;
  //    std::cout << "Overlapping indices are: ";
  for (auto v : *tmpBucketIndices) {
    if (NodesToMerge[1]->occursHowOftenInBucket[v] >= 1) {
      //            std::cout << v << " ";
      overlappingIndices.push_back(v);
    }
  }
  //    std::cout << std::endl;
  //    std::cout << "They have " << overlappingIndices.size() << " overlapping
  //    indices." << std::endl;

  std::vector<SoftClauseNodes *> overlappingSCN;
  for (auto softClauseNodes : _processingSoftClauseTree) {
    size_t overlaps = 0;
    //        std::cout << std::endl;
    for (auto v : overlappingIndices) {
      if (softClauseNodes->occursHowOftenInBucket[v] >= 1 &&
          softClauseNodes->highestBucket >= v) {
        overlaps++;
        //                std::cout << "v: " << v << "  howOften: " <<
        //                softClauseNodes->occursHowOftenInBucket[v] << "
        //                overlaps: " << overlaps << std::endl;
      } else {
        break;
      }
    }
    if (overlaps == overlappingIndices.size()) {
      overlappingSCN.push_back(softClauseNodes);
      //            softClauseNodes->dumpStructure(0);
    }
  }
  //    std::cout << "Overlapping SoftClauseNodes: " << overlappingSCN.size() <<
  //    std::endl;

  return static_cast<int32_t>(overlappingSCN.size());

  //    std::cout << std::endl;
  for (auto w : NodesToMerge) {
    std::cout << "DumpStructure: " << std::endl;
    w->dumpStructure(1);
    std::cout << "Done!" << std::endl;
  }
  std::cout << std::endl;
  std::cout << std::endl;
  exit(1);
}

int32_t Cascade::CalculateNumberOfPossibleSubmerges(
    std::vector<uint16_t> *tmpBucketIndices,
    std::vector<SoftClauseNodes *> NodesToMerge, int32_t depth) {
  int32_t incomingDepth(depth);
  int32_t possibleSubmerges(0);
  int32_t highestDepth(0);
  std::vector<uint16_t> lastNodeIndices;
  // std::cout << "SNI: ";

  for (auto v : *tmpBucketIndices) {
    if (!(NodesToMerge.back()->occursHowOftenInBucket[v] > 0 &&
          v < NodesToMerge.back()->highestBucket + 1))
      continue;
    lastNodeIndices.push_back(v);
    // std::cout << v << ", ";
  }

  for (uint32_t j = 0; j < _processingSoftClauseTree.size(); j++) {
    // checking if j'th element is one of the yet merged ones.
    if (std::find(NodesToMerge.begin(), NodesToMerge.end(),
                  _processingSoftClauseTree[j]) != NodesToMerge.end()) {
      // std::cout << "sizeofNodesToMerge: " << NodesToMerge.size() << "  Index:
      // " << j << std::endl;
      continue;
    }
    // std::cout << "j: " << j << std::endl;
    // if (_processingSoftClauseTree[j] == NodesToMerge.back())
    //    continue;
    int32_t occurencesInNode(0);
    for (auto v : lastNodeIndices) {
      if (_processingSoftClauseTree[j]->highestBucket < v) continue;
      // std::cout << "j: " << j << " v: " << v << std::endl;
      if (_processingSoftClauseTree[j]->occursHowOftenInBucket[v] > 0)
        occurencesInNode++;
    }
    if (depth > 0 && occurencesInNode > 1) {
      NodesToMerge.push_back(_processingSoftClauseTree[j]);
      // here we are at depth 3! - go one depth deeper and Calculate then the
      // number of possible submerges!
      if (incomingDepth == 2)
        return (incomingDepth + 1 +
                CalculateNumberOfPossibleSubmerges(&lastNodeIndices,
                                                   NodesToMerge, 0));

      int32_t newDepth = CalculateNumberOfPossibleSubmerges(
          &lastNodeIndices, NodesToMerge, incomingDepth + 1);
      // std::cout << "newDepth: " << newDepth << std::endl;
      if (newDepth > highestDepth) highestDepth = newDepth;
    } else if (occurencesInNode > 1) {
      // std::cout << "Index: " << j << "  oIN: " << occurencesInNode <<
      // std::endl;
      possibleSubmerges++;
    }
  }
  if (highestDepth == 0) highestDepth = incomingDepth;
  if (depth > 0)
    return highestDepth;
  else
    return possibleSubmerges;
}

void Cascade::MergeNodes(std::vector<uint16_t> *tmpBucketIndices,
                         int32_t nodeIndexToMergeWith, int32_t maxNodeIndex) {
  // Exception - there is no node to merge with!!
  if (nodeIndexToMergeWith == -1) {
    _softClauseTree.push_back(_processingSoftClauseTree[maxNodeIndex]);
    _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() +
                                    maxNodeIndex);
    return;
  }

  // Exception - no node to merge with, but with other merges it can become
  // possible again.
  else if (nodeIndexToMergeWith == -2) {
    _processingPercentOffTree.insert(
        std::make_pair(_processingSoftClauseTree[maxNodeIndex]->size(),
                       _processingSoftClauseTree[maxNodeIndex]));
    _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() +
                                    maxNodeIndex);
    return;
  }

  // calc new weight for the two subweights to merge
  uint64_t newWeightforBoth(0);
  for (auto v : *tmpBucketIndices) {
    if (_processingSoftClauseTree[nodeIndexToMergeWith]
            ->occursHowOftenInBucket[v] &&
        v < _processingSoftClauseTree[nodeIndexToMergeWith]->highestBucket +
                1) {
      newWeightforBoth += pow(_base, v);
    }
  }

  Bucket *tmpBucket = new Bucket(_dgpw, this);
  tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[maxNodeIndex]);
  tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[nodeIndexToMergeWith]);
  SoftClauseNodes *sCNode =
      new SoftClauseNodes(tmpBucket, newWeightforBoth, _base);
  // at least in two buckets - otherwise no merge
  _processingSoftClauseTree.push_back(sCNode);

  // erase the old SoftClauseNodes.
  uint16_t nodeIndex =
      maxNodeIndex > nodeIndexToMergeWith ? maxNodeIndex : nodeIndexToMergeWith;

  for (uint16_t i = 0; i < 2; i++) {
    if (i == 1)
      nodeIndex = maxNodeIndex > nodeIndexToMergeWith ? nodeIndexToMergeWith
                                                      : maxNodeIndex;

    if (_processingSoftClauseTree[nodeIndex]->weight - newWeightforBoth == 0) {
      _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() +
                                      nodeIndex);
    } else {
      _processingSoftClauseTree[nodeIndex]->setWeight(
          _processingSoftClauseTree[nodeIndex]->weight - newWeightforBoth);
      // only in one bucket - then move it to the from processingSCT to the SCT
      if (_processingSoftClauseTree[nodeIndex]->inHowManyBuckets == 1) {
        _softClauseTree.push_back(_processingSoftClauseTree[nodeIndex]);
        _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() +
                                        nodeIndex);
      }
    }
  }

  // if the new weight is bigger than a weight from the _laterToMergeMultimap
  // move the corresponding SCN to the _processing SCT.
  for (std::multimap<uint32_t, SoftClauseNodes *>::iterator it =
           _processingPercentOffTree.begin();
       it != _processingPercentOffTree.end(); ++it) {
    if (!AtLeastTwoBucketsInCommon(it->second,
                                   _processingSoftClauseTree.back()))
      continue;

    int32_t maxSize = (it->first > _processingSoftClauseTree.back()->size())
                          ? it->first
                          : _processingSoftClauseTree.back()->size();
    int32_t minSize = (it->first > _processingSoftClauseTree.back()->size())
                          ? _processingSoftClauseTree.back()->size()
                          : it->first;
    uint32_t calcPercentOff =
        100 - (uint32_t)(((double)minSize / (double)maxSize) * 100);

    // Is there a chance of merging again?
    if (calcPercentOff <= _setting->percentOff) {
      if (_setting->verbosity > 3) {
        std::cout << "_processingPercentOffTree.size(); "
                  << _processingPercentOffTree.size() << std::endl;
        std::cout << "percentOff: " << calcPercentOff
                  << "   _setting->_percentOff: " << _setting->percentOff
                  << "   minSize: " << minSize << "  maxSize: " << maxSize
                  << std::endl;
      }

      _processingSoftClauseTree.push_back(it->second);
      _processingPercentOffTree.erase(it);
      _howOftenReinsertedFromProcessingPercentOffTree++;
    }
    // because of sorted Multimap - there is no chance that calcPercentOff gets
    // biggerEqual than _dgpw->groupPercentOff again.
    else if (it->first > _processingSoftClauseTree.back()->size()) {
      break;
    }
  }

  if (_setting->verbosity < 4) return;

  std::cout << std::setw(30) << "new Weight: " << newWeightforBoth << std::endl;
}

void Cascade::FillBuckets() {
  // create the bucket structure.
  for (uint16_t ind = 0; ind <= _numberOfBuckets; ind++) {
    Bucket *tmpBucket = new Bucket(_dgpw, this, ind);
    _structure.push_back(tmpBucket);
  }

  // fill each Bucket according to the _softClauseTree structure
  for (auto sCNode : _softClauseTree) {
    for (uint16_t ind = 0; ind <= sCNode->highestBucket; ind++) {
      //            if (sCNode->occursHowOftenInBucket[ind] == 0)
      //                continue;
      //
      for (uint32_t howOften = 0;
           howOften < sCNode->occursHowOftenInBucket[ind]; howOften++) {
        _structure[ind]->AddSoftClauseNode(sCNode);
      }
    }
  }

  _highestBucketMultiplicator =
      static_cast<uint64_t>(pow(_base, _numberOfBuckets));

  if (_setting->verbosity < 4) return;

  std::cout << "Buckets are filled with SoftClauseNodes!" << std::endl;
  DumpBucketStructure(false, 5);
}

void Cascade::AddTaresToBuckets() {
  assert(_structure.size() > 0);

  //    uint32_t addNoTareToLastBucket = ( _dgpw->_cascadeDivider > 0 ||
  //    _onlyByTares ) ? 0 : 1;
  uint32_t addNoTareToLastBucket = 1;

  if (_setting->verbosity > 3) {
    std::cout << "caddNoTareToLastBucket: " << addNoTareToLastBucket
              << std::endl;
    std::cout << "cTare(s): ";
  }
  // add at the last position _base many variables T to each trigger except the
  // last one.
  for (uint32_t i = 0; i < _structure.size() - addNoTareToLastBucket; i++) {
    // assert( _dgpw->_sorterTree[i].size() <= 1 );
    if (!_structure.empty()) {
      for (uint32_t j = 0; j < _setting->base - 1; j++) AddTare(i);
    }
  }

  //    std::cout << "Tares: ";
  //    for (uint32_t i = 0; i < _structure.size() - addNoTareToLastBucket; i++)
  //    {
  //        //assert( _dgpw->_sorterTree[i].size() <= 1 );
  //        if( !_structure.empty() )
  //        {
  //            std::cout << _structure[i]->_tares[0] << ", ";
  //        }
  //    }
  //    std::cout << std::endl;

  if (_setting->verbosity < 2) return;

  std::cout << std::endl << "c Tares are added to Structure!" << std::endl;
}

void Cascade::EncodeTopBuckets() {
  TimeMeasurement timeEncoding(&_dgpw->_timeVariables->encoding, true);
  for (auto bucket : _structure) {
    bucket->CalculateNumberOfClauses(true, true, false);

    if (_setting->partitionStrategy == GROUPBYWEIGHTADDATLAST) {
      bucket->EncodeTopAddAtLast();
    } else {
      bucket->EncodeTop();
    }

    bucket->CalculateNumberOfClauses(true, false, true);
  }

  if (_setting->verbosity < 1) return;
  if (_setting->verbosity > 3)
    std::cout << std::endl << "Top Buckets are encoded!" << std::endl;
}

void Cascade::EncodeBottomBuckets() {
  TimeMeasurement timeEncoding(&_dgpw->_timeVariables->encoding, true);
  for (uint16_t ind = 1; ind < _structure.size(); ind++) {
    _structure[ind]->CalculateNumberOfClauses(false, true, false);

    _structure[ind]->MergeSorterWith(
        _structure[ind - 1]->GetEveryNthOutput(_base));

    _structure[ind]->CalculateNumberOfClauses(false, false, true);
  }

  _structure.back()->_isLastBucket = true;

  if (_setting->verbosity < 1) return;

  DumpNumberOfBucketsAndClauses();

  if (_setting->verbosity > 3)
    std::cout << "Bottom Buckets are encoded!" << std::endl << std::endl;
}

void Cascade::UnionBucketsIntoLast() {
  // Push all Buckets to one final Bucket == _structure.back()
  for (uint16_t ind = 1; ind < _structure.size(); ind++) {
    _structure[ind - 1]->_nthOutputTaken = _base;
    // std::cout << "_base " << _base << std::endl;
    _structure[ind]->_subBuckets.push_back(_structure[ind - 1]);
  }

  if (_setting->verbosity > 3)
    std::cout << "All Buckets in the last bucket!" << std::endl;
}

void Cascade::CreateTotalizerEncodeTree() {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  TimeMeasurement timeEncodeTree(&_dgpw->_timeVariables->createTree, true);

  _structure.back()->_isLastBucket = true;
  _structure.back()->CreateTotalizerEncodeTree();

  if (_setting->createGraphFile != "")
    _structure.back()->_sorter->_outputTree->DumpOutputTree(
        _setting->createGraphFile + std::to_string(_structure.back()->size()) +
            ".tgf",
        false);

  if (_setting->verbosity > 0)
    std::cout << "c #max sorter depth......: "
              << _structure.back()->_sorter->_outputTree->_depth + 1
              << std::endl;
  //    std::cout << "SIZE OF TREE: " <<
  //    _structure.back()->_sorter->_outputTree->_size << std::endl;
  if (_setting->verbosity < 1) return;

  if (_setting->verbosity > 3)
    std::cout << std::endl << "Totalizer Tree encoded!" << std::endl;
}

void Cascade::DumpNumberOfBucketsAndClauses() {
  std::cout << std::endl;
  // Dump BucketID
  DumpNumberOfBucketEntriesOrClauses(false, false, false, false, false, false);
  std::cout << std::endl;

  // Dump TopBucketEntries
  DumpNumberOfBucketEntriesOrClauses(true, false, false, false, false, false);
  // Dump BottomBucketEntries
  DumpNumberOfBucketEntriesOrClauses(false, true, false, false, false, false);
  std::cout << std::endl;

  // Dump Estimated TopBinaryClauses
  DumpNumberOfBucketEntriesOrClauses(true, false, true, false, true, false);
  // Dump Calculated TopBinaryClauses
  DumpNumberOfBucketEntriesOrClauses(true, false, false, true, true, false);
  std::cout << std::endl;

  // Dump Estimated TopTernaryClauses
  DumpNumberOfBucketEntriesOrClauses(true, false, true, false, false, true);
  // Dump Calculated TopTernaryClauses
  DumpNumberOfBucketEntriesOrClauses(true, false, false, true, false, true);
  std::cout << std::endl;

  // Dump Estimated Bottom Binary Clauses
  DumpNumberOfBucketEntriesOrClauses(false, true, true, false, true, false);
  // Dump Calculated Bottom Binary Clauses
  DumpNumberOfBucketEntriesOrClauses(false, true, false, true, true, false);
  std::cout << std::endl;

  // Dump Estimated Bottom Ternary Clauses
  DumpNumberOfBucketEntriesOrClauses(false, true, true, false, false, true);
  // Dump Calculated Bottom Ternary Clauses
  DumpNumberOfBucketEntriesOrClauses(false, true, false, true, false, true);
  std::cout << std::endl;
}

void Cascade::DumpNumberOfBucketEntriesOrClauses(bool top, bool bottom,
                                                 bool estimated,
                                                 bool calculated, bool binary,
                                                 bool ternary) {
  if (!top && !bottom && !estimated && !calculated)
    std::cout << std::setw(35) << "Bucket ID: ( ";
  else if (top && !estimated && !calculated)
    std::cout << std::setw(35) << "TopBucket Entries: ( ";
  else if (bottom && !estimated && !calculated)
    std::cout << std::setw(35) << "BottomBucket Entries: ( ";
  else if (estimated && top && binary)
    std::cout << std::setw(35) << "Estimated BinaryTopClauses: ( ";
  else if (estimated && top && ternary)
    std::cout << std::setw(35) << "Estimated TernaryTopClauses: ( ";
  else if (calculated && top && binary)
    std::cout << std::setw(35) << "Calculated BinaryTopClauses: ( ";
  else if (calculated && top && ternary)
    std::cout << std::setw(35) << "Calculated TernaryTopClauses: ( ";
  else if (estimated && bottom && binary)
    std::cout << std::setw(35) << "Estimated BinaryBottomClauses: ( ";
  else if (estimated && bottom && ternary)
    std::cout << std::setw(35) << "Estimated TernaryBottomClauses: ( ";
  else if (calculated && bottom && binary)
    std::cout << std::setw(35) << "Calculated BinaryBottomClauses: ( ";
  else if (calculated && bottom && ternary)
    std::cout << std::setw(35) << "Calculated TernaryBottomClauses: ( ";

  uint32_t actualValue(0);
  uint32_t sum(0);

  // for (auto bucket : stuint16_t ind = 0; ind < _structure.size(); ind++)
  for (uint16_t ind = 0; ind < _structure.size(); ind++) {
    if (!top && !bottom && !estimated && !calculated)
      actualValue = ind;
    else if (top && !estimated && !calculated)
      actualValue = _structure[ind]->_topEntries;
    else if (bottom && !estimated && !calculated)
      actualValue = static_cast<uint32_t>(_structure[ind]->_outputs->size());
    else if (estimated && top && binary)
      actualValue = _structure[ind]->_binaryTopClEstimated;
    else if (estimated && top && ternary)
      actualValue = _structure[ind]->_ternaryTopClEstimated;
    else if (calculated && top && binary)
      actualValue = _structure[ind]->_binaryTopCl;
    else if (calculated && top && ternary)
      actualValue = _structure[ind]->_ternaryTopCl;
    else if (estimated && bottom && binary)
      actualValue = _structure[ind]->_binaryBottomClEstimated;
    else if (estimated && bottom && ternary)
      actualValue = _structure[ind]->_ternaryBottomClEstimated;
    else if (calculated && bottom && binary)
      actualValue = _structure[ind]->_binaryBottomCl;
    else if (calculated && bottom && ternary)
      actualValue = _structure[ind]->_ternaryBottomCl;

    std::cout << std::setw(5) << actualValue;
    if (ind < _structure.size() - 1) std::cout << ", ";
    sum += actualValue;
  }
  if (!top && !bottom && !estimated && !calculated)
    std::cout << " )" << std::endl;
  else
    std::cout << " )   sum =" << std::setw(6) << sum << std::endl;
}

uint32_t Cascade::SolveAllTares() {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  TimeMeasurement timeSolvingTares(&_dgpw->_timeVariables->solvingTares, true);
  uint32_t currentresult(UNKNOWN);
  std::vector<uint32_t> collectedAssumptions;

  if (_setting->verbosity > 0) {
    std::cout << std::endl
              << std::setw(90) << "Minimize all Tares from now on!"
              << std::endl;
    if (_setting->verbosity > 1) {
      std::cout << "----------------temporary solve "
                   "structure--------------------------------------------------"
                << std::endl;
    }
  }

  uint16_t sizeMinus = 2;
  if (_onlyByTares) sizeMinus = 1;

  // start with second last bucket!
  for (int32_t ind = static_cast<int32_t>(_structure.size() - sizeMinus);
       ind >= 0; ind--) {
    // assert(_estimatedWeightBoundaries[1] >= _satWeight);
    collectedAssumptions.push_back(_structure[ind]->_tares[0] << 1);
    currentresult = _dgpw->Solve(collectedAssumptions);

    //        if (currentresult == ANTOM_SAT && _multipleCascade != nullptrptr)
    //            _multipleCascade->CalculateWeightBoundaries(_structure[ind]->_multiplicator);
    //        else if (currentresult == ANTOM_UNSAT && _multipleCascade !=
    //        nullptrptr)
    //            _multipleCascade->CalculateWeightBoundaries(-_structure[ind]->_multiplicator);

    if (currentresult == SATISFIABLE) {
      if (_setting->verbosity > 3) std::cout << "SAT" << std::endl;
//            _dgpw->_lastModel = _dgpw->Model();
#ifndef NDEBUG
      bool rst = _dgpw->AddUnit(collectedAssumptions.back());
      assert(rst);
#else
      _dgpw->AddUnit(collectedAssumptions.back());
#endif
      assert(collectedAssumptions.size() == 1);
      _dgpw->CalculateOverallOptimum(_satWeight, true);

    } else if (currentresult == UNSAT) {
      if (_setting->verbosity > 3) {
        std::cout << "UNSAT" << std::endl;
      }
#ifndef NDEBUG
      bool rst = _dgpw->AddUnit(collectedAssumptions.back() ^ 1);
      assert(rst);
#else
      _dgpw->AddUnit(collectedAssumptions.back() ^ 1);
#endif
    }
    if (_setting->verbosity > 2) {
      std::cout << std::setw(50) << "Current Bucket Multiplicator: "
                << _structure[ind]->_multiplicator << std::endl;
      if (currentresult == SATISFIABLE) {
        std::cout << std::setw(34) << "All Tares of Bucket " << ind
                  << " could be set! " << std::setw(40) << "ANTOM_SAT"
                  << std::endl;
        DumpModelOfTares(5);
      } else if (currentresult == UNSAT)
        std::cout << std::setw(31) << "At least one Tare of Bucket " << ind
                  << " couldn't be set! " << std::setw(40) << "ANTOM_UNSAT"
                  << std::endl;
      else
        std::cout << std::setw(88) << "TIMEOUT: " << currentresult << std::endl;
      currentresult = SATISFIABLE;
    }

    if (_setting->verbosity > 3)
      std::cout << "-----------------------------------------------------------"
                   "--------------"
                   "------------------"
                << std::endl;
    collectedAssumptions.pop_back();
  }

  return currentresult;
}

uint32_t Cascade::SolveTares(bool onlyWithAssumptions,
                             bool solveTareInLastBucketToo) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;
  //    std::cout << std::setw(50) << "Weight boundaries: " <<
  //    _estimatedWeightBoundaries[0] << " / " << _estimatedWeightBoundaries[1]
  //    << std::endl;
  if (_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] == 0 ||
      _structure.size() == 1) {
    return SATISFIABLE;
  }

  TimeMeasurement timeSolvingTares(&_dgpw->_timeVariables->solvingTares, true);
  uint32_t currentresult(UNKNOWN);

  if (_setting->weightPlusOne)
    return SolveTareWeightPlusOne(onlyWithAssumptions);

  if (onlyWithAssumptions) {
    std::cout
        << "c onlyWithAssumptions not yet implemented, to use that use the "
           "solveTareWeightPlusOne option!!!";
  }
  assert(!onlyWithAssumptions);

  if (_setting->verbosity > 0) {
    if (_setting->verbosity > 4) std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << std::endl
              << std::setw(90) << "Minimize Tares from now on!" << std::endl;
    if (_setting->verbosity > 1) {
      std::cout << "-----------------------------------------------------------"
                   "--------------"
                   "------------------"
                << std::endl;
    }
  }

  uint16_t sizeMinus = 2;
  if (solveTareInLastBucketToo || _onlyByTares) sizeMinus = 1;

  // start with second last bucket!
  for (int16_t ind = static_cast<int16_t>(_structure.size() - sizeMinus);
       ind >= 0; ind--) {
    uint64_t diffEstimatedCurWeight;

    diffEstimatedCurWeight = _estimatedWeightBoundaries[1] - _dgpw->_satWeight;

    currentresult = _structure[ind]->SolveTares(diffEstimatedCurWeight);

    if (_setting->verbosity > 1) {
      std::cout << std::setw(50)
                << "Weight boundaries: " << _estimatedWeightBoundaries[0]
                << " / " << _estimatedWeightBoundaries[1] << std::endl;
      if (_setting->verbosity > 2) {
        std::cout << std::setw(50)
                  << "Current SAT Weight dgpw: " << _dgpw->_satWeight
                  << std::endl;
        std::cout << std::setw(50) << "Current Bucket Multiplicator: "
                  << _structure[ind]->_multiplicator << std::endl;
        if (currentresult == SATISFIABLE) {
          std::cout << std::setw(34) << "All Tares of Bucket " << ind
                    << " could be set! " << std::setw(40) << "ANTOM_SAT"
                    << std::endl;
          DumpModelOfTares(5);
        } else if (currentresult == UNSAT)
          std::cout << std::setw(31) << "At least one Tare of Bucket " << ind
                    << " couldn't be set! " << std::setw(40) << "ANTOM_UNSAT"
                    << std::endl;
        else
          std::cout << std::setw(88) << "TIMEOUT: " << currentresult
                    << std::endl;
      }

      if (ind != 0)
        std::cout << "---------------------------------------------------------"
                     "------------"
                     "----------------------"
                  << std::endl;
      else
        std::cout << "========================================================="
                     "============"
                     "======================"
                  << std::endl;
    }

    // CASE TIMEOUT
    if (currentresult == UNKNOWN /*|| _control->ReachedLimits()*/) {
      return UNKNOWN;
    }
    //        std::cout << _estimatedWeightBoundaries[0] << std::endl;
    //        std::cout << _dgpw->_satWeight << std::endl;
    assert(_estimatedWeightBoundaries[0] <=
           static_cast<int64_t>(_dgpw->_satWeight));
    // assert(static_cast<int64_t>(_dgpw->_satWeight) <=
    // _estimatedWeightBoundaries[1]);
  }
  if (_setting->verbosity < 3) return SATISFIABLE;
  //    DumpModelOfTares(3);
  if (_setting->verbosity > 0)
    std::cout << "All Tares are solved!" << std::endl << std::endl;
  return SATISFIABLE;
}

void Cascade::CalculateBucketEntries() {
  uint32_t sumOfSizes(0);
  uint32_t maxBucketEntries(0);

  for (auto bucket : _structure) {
    sumOfSizes += bucket->size();
    if (bucket->size() > maxBucketEntries) maxBucketEntries = bucket->size();
  }
  if (_setting->verbosity < 1) return;
  std::cout << "c #buckets...............: " << _structure.size() << std::endl;
  std::cout << "c max Bucket entries.....: " << GetMaxBucketSize() << std::endl;
  std::cout << "c average Bucket entries.: " << sumOfSizes / _structure.size()
            << std::endl;
}

uint64_t Cascade::CountSatisfiedSoftClauses(
    std::vector<SoftClause *> softclauses, const std::vector<uint32_t> &model,
    bool addWeight) {
  uint64_t result(0);
  // Proceed all soft clauses
  for (uint32_t i = 0; i != softclauses.size(); ++i) {
    uint32_t relaxlit = softclauses[i]->relaxationLit;

    // Just proceed satisfied triggers
    if (_dgpw->Model(relaxlit >> 1) == relaxlit) {
      std::vector<uint32_t> clause(softclauses[i]->clause);
      uint32_t pos = 0;
      for (; pos != clause.size(); ++pos) {
        // clause satisfied without trigger?
        if (_dgpw->Model(clause[pos] >> 1) == clause[pos]) {
          if (addWeight)
            result += softclauses[i]->weight;
          else
            result += 1;
          break;
        }
      }
    } else if (_dgpw->Model(relaxlit >> 1) != 0) {
      assert(_dgpw->Model(relaxlit >> 1) == (relaxlit ^ 1));
      if (addWeight)
        result += softclauses[i]->weight;
      else
        result += 1;
    }
  }
  return result;
}

uint64_t Cascade::CountSatisfiedSoftClauses(
    Bucket *bucket, const std::vector<uint32_t> &model) {
  bool addWeight;
  std::vector<SoftClause *> softclauses;

  if (bucket == nullptr || bucket->_isLastBucket) {
    addWeight = true;
    softclauses = _softClauses;
  } else {
    addWeight = false;
    softclauses = *bucket->_softClauses;
    //        std::cout << "How many Softclauses: " << softclauses.size() <<
    //        std::endl;
  }

  uint64_t result = CountSatisfiedSoftClauses(softclauses, model, addWeight);

  //    std::cout << "Fulfilled Softclauses: " << result << std::endl;
  if (addWeight) _satWeight = result > _satWeight ? result : _satWeight;

  return result;
}

void Cascade::CalculateTotalBucketEntries(
    std::vector<SoftClauseNodes *> *SCTree, bool add) {
  _totalBucketEntriesperWeight = 0;
  _totalBucketOccurrences = 0;
  if (!add) {
    for (uint32_t i = 0; i <= _numberOfBuckets; ++i) {
      if (i >= _totalBucketEntries.size())
        _totalBucketEntries.push_back(0);
      else
        _totalBucketEntries[i] = 0;
    }
  }
  for (uint32_t i = 0; i < SCTree->size(); ++i) {
    _totalBucketEntriesperWeight += (*SCTree)[i]->inHowManyBuckets;
    _totalBucketOccurrences +=
        (*SCTree)[i]->GetOccurrences() * (*SCTree)[i]->size();
    for (uint32_t j = 0; j <= (*SCTree)[i]->highestBucket; ++j) {
      if ((*SCTree)[i]->hasSubBuckets &&
          (*SCTree)[i]->occursHowOftenInBucket[j])
        _totalBucketEntries[j] +=
            (*SCTree)[i]->occursHowOftenInBucket[j] * (*SCTree)[i]->size();
      else
        _totalBucketEntries[j] += (*SCTree)[i]->occursHowOftenInBucket[j];

      //            if ((*SCTree)[i]->occursHowOftenInBucket[j] >= 1)
      //                _totalBucketEntries[j]++;
      //            std::cout << "_totalBucketEntries[" << j << "]: " <<
      //            _totalBucketEntries[j] << std::endl;
    }
  }
  //    std::cout << "_totalBucketEntriesperWeight: " <<
  //    _totalBucketEntriesperWeight << std::endl; std::cout <<
  //    "_totalBucketOccurrences: " << _totalBucketOccurrences << std::endl;
}

bool Cascade::AtLeastTwoBucketsInCommon(SoftClauseNodes *SCN1,
                                        SoftClauseNodes *SCN2) {
  uint32_t minSize = (SCN1->occursHowOftenInBucket.size() <
                      SCN2->occursHowOftenInBucket.size())
                         ? SCN1->occursHowOftenInBucket.size()
                         : SCN2->occursHowOftenInBucket.size();
  uint32_t overlappings(0);
  for (uint32_t ind = 0; ind < minSize; ind++) {
    if (SCN1->occursHowOftenInBucket[ind] > 0 &&
        SCN2->occursHowOftenInBucket[ind] > 0) {
      overlappings++;
      if (overlappings > 1) return true;
    }
  }
  return false;
}

void Cascade::DumpSCNodeStructure(std::vector<SoftClauseNodes *> *dumpingSCTree,
                                  uint16_t verbosity) {
  if (_setting->verbosity < verbosity || dumpingSCTree->size() > 1000) return;
  std::cout << std::endl;
  std::cout << std::setw(4) << "" << std::setw(8) << "#" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "" << std::setw(8) << "O" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "" << std::setw(8) << "c" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "" << std::setw(8) << "c" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "#" << std::setw(8) << "u" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "B" << std::setw(8) << "r" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "u" << std::setw(8) << "r" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "c" << std::setw(8) << "e" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "k" << std::setw(8) << "n" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "e" << std::setw(8) << "c" << std::setw(15) << ""
            << std::setw(3) << "|"
            << "   " << _base << " up to the power of:" << std::endl;
  std::cout << std::setw(4) << "t" << std::setw(8) << "e" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "s" << std::setw(8) << "s" << std::setw(15)
            << "Weight" << std::setw(3) << "|";
  for (uint32_t i = 0; i <= _numberOfBuckets; ++i) {
    std::cout << std::setw(4) << i;
  }
  std::cout << std::endl;
  std::cout << "--#occurrences---------------|";

  for (uint32_t i = 0; i <= _numberOfBuckets; ++i) {
    std::cout << "----";
  }
  std::cout << "----" << std::endl;

  for (uint32_t i = 0; i != (*dumpingSCTree).size(); ++i) {
    (*dumpingSCTree)[i]->dumpStructure(true, i);
  }
  std::cout << "-----------------------------|";
  for (uint32_t i = 0; i <= _numberOfBuckets; ++i) {
    std::cout << "----";
  }
  std::cout << "----" << std::endl;

  std::cout << std::setw(4) << _totalBucketEntriesperWeight << std::setw(8)
            << _totalBucketOccurrences << std::setw(15) << _sumOfSoftWeights
            << std::setw(3) << "|";
  for (uint32_t i = 0; i < _totalBucketEntries.size(); ++i) {
    std::cout << std::setw(4) << _totalBucketEntries[i];
  }
  std::cout << std::endl;
}

void Cascade::DumpBucketStructure(bool onlyLastBucket, uint16_t verbosity) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_setting->verbosity < verbosity) return;
  uint16_t depth(0);
  uint16_t maxDepth(0);
  std::cout << "(SoftClauseEntries, SorterEntries, Multiplicator)" << std::endl;
  if (onlyLastBucket) {
    std::cout << "Last Bucket: " << std::endl;
    depth = _structure.back()->DumpAndGetMaxDepth(0);
    std::cout << "LastBucket has a depth of: " << depth << std::endl;
    std::cout << std::endl;
    return;
  }
  for (uint32_t i = 0; i < _structure.size(); i++) {
    std::cout << "Bucket " << i << ":" << std::endl;
    depth = _structure[i]->DumpAndGetMaxDepth(0);
    if (depth > maxDepth) maxDepth = depth;
    std::cout << "Bucket " << i << " has a depth of: " << depth << std::endl;
    std::cout << std::endl;
  }
  std::cout << "Max depth of all buckets is: " << maxDepth << std::endl
            << std::endl;
}

uint64_t Cascade::GetHighestMultiplicator() {
  assert(_highestBucketMultiplicator != 0);
  return _highestBucketMultiplicator;
}

uint32_t Cascade::GetMaxBucketSize() {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  uint32_t maxSize = 0;
  for (auto bucket : _structure)
    maxSize = bucket->size(true) * bucket->_nthOutputTaken > maxSize
                  ? bucket->size(true) * bucket->_nthOutputTaken
                  : maxSize;
  //    for (uint32_t ind=0; ind < _structure.size(); ind++)
  //    {
  //        std::cout << "Bucket " << ind << " has a size of " <<
  //        _structure[ind]->size(true) * _structure[ind]->_nthOutputTaken << "
  //        and every " << _structure[ind]->_nthOutputTaken << " output is
  //        taken." << std::endl; maxSize = ((_structure[ind]->size(true) *
  //        _structure[ind]->_nthOutputTaken) > maxSize) ?
  //        (_structure[ind]->size(true) * _structure[ind]->_nthOutputTaken) :
  //        maxSize;
  //    }
  return maxSize;
}

bool Cascade::AddNewBucketsTillSizeBoundary(uint32_t maxSize,
                                            bool onlySolveWithTares,
                                            bool addTareToLastBucket) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_setting->verbosity > 3)
    std::cout << "Add another bucket because of size boundary: " << maxSize
              << std::endl;

  while (_structure.back()->size(true) > maxSize) {
    if (!AddNewBucketsTillMultiplicatorMatches(_highestBucketMultiplicator * 2,
                                               onlySolveWithTares,
                                               addTareToLastBucket))
      return false;
  }
  return true;
}

bool Cascade::AddNewBucketsTillMultiplicatorMatches(uint64_t maxMultiplicator,
                                                    bool onlySolveWithTares,
                                                    bool addTareToLastBucket) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_setting->verbosity > 3)
    std::cout << "_highestBucketMultiplicator: " << _highestBucketMultiplicator
              << std::endl;
  assert(_structure.back()->_multiplicator <= maxMultiplicator);

  while (_structure.back()->_multiplicator < maxMultiplicator) {
    if (_setting->verbosity > 3) {
      std::cout << "_structure.back()->size(true): "
                << _structure.back()->size(true) << std::endl;
      std::cout << "onlySolveWithTares: " << onlySolveWithTares << std::endl;
    }

    // This bucket is too small to connect with another bucket.
    // It contains only tare and one result from the bucket before!
    // Calculate the results only by the tares!
    // TOBI: VERIFY!!! - Should be possible to calculate all results only by
    // tares
    if (_structure.back()->size(true) == 1 && onlySolveWithTares) {
      //            if (_setting->verbosity > 3)
      //                std::cout << "c MORE THAN TWO ENTRIES..: TRUE" <<
      //                std::endl;
      if (_setting->verbosity > 3)
        std::cout << "Add one tare to last bucket with size one!" << std::endl;
      AddTare(_structure.size() - 1);

      // here is the last position the bucket can be dumped
      // otherwise the subbuckets are encoded and not dumpable anymore :-)
      DumpBucketStructure(true, 4);

      CreateTotalizerEncodeTree();
      CalculateBucketEntries();

      // at least the tare should be 0, then by solving max tares to 1 is asked
      // for!
      uint32_t unitClauseVar =
          (_structure.back()->_sorter->GetOrEncodeOutput(0) << 1) ^ 1;

      if (_setting->verbosity > 3)
        std::cout << "Unit Clause for first entry in last Bucket Added: "
                  << unitClauseVar << std::endl;
      _dgpw->AddUnit(unitClauseVar);
      return false;
    } else {
      // if (_structure.back()->size(true) >= 1)
      AddAdditionalBucket();
    }

    // at least one element + one tare!!
    assert(_structure.back()->size(true) > 0);
  }

  if (addTareToLastBucket && _structure.back()->_tares.empty()) {
    if (_setting->verbosity > 3)
      std::cout << std::endl
                << "There is no tare at the last bucket! Add one!" << std::endl;
    AddTare(_structure.size() - 1);
  }

  return true;
}

void Cascade::AddTare(unsigned long position) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  uint32_t tare(_dgpw->NewVariable());
  _structure[position]->AddTare(tare);
  _tareWeight += _structure[position]->_multiplicator;
  if (_setting->verbosity < 3) return;

  std::cout << "   Tare added: " << tare << std::endl;
}

void Cascade::AddAdditionalBucket() {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_setting->verbosity > 3) {
    std::cout << "_highestBucketMultiplicator: " << _highestBucketMultiplicator
              << std::endl;
    std::cout << "Add another Bucket with multiplicator: "
              << _highestBucketMultiplicator;
  }
  _structure.back()->_isLastBucket = false;
  if (_structure.back()->_tares.empty()) {
    if (_setting->verbosity > 3)
      std::cout << std::endl
                << "There is no tare at the last bucket! Add one!" << std::endl;
    AddTare(_structure.size() - 1);
  }

  _structure.push_back(new Bucket(_dgpw, this));
  _structure.back()->_isLastBucket = true;
  _highestBucketMultiplicator =
      static_cast<uint64_t>(pow(_base, _structure.size() - 1));
  _structure.back()->_multiplicator = _highestBucketMultiplicator;
  _structure[_structure.size() - 2]->_nthOutputTaken = _base;
  _structure.back()->_subBuckets.push_back(_structure[_structure.size() - 2]);
  if (_structure.back()->_subBuckets.back()->_localMaxPos !=
      static_cast<uint32_t>(-1))
    _structure.back()->_localMaxPos = static_cast<uint32_t>(
        floor(_structure.back()->_subBuckets.back()->_localMaxPos / 2));
}

std::vector<uint32_t> Cascade::CalculateAssumptionsFor(int64_t weight,
                                                       int32_t startingPos) {
  if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (startingPos == -1) {
    return {};
  }

  assert(weight <= _estimatedWeightBoundaries[1]);
  assert(weight >= _estimatedWeightBoundaries[0]);
  assert(weight >= static_cast<int64_t>(_dgpw->_satWeight));

  //    std::cout << "c Calculate assumptions for weight: " << weight <<
  //    std::endl;

  std::vector<uint32_t> collectedAssumptions;

  int64_t upperDiff = _estimatedWeightBoundaries[1] - weight;
  int64_t lowerDiff = weight - _estimatedWeightBoundaries[0];

  // start with second last bucket!
  for (int32_t ind = startingPos; ind >= 0; --ind) {
    int64_t actualMult = static_cast<int64_t>(_structure[ind]->_multiplicator);

    if (_setting->verbosity > 5) {
      std::cout << std::setw(50) << "c Actual Diffs: " << lowerDiff << " / "
                << upperDiff << std::endl;
      std::cout << std::setw(50) << "c actualMult: " << actualMult << std::endl;
    }

    assert(upperDiff + lowerDiff == 2 * actualMult - 1);

    if (lowerDiff >= actualMult) {
      assert(upperDiff < actualMult);
      collectedAssumptions.push_back(_structure[ind]->_tares[0] << 1);
      lowerDiff -= actualMult;
    } else if (upperDiff >= actualMult) {
      assert(lowerDiff < actualMult);
      collectedAssumptions.push_back(_structure[ind]->_tares[0] << 1 ^ 1);
      upperDiff -= actualMult;
    } else {
      assert(false);
    }
  }

  for (auto assumption : _fixedTareAssumption) {
    collectedAssumptions.push_back(assumption);
  }

  if (_setting->verbosity < 3) return collectedAssumptions;

  std::cout << std::setw(51) << "Assumptions: (";
  for (auto assumption : collectedAssumptions) {
    std::cout << assumption << ", ";
  }
  std::cout << ")" << std::endl;
  return collectedAssumptions;
}

int32_t Cascade::SetUnitClauses(int32_t startingPos) {
  if (_setting->verbosity > 6) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
  }

  // start with second last bucket!
  for (int32_t ind = startingPos; ind >= 0; ind--) {
    int64_t actualMult = static_cast<int64_t>(_structure[ind]->_multiplicator);

    if (_setting->verbosity > 5) {
      std::cout << std::setw(50)
                << "Weight boundaries: " << _estimatedWeightBoundaries[0]
                << " / " << _estimatedWeightBoundaries[1] << std::endl;
      std::cout << std::setw(50) << "actualMult: " << actualMult << std::endl;
    }

    assert(_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] ==
           2 * actualMult - 1);
    // The actual result implies that the currentTare and maybe more tares can
    // be set directly to TRUE
    if (_estimatedWeightBoundaries[1] -
            static_cast<int64_t>(_dgpw->_satWeight) <
        actualMult) {
//            std::cout << _structure[ind]->_tares[0]<< std::endl;
#ifndef NDEBUG
      bool rst = _dgpw->AddUnit(_structure[ind]->_tares[0] << 1);
      assert(rst);
#else
      _dgpw->AddUnit(_structure[ind]->_tares[0] << 1);
#endif
      _estimatedWeightBoundaries[0] += actualMult;
      continue;
    }
    // Corner Case!
    // If new result is already larger as maximal possible weight, we do not
    // need to try,
    // -> the corresponding tare can be directly set to FALSE
    else if (_estimatedWeightBoundaries[1] - actualMult >=
             static_cast<int64_t>(_dgpw->_sumOfSoftWeights)) {
#ifndef NDEBUG
      bool rst = _dgpw->AddUnit((_structure[ind]->_tares[0] << 1) ^ 1);
      assert(rst);
#else
      _dgpw->AddUnit((_structure[ind]->_tares[0] << 1) ^ 1);
#endif

      _estimatedWeightBoundaries[1] -= actualMult;
      continue;
    }
    // set assumption vector.
    else {
      if (_setting->verbosity > 2) {
        std::cout << std::setw(50)
                  << "Weight boundaries: " << _estimatedWeightBoundaries[0]
                  << " / " << _estimatedWeightBoundaries[1] << std::endl;
        std::cout << std::setw(50) << "actualMult: " << actualMult << std::endl;
      }
      return ind;
    }
  }

  if (_estimatedWeightBoundaries[1] - _dgpw->_satWeight > 0)
    return 0;
  else
    // result is found.
    return -1;
}

uint32_t Cascade::SolveTareWeightPlusOne(bool onlyWithAssumptions) {
  TimeMeasurement timeSolvingTares(&_dgpw->_timeVariables->solvingTares, true);
  uint32_t currentresult(SATISFIABLE);
  std::vector<uint32_t> collectedAssumptions;
  std::vector<uint32_t> collectedAssumptionsMinusOne;
  std::vector<uint32_t> lastCollectedAssumptions;

  if (_setting->verbosity > 0) {
    if (_setting->verbosity > 6) std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << std::endl
              << std::setw(90) << "Try to set tares for actual weight + 1!"
              << std::endl;
    if (_setting->verbosity > 1)
      std::cout << "-----------------------------------------------------------"
                   "--------------"
                   "------------------"
                << std::endl;
  }

  int32_t startingPos = _structure.size() - 2;

  if (_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] >
      static_cast<int64_t>(_highestBucketMultiplicator))
    startingPos++;
  //    else if (_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0]
  //    == 1 && _highestBucketMultiplicator <= 2)
  //        //border case, only 2 buckets... est weight diffs = 1
  //        startingPos++;
  else
    assert(_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] >=
               static_cast<int64_t>(_highestBucketMultiplicator / 2) ||
           _estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] == 1);

  //    std::cout << "startingPos: " << startingPos << std::endl;

  while (currentresult == SATISFIABLE) {
    if (!onlyWithAssumptions) {
      startingPos = SetUnitClauses(startingPos);
      //            std::cout << "SP: " << startingPos << std::endl;
    }

    if (startingPos == -1 || static_cast<int64_t>(_dgpw->_satWeight) ==
                                 _estimatedWeightBoundaries[1]) {
      collectedAssumptions = CalculateAssumptionsFor(
          static_cast<int64_t>(_dgpw->_satWeight), startingPos);
      break;
    }

    //        collectedAssumptionsMinusOne =
    //        CalculateAssumptionsFor(static_cast<int64_t>(_dgpw->_satWeight),
    //        startingPos); lastCollectedAssumptions = collectedAssumptions;
    collectedAssumptions = CalculateAssumptionsFor(
        static_cast<int64_t>(_dgpw->_satWeight) + 1, startingPos);

    currentresult = _dgpw->Solve(collectedAssumptions);
    //        std::cout << "tried SATWeight: " << _dgpw->_satWeight + 1 <<
    //        std::endl;
    if (_setting->verbosity > 1)
      std::cout << "-----------------------------------------------------------"
                   "--------------"
                   "------------------"
                << std::endl;

    //        std::cout << "CURRENT RESULT: " << currentresult << std::endl;
    if (currentresult == SATISFIABLE) {
      _dgpw->CalculateOverallOptimum(0, true);
      //            std::cout << "real SATWeight:  " << _dgpw->_satWeight <<
      //            std::endl;
    }
  }
  //    collectedAssumptions =
  //    CalculateAssumptionsFor(static_cast<int64_t>(_dgpw->_satWeight) - 1,
  //    startingPos);

  if (currentresult == SATISFIABLE) {
    std::cout << "c SAT AFTER SOLVING TARES!" << std::endl;
    //        for (auto unitClause : collectedAssumptions) {
    //            _dgpw->AddUnit(unitClause);
    //        }
  } else {
    //        std::cout << "collect assumptions for: " << _dgpw->_satWeight <<
    //        std::endl;
    if (_setting->verbosity > 0)
      std::cout
          << "c UNSAT AFTER SOLVING TARES!, solve again with right assumptions!"
          << _dgpw->_satWeight << std::endl;
    collectedAssumptions = CalculateAssumptionsFor(
        static_cast<int64_t>(_dgpw->_satWeight), startingPos);

    //        assert(_dgpw->Solve(collectedAssumptions)==SATISFIABLE);
    //        if (onlyWithAssumptions && _dgpw->Solve(collectedAssumptions) !=
    //        10)
    //        {
    //            assert(false);
    //        }
    //        std::cout << "calc!1" << std::endl;
    //        _dgpw->CalculateOverallOptimum(0,true);
    //        collectedAssumptions =
    //        CalculateAssumptionsFor(static_cast<int64_t>(_dgpw->_satWeight) +
    //        1, startingPos);
    //        assert(_dgpw->Solve(collectedAssumptions)!=SATISFIABLE);
  }
  for (auto unitClause : collectedAssumptions) {
    //        _fixedTareAssumption.clear();
    if (onlyWithAssumptions) {
      _fixedTareAssumption.push_back(unitClause);
      //            std::cout << "TareAssumptions: " << unitClause << std::endl;
    } else {
      _dgpw->AddUnit(unitClause);
      //            std::cout << "AddUnit: " << unitClause << std::endl;
    }
  }
  if (onlyWithAssumptions && _dgpw->Solve(_fixedTareAssumption) != 10) {
    std::cout << "c Wrong result!" << std::endl;
    assert(false);
  } else if (!onlyWithAssumptions && _dgpw->Solve() != 10) {
    std::cout << "c Wrong Solve Result!" << std::endl;
    assert(false);
  }
  _dgpw->CalculateOverallOptimum(0, true);
  //    std::cout << "SolveValue: " << _dgpw->Solve() << std::endl;
  //    std::cout << "SATWeight: " << _dgpw->_satWeight << std::endl;
  //    _dgpw->CalculateOverallOptimum(0,true);
  //    std::cout << "SATWeight: " << _dgpw->_satWeight << std::endl;
  //    collectedAssumptions =
  //    CalculateAssumptionsFor(static_cast<int64_t>(_dgpw->_satWeight) + 1,
  //    startingPos); currentresult = _dgpw->Solve(collectedAssumptions);
  //     assert(_dgpw->Solve() == SATISFIABLE);

  //    std::cout << "SOLVE TARE WEIGHT + 1 - currentresult: " << currentresult
  //    << std::endl; assert(currentresult == SATISFIABLE);

  if (_setting->verbosity > 1)
    std::cout << "============================================================="
                 "================"
                 "=============="
              << std::endl;

  if (currentresult == UNKNOWN /*|| _control->ReachedLimits()*/) {
    _dgpw->_resultUnknown = true;
    std::cout << std::setw(88) << "TIMEOUT: " << currentresult << std::endl;
    return UNKNOWN;
  }

  if (_setting->verbosity > 1)
    std::cout << "All Tares are solved!" << std::endl << std::endl;
  return SATISFIABLE;
}

}  // namespace DGPW
