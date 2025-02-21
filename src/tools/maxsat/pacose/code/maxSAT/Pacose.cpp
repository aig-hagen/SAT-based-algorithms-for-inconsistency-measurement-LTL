/********************************************************************************************
SATSolverProxy.cpp -- Copyright (c) 2020, Tobias Paxian
    Parts are taken from the QMaxSAT 2017 MaxSAT evaluation version.

Permission is hereby granted, free of charge, to any person obtaining a copy of
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

#include "Pacose.h"
#include "Settings.h"

// to use original antom code
#include "DGPW/dgpw.h"
#include "Encodings.h"
#include "Greedyprepro.h"
#include "Pacose.h"
#include "SATSolverProxy.h"
#include "Softclause.h"
#include "timevariables.h"

#include <math.h>    // std::cout
#include <iostream>  // std::cout
//#include <sstream>      // std::setw
#include <stdlib.h>        // strtoull
#include <sys/resource.h>  // timing
#include <sys/time.h>      // timing
#include <algorithm>       // std::sort
#include <cassert>         // assertions
#include <cstring>         // fast parser
#include <fstream>         // parseWcnfFile std::ifstream
#include <iomanip>         // std::setw
#include <set>             // comparison of weight difference
#include <sstream>         // parseWcnfFile std::istringstream

Pacose::Pacose()
    : _settings()
      //    _satSolver(SATSolverProxy()),
      ,
      _nbVars(0),
      _nbClauses(0),
      _top(0),
      _originalSoftClauses(),
      _hasHardClauses(false),
      _cpuLimit(INT32_MAX),
      _memLimit(INT32_MAX),
      _nbOfOrigVars(0),
      _sumOfSoftWeights(0),
      _overallSoftWeights(0),
      _localUnSatWeight(INT64_MAX),
      _minWeight(0),
      _maxWeight(0),
      _GCD(1),
      _encodings(new Encodings(&_settings)),
      _negRelaxLit(0) {}

Pacose::~Pacose() {}

void Pacose::InitSatSolver(SATSolverType solverType) {
  _satSolver = SATSolverProxy::InitSATSolver(solverType, 1);
  //    _satSolver = SATSolverProxy();
  //            SATSolverProxy::InitSATSolver(solverType, 1);
  if (_settings.reconf != 99) {
    _satSolver->SetReconf(_settings.reconf);
  }
}

static bool isEndOfLine(char *p)

{
  return *p == '\0';
}

static void skipWhiteSpace(char **p)

{
  while (**p == ' ') {
    ++(*p);
  }
}

static long long int fast_atoi(char **p) {
  long long int x = 0;

  bool neg = false;

  if (**p == '-') {
    neg = true;

    ++(*p);
  }

  while (**p >= '0' && **p <= '9') {
    x = (x * 10) + (**p - '0');

    ++(*p);
  }

  return neg ? -x : x;
}

static long long int getInt(char **p) {
  skipWhiteSpace(p);

  return fast_atoi(p);
}

static bool match(char **p, const char *str) {
  for (; *str != '\0'; ++str, ++(*p))

  {
    if (*str != **p) return false;
  }

  return true;
}

bool Pacose::parseWcnfFile(std::string wcnfFile) {
  double parseTimeStart;
  struct rusage resources;
  getrusage(RUSAGE_SELF, &resources);
  parseTimeStart =
      resources.ru_utime.tv_sec + 1.e-6 * (double)resources.ru_utime.tv_usec;

  std::ifstream source;
  source.open(wcnfFile.c_str(), std::ifstream::in);
  assert(source.good());

  std::vector<unsigned> sClause;
  std::string line;

  while (getline(source, line) && source.good()) {
    while (line[0] == ' ' || line[0] == '\t') {
      line.erase(line.begin());
    }

    // comment and invalid lines
    if ((line.length() > 0 && !isdigit(line[0]) &&
         line[0] != 'p')  //  && line[0] != ' ' && line[0] != '\t'
        || line.length() == 0) {
      continue;
    }

    char c_array[line.length() + 1];

    strncpy(c_array, line.c_str(), line.length());

    c_array[line.length()] = '\0';

    char *cLine = c_array;

    // information about the WCNF structure

    if (*cLine == 'p') {
      if (!match(&cLine, "p")) {
        std::cout << "c Header must start with p!" << std::endl;

        exit(1);
      }

      skipWhiteSpace(&cLine);

      if (!match(&cLine, "wcnf")) {
        std::cout << "c File is not a MaxSAT instance!" << std::endl;
        exit(1);
      }

      skipWhiteSpace(&cLine);

      _nbVars = static_cast<unsigned>(getInt(&cLine));

      _nbClauses = static_cast<unsigned>(getInt(&cLine));

      if (isEndOfLine(cLine)) {
        _top = static_cast<uint64_t>(-1);
      } else {
        _top = static_cast<uint64_t>(getInt(&cLine));
      }

      std::cout << "c Number of Variables....: " << _nbVars << std::endl;
      std::cout << "c Number of Clauses......: " << _nbClauses << std::endl;
      std::cout << "c HardClause weight......: " << _top << std::endl;

      _satSolver->NewVariables(_nbVars + 1);

      continue;
    }

    uint64_t weight = static_cast<uint64_t>(getInt(&cLine));
    if (weight == 0) {
      continue;
    }

    skipWhiteSpace(&cLine);

    int literal;

    sClause.clear();
    if (weight >= _top) {
      _hasHardClauses = true;

      _satSolver->NewClause();

      while (!isEndOfLine(cLine)) {
        literal = getInt(&cLine);

        if (literal == 0) {
          break;
        }

        assert(static_cast<unsigned>(abs(literal)) <= _nbVars);

        _satSolver->AddLiteral(SignedToUnsignedLit(literal));
      }

      if (!_satSolver->CommitClause()) {
        // std::cout << "CommitClause went wrong!" << std::endl;

        assert(false);

        return false;
      }

      _satSolver->ResetClause();

    } else {
      while (!isEndOfLine(cLine)) {
        literal = getInt(&cLine);

        if (literal == 0) {
          break;
        }

        assert(static_cast<unsigned>(abs(literal)) <= _nbVars);

        sClause.push_back(SignedToUnsignedLit(literal));
      }

      AddSoftClause(sClause, weight);
    }
  }

  std::cout << "c Number of SoftClauses..: " << _originalSoftClauses.size()
            << std::endl;

  getrusage(RUSAGE_SELF, &resources);

  std::cout << "c parse time.............: "
            << (resources.ru_utime.tv_sec +
                1.e-6 * (double)resources.ru_utime.tv_usec) -
                   parseTimeStart
            << std::endl;

  return true;
}

unsigned Pacose::SignedToUnsignedLit(int literal)

{
  // if literal < 0 then lit=(2*literal)+1 else lit=2*literal

  return (static_cast<unsigned>(abs(literal)) << 1) ^ (literal < 0);
}

void Pacose::AddSoftClause(std::vector<unsigned> &clause, uint64_t weight) {
  unsigned relaxLit = static_cast<unsigned>(_satSolver->NewVariable() << 1);
  //    std::cout << "RL, weight: << " << relaxLit << ", " << weight <<
  //    std::endl;
  SoftClause *SC = new SoftClause(relaxLit, clause, weight);
  _originalSoftClauses.push_back(SC);
  clause.push_back(relaxLit);

  //  for (auto lit : clause) {
  //    std::cout << lit << ", ";
  //  }
  //  std::cout << std::endl;
  _satSolver->AddClause(clause);
  //  if (!_satSolver->AddClause(clause)) {
  //    assert(false);
  //    //        exit(1);
  //  }
}

void Pacose::AddSoftClauseTo(std::vector<SoftClause *> *softClauseVector,
                             std::vector<unsigned> &clause, uint64_t weight) {
  unsigned relaxLit = static_cast<unsigned>(_satSolver->NewVariable() << 1);
  SoftClause *SC = new SoftClause(relaxLit, clause, weight);
  softClauseVector->push_back(SC);
  clause.push_back(relaxLit);

  if (!_satSolver->AddClause(clause)) {
    assert(false);
    //        exit(1);
  }
}

// koshi 2013.04.05, 2013.05.21, 2013.06.28, 2013.07.01, 2013.10.04
// koshi 20140121
void Pacose::HeuristicQMaxSAT(long long int sum, long long int k) {
  printf("c auto-mode for generating cardinality constraints\n");
  int logk = 0;
  int logsum = 0;
  // calculate 2nd logarithm of k
  for (long long int ok = k; ok > 0; ok = ok >> 1) logk++;
  // calculate 2nd logarithm of sum of weights
  for (long long int osum = sum; osum > 0; osum = osum >> 1) logsum++;
  printf("c logk = %d, logsum = %d\n", logk, logsum);
  if (logk + logsum < 15) {
    // Bailleux
    _settings.SetEncoding(BAILLEUX);
    _settings.SetCompression(0);
    printf("c Bailleux's encoding (comp=0)\n");
  } else if (k < 3) {
    // Warners
    _settings.SetEncoding(WARNERS);
    _settings.SetCompression(1);
    printf("c Warners' encoding (comp=1)\n");
  } else if (logsum < 17) {
    // Ogawa
    _settings.SetEncoding(OGAWA);
    _settings.SetCompression(0);
    printf("c Ogawa's encoding (comp=0)\n");
  } else {
    _settings.SetEncoding(WARNERS);
    _settings.SetCompression(1);
    printf("c Warners' encoding (comp=1)\n");
  }
}

void Pacose::wbSortAndFilter(long long int UnSATWeight) {
  unsigned long sizeBefore = _actualSoftClauses->size();

  //  unsigned long weight = 0;
  //  for (auto SC : _originalSoftClauses) weight += SC->originalWeight;

  //  std::cout << "sizeBefore: " << sizeBefore << std::endl;
  //  std::cout << "_originalSoftClauses.size(): " <<
  //  _originalSoftClauses.size()
  //            << "  weight: " << weight << std::endl;

  //  std::cout << "c remaining SCs..........: " << _actualSoftClauses->size()
  //            << std::endl;

  for (unsigned i = 0; i < (*_actualSoftClauses).size(); i++) {
    //    std::cout << "AW: " << (*_actualSoftClauses)[i]->weight << std::endl;
    // @TOBI think through - in my optinion it should be smaller equal!
    if (!((*_actualSoftClauses)[i]->weight <= UnSATWeight)) {
      // not satisfiable anymore - add negated unit clause of blockings
      _satSolver->ResetClause();
      _satSolver->NewClause();
      unsigned int ulit = (*_actualSoftClauses)[i]->relaxationLit ^ 1;
      _satSolver->AddLiteral(&ulit);
      _satSolver->CommitClause();
      //      std::cout << "REMOVE SC!!! " << (*_actualSoftClauses)[i]->weight
      //                << std::endl;
      _actualSoftClauses->erase(_actualSoftClauses->begin() + i);
      i--;
    }
    /*else
    // if weight is still fulfillable push it
    tmpSoftClauses.push_back((*_actualSoftClauses)[i]);
    } */
  }
  //  std::cout << "  _actualSoftClauses.size(): " <<
  //  _actualSoftClauses->size(); std::cout << "_originalSoftClauses.size(): "
  //  << _originalSoftClauses.size(); std::cout << std::endl;
  //    _actualSoftClauses->assign(tmpSoftClauses.begin(),
  //    tmpSoftClauses.end());

  //    tmpSoftClauses.moveTo(_softClauses);

  if (_actualSoftClauses->size() < sizeBefore) {
    //    std::cout << "sizeAfter: " << _actualSoftClauses->size() << std::endl;
    if (_settings.verbosity > 0) {
      std::cout << "c removed SCs by wb......: "
                << sizeBefore - _actualSoftClauses->size() << std::endl;
      std::cout << "c remaining SCs..........: " << _actualSoftClauses->size()
                << std::endl;
    }
    //    unsigned long weight = 0;
    //    for (auto SC : _originalSoftClauses) weight += SC->originalWeight;

    //    std::cout << "_originalSoftClauses.size(): " <<
    //    _originalSoftClauses.size()
    //              << "  weight: " << weight << std::endl;

    if (_satSolver->Solve() != SATISFIABLE) {
      std::cout << "ERROR: Solver call should've been satisfiable!"
                << std::endl;
    }

    CalculateSATWeight();
  }
}

void Pacose::genCardinals(
    long long int tmpUnSATWeight,
    long long int &divisor,  // koshi 2013.10.04
    std::vector<unsigned> &lits, std::vector<unsigned> &linkingVar,
    std::vector<long long int> &linkingWeight,  // uemura 20161202
    std::vector<long long int> &divisors,       // uemura 20161128
    std::vector<std::vector<unsigned>> &linkingVars,
    std::vector<std::vector<long long int>> &linkingWeights,  // uemura 20161128
    int compression) {
  // simple inprocessor, after first solve call
  // filter out weights which cannot be fulfilled anymore!
  // PAX: copys weights and blockings into sweits and sblockings.

  //  std::cout << tmpUnSATWeight << std::endl;
  //  wbSortAndFilter(tmpUnSATWeight);

  //    long long int sum = sumWeight(); // koshi 20140124
  // attention! Do not update original sum of softweights!
  unsigned long sum = 0;
  for (int i = 0; i < _actualSoftClauses->size(); i++) {
    sum += (*_actualSoftClauses)[i]->weight;
  }

  if (_settings.verbosity > 0)
    std::cout << "c Sum of weights = " << sum
              << " unSatWeight: " << tmpUnSATWeight << std::endl;

  //  if (_encoding == HEURISTICQMAXSAT) {
  //    // koshi 20140324 auto mode
  //    HeuristicQMaxSAT(sum, tmpUnSATWeight);
  //  }

  if (_actualSoftClauses->size() == 0) {
    linkingVar.clear();
  } else {
    // original QMaxSAT Encodings need weights and blockings for their recursive
    // calls
    //      std::cout << "NOT DGPW" << std::endl;
    _blockings.clear();
    _weights.clear();
    for (int i = 0; i < _actualSoftClauses->size(); i++) {
      _blockings.push_back((*_actualSoftClauses)[i]->relaxationLit);
      _weights.push_back((*_actualSoftClauses)[i]->weight);
      //      std::cout << _weights.back() << ", " << _blockings.back() <<
      //      std::endl;
    }

    //    std::cout << "ASC: " << _actualSoftClauses->size() << std::endl;
    //    std::cout << "BLO: " << _blockings.size() << std::endl;
    //    std::cout << "WEI: " << _weights.size() << std::endl;
    //    std::cout << "SUM: " << sum << std::endl;
    //    std::cout << "USW: " << tmpUnSATWeight << std::endl;
    //    std::cout << "COM: " << _settings.GetCompression() << std::endl;
    //    std::cout << "LTS: " << lits.size() << std::endl;
    //    std::cout << "LVS: " << linkingVar.size() << std::endl;
    //    std::cout << "LWS: " << linkingWeight.size() << std::endl;
    //    std::cout << "LWSS: " << linkingWeights.size() << std::endl;

    // koshi 20140124 20140129
    // koshi 2013.06.28

    // why using old sum of Softweights and not the newly calculated sum?
    switch (_encoding) {
      case WARNERS:
        _encodings->genWarners0(_weights, _blockings, sum, tmpUnSATWeight,
                                compression, *_satSolver, lits, linkingVar);
        break;
      case BAILLEUX:
        _encodings->genBailleux0(_weights, _blockings, sum, tmpUnSATWeight,
                                 compression, *_satSolver, lits, linkingVar);
        break;
      case ASIN:
        _encodings->genAsin(_weights, _blockings, sum, tmpUnSATWeight,
                            compression, *_satSolver, lits, linkingVar);
        break;
      case OGAWA:
        _encodings->genOgawa0(_weights, _blockings, sum, tmpUnSATWeight,
                              divisor, compression, *_satSolver, lits,
                              linkingVar);
        break;
      case BAILLEUXW2:
        _encodings->genBailleuxW20(_weights, _blockings, sum, tmpUnSATWeight,
                                   compression, *_satSolver, lits, linkingVar,
                                   linkingWeight);
        break;
      case WMTO:
        _encodings->genKWMTO0(_weights, _blockings, sum, tmpUnSATWeight,
                              divisors, *_satSolver, lits, linkingVars,
                              linkingWeights);
        break;
      case MRWTO:
        _encodings->genMRWTO0(_weights, _blockings, sum, tmpUnSATWeight,
                              divisors, *_satSolver, lits, linkingVars,
                              linkingWeights, _encoding);
        break;
      case MRWTO2:
        _encodings->genMRWTO0(_weights, _blockings, sum, tmpUnSATWeight,
                              divisors, *_satSolver, lits, linkingVars,
                              linkingWeights, _encoding);
        break;
      case DGPW18:
        std::cout << "DGPW18 NOT SUPPORTED!" << std::endl;
        break;
      default:
        std::cout << "Encoding not yet defined!";
        exit(EXIT_FAILURE);
    }
  }
}

unsigned Pacose::SolveQMax(std::vector<SoftClause *> *tmpSoftClauses,
                           EncodingType *encodingType) {
  if (_settings.verbosity > 5) std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (tmpSoftClauses == nullptr) tmpSoftClauses = &_originalSoftClauses;

  if (encodingType == nullptr) encodingType = &_encoding;

  if (*encodingType == DGPW18 || *encodingType == HEURISTICQMAXSAT) {
    std::cout << "This encoding cannot be solved with QMaxSAT!" << std::endl;
    return UNKNOWN;
  }

  _nbOfOrigVars = _satSolver->GetNumberOfVariables();

  unsigned long localSumOfSoftWeights = 0;
  for (int i = 0; i < tmpSoftClauses->size(); i++) {
    localSumOfSoftWeights += (*tmpSoftClauses)[i]->weight;
    //    std::cout << (*tmpSoftClauses)[i]->weight << " ";
  }
  _encodings = new Encodings(&_settings);

  // Heuristic to choose which compression of warners!
  int compression = _settings._compression;
  if (compression == -1) {
    assert(*encodingType == WARNERS);
    if (tmpSoftClauses->size() > 1500) {
      compression = 1;
    } else {
      compression = 99;
    }
    if (_settings.verbosity > 0)
      std::cout << "c compressionRate set..: " << compression << std::endl;
  }
  //  std::cout << std::endl;
  unsigned long answer = localSumOfSoftWeights;
  unsigned long oldanswer = localSumOfSoftWeights;
  if (_settings.verbosity > 0) {
    std::cout << "c localSumOfSoftWeights: " << localSumOfSoftWeights
              << std::endl;
    std::cout << "c Local No Soft Clauses: " << tmpSoftClauses->size()
              << std::endl;
  }

  // this vector is used to create clauses!
  std::vector<unsigned> lits;
  // loop count
  int lcnt = 0;
  std::vector<unsigned> linkingVar;
  std::vector<long long int> linkingWeight;  // uemura 20161202
  //    bool *mmodel = new bool[static_cast<unsigned long>(_nbOfOrigVars)];
  //    //uemura 20161128
  lits.clear();
  linkingVar.clear();
  linkingWeight.clear();

  long long int divisor = 1;  // koshi 2013.10.04

  std::vector<long long int>
      ndivisor;  // mrwto用の複数の基数を保存する変数 uemura 2016.12.05
  ndivisor.clear();

  std::vector<std::vector<unsigned>>
      linkingVarMR;  // uemura 2016.12.05 for mrwto
  linkingVarMR.clear();
  std::vector<std::vector<long long int>>
      linkingWeightMR;  // uemura 2016.12.05 for mrwto
  linkingWeightMR.clear();

  std::vector<long long int> cc;  // cardinality constraints
  cc.clear();
  int ccSizeOld = -1;

  // koshi 20140701        lbool ret = S.solveLimited(dummy);
  unsigned ret = UNKNOWN;
  _encodings->_relaxLit = 0;

  while ((ret = _satSolver->Solve()) == SATISFIABLE) {  // koshi 20140107

    //    std::cout << "c noOfClauses AfterSolve: " <<
    //    _satSolver->GetNumberOfClauses()
    //              << std::endl;
    if (_encodings->_relaxLit != 0) {
      _satSolver->ClearAssumption();
      _satSolver->ResetClause();
      _satSolver->NewClause();
      unsigned rl = _encodings->_relaxLit ^ 1;
      _satSolver->AddLiteral(&rl);
      _satSolver->CommitClause();
      //      std::cout << "c RELAXLIT added: " << (_encodings->_relaxLit ^ 1)
      //                << std::endl;
      _encodings->_relaxLit = 0;
    }
    assert(_satSolver->Solve() == 10);
    //        std::cout << "ret: " << ret << std::endl;
    lcnt++;

    //    for (unsigned i = 0; i < _actualSoftClauses->size(); i++) {
    //      int varnum = (*_actualSoftClauses)[i]->relaxationLit >> 1;
    //      if ((*_actualSoftClauses)[i]->relaxationLit & 1) {
    //        if (!_satSolver->GetModel(varnum)) {
    //          answerNew += (*_actualSoftClauses)[i]->weight;
    //        }
    //      } else {
    //        if (_satSolver->GetModel(varnum)) {
    //          answerNew += (*_actualSoftClauses)[i]->weight;
    //        }
    //      }
    //    }
    long long int answerNew;
    if (lcnt != 1) {
      CalculateSATWeight();
      answerNew = CalculateLocalSATWeight();
      answerNew = localSumOfSoftWeights - answerNew;
    } else {
      answerNew = _localUnSatWeight;
    }

    if (lcnt > 2) {
      assert(answerNew < answer);
    }
    if (lcnt == 1 &&
        answerNew != 0) {  // first model: generate cardinal constraints

      unsigned ncls = _satSolver->GetNumberOfClauses();

      //      printf("c linkingVar.size() = %zu\n", linkingVar.size());
      // uemura 20161128
      genCardinals(answerNew, divisor, lits, linkingVar, linkingWeight,
                   ndivisor, linkingVarMR, linkingWeightMR, compression);

      localSumOfSoftWeights = 0;
      for (int i = 0; i < tmpSoftClauses->size(); i++) {
        localSumOfSoftWeights += (*tmpSoftClauses)[i]->weight;
        //    std::cout << (*tmpSoftClauses)[i]->weight << " ";
      }
      //      printf("c linkingVar.size() = %zu\n", linkingVar.size());
      //      std::cout << "Clauses before: " << ncls << std::endl;
      //      std::cout << "Clauses after: " << _satSolver->GetNumberOfClauses()
      //                << std::endl;
      //      std::cout << "Clauses (new): " << _satSolver->GetNumberOfClauses()
      //      - ncls
      //                << std::endl;
      if (_settings.verbosity > 0)
        std::cout << "c Cardinals are generated!" << std::endl;

      if (_settings.verbosity > 0) {
        // uemura 20161129
        if (_encoding == WMTO || _encoding == MRWTO || _encoding == MRWTO2) {
          printf("c ");
          for (int i = 0; i < linkingVarMR.size(); i++) {
            printf("c linkingVar[%d].size = %zu, ", i, linkingVarMR[i].size());
          }
          printf("\n");
        } else {
          printf("c linkingVar.size() = %zu\n", linkingVar.size());
        }
        printf("c Cardinality Constraints: %d variables and %d clauses\n",
               _satSolver->GetNumberOfVariables() - _nbOfOrigVars,
               _satSolver->GetNumberOfClauses() - ncls);
      }
    }

    if (answerNew > 0) {
      unsigned long nofCl = _satSolver->GetNumberOfClauses();
      if (_encoding == WMTO || _encoding == MRWTO || _encoding == MRWTO2) {
        _encodings->lessthanMR(linkingVarMR, linkingWeightMR, answer, answerNew,
                               ndivisor, cc, *_satSolver, lits, _encoding);
      } else {
        if (_encoding == BAILLEUX && lcnt == 1) answer = linkingVar.size();

        ccSizeOld = cc.size();
        //        std::cout << "ccSizeOld: " << cc.size() << std::endl;

        _encodings->lessthan(linkingVar, linkingWeight, answer, answerNew,
                             divisor, cc, *_satSolver, _encoding);
      }
      if (_settings.verbosity > 0)
        std::cout << "c Clauses of encoding: "
                  << _satSolver->GetNumberOfClauses() - nofCl << std::endl;
      //      std::cout << "ccSizeNew: " << cc.size() << std::endl;
      oldanswer = answer;
      answer = answerNew;
    } else {
      answer = answerNew;
      ret = UNSAT;  // koshi 20140124
      break;
    }
    //    std::cout << "noOfClauses BeforeNextSolve: "
    //              << _satSolver->GetNumberOfClauses() << std::endl;
    //    std::cout << "oldAnswer: " << oldanswer << "  answer: " << answer
    //              << std::endl;
  }

  // koshi 20140124
  if (ret == UNSAT) {
    if (lcnt > 0) {
      if (_settings.verbosity > 0) printf("c local opt found\n");

      if (answer != 0 && oldanswer == answer + 1 && lcnt > 1) {
        assert(_satSolver->Solve() == 20);
        _satSolver->ClearAssumption();
        assert(_satSolver->Solve() == 10);
        //        unsigned lastResult = _satSolver->Solve();
        //        std::cout << "SolveBefore: " << lastResult <<
        //        std::endl;
      } else if (answer != 0) {
        assert(_satSolver->Solve() == 20);
        _satSolver->ClearAssumption();
        assert(_satSolver->Solve() == 10);
        //        std::cout << "SOLVE RESULT: " << _satSolver->Solve() <<
        //        std::endl;
        // deactivate last assumption!
        //        _satSolver->ClearAssumption();
        //        unsigned lastResult = _satSolver->Solve();
        //        std::cout << "SolveBefore: " << lastResult << std::endl;
        //        assert(lastResult == SATISFIABLE);
        if (_encodings->_relaxLit != 0) {
          _satSolver->ResetClause();
          _satSolver->NewClause();
          unsigned rl = _encodings->_relaxLit;
          //        std::cout << "AddRL: " << _encodings->_relaxLit <<
          //        std::endl;
          _satSolver->AddLiteral(&rl);
          _satSolver->CommitClause();
        }
        _encodings->_relaxLit = 0;
        //        std::cout << "SolveBefore2: " << _satSolver->Solve() <<
        //        std::endl; printf("local o %lld\n",
        //               localSumOfSoftWeights - CalculateLocalSATWeight());
        //        std::cout << "noOfClauses: " <<
        //        _satSolver->GetNumberOfClauses()
        //                  << std::endl;
        if (cc.size() > ccSizeOld) {
          if (_settings.verbosity > 0)
            std::cout << "c ccSize is going to be decreased Old/New: "
                      << ccSizeOld << ", " << cc.size() << std::endl;
          while (cc.size() > ccSizeOld) {
            cc.pop_back();
          }
        }
        //        std::cout << "ccSizeOld: " << cc.size() << std::endl;
        _satSolver->ClearAssumption();
        unsigned long nofCl = _satSolver->GetNumberOfClauses();
        if (_encoding == WMTO || _encoding == MRWTO || _encoding == MRWTO2) {
          _encodings->lessthanMR(linkingVarMR, linkingWeightMR, oldanswer,
                                 answer + 1, ndivisor, cc, *_satSolver, lits,
                                 _encoding);
        } else {
          if (_encoding == BAILLEUX && lcnt == 1) answer = linkingVar.size();

          _encodings->lessthan(linkingVar, linkingWeight, oldanswer, answer + 1,
                               divisor, cc, *_satSolver, _encoding);
        }
        //        std::cout << "Newly Added clauses: "
        //                  << _satSolver->GetNumberOfClauses() - nofCl <<
        //                  std::endl;
        //        //        std::cout << "ccSizeNew: " << cc.size() <<
        //        std::endl; unsigned lastResult = _satSolver->Solve();

        //        std::cout << "SolveAfterWithNewAssumptions: " << lastResult
        //                  << std::endl;
        assert(_satSolver->Solve() == SATISFIABLE);
        if (_encodings->_relaxLit != 0) {
          _satSolver->ResetClause();
          _satSolver->NewClause();
          unsigned rl = _encodings->_relaxLit ^ 1;

          _satSolver->AddLiteral(&rl);
          //        std::cout << "AddRL: " << _encodings->_relaxLit <<
          //        std::endl;
          _satSolver->CommitClause();
        }
        _encodings->_relaxLit = 0;
        //        CalculateLocalSATWeight();
        _satSolver->ClearAssumption();
        //        lastResult = _satSolver->Solve();
        //        std::cout << "SolveAfter2: " << lastResult << std::endl;
        assert(_satSolver->Solve() == SATISFIABLE);
        //        CalculateLocalSATWeight();
        //        CalculateSATWeight();
      } else {
        if (_settings.verbosity > 0)
          std::cout << "c ANSWER IS 0!, Add all SCs" << std::endl;
        //        unsigned lastResult = _satSolver->Solve();
        //        std::cout << "SolveAfter2: " << lastResult << std::endl;
        //        assert(lastResult == SATISFIABLE);

        //        case answer == 0 -- all SCs are SAT
        for (auto SC : *_actualSoftClauses) {
          _satSolver->ResetClause();
          _satSolver->NewClause();
          unsigned rlit = SC->relaxationLit ^ 1;
          _satSolver->AddLiteral(&rlit);
          _satSolver->CommitClause();
        }
        _satSolver->ClearAssumption();
        //        lastResult = _satSolver->Solve();
        //        std::cout << "SolveAfterAddingAllRelaxLits: " << lastResult
        //                  << std::endl;
        //        assert(lastResult == SATISFIABLE);
      }
    } else {
      if (_settings.verbosity > 0) printf("s Hard clauses are UNSATISFIABLE\n");
    }
  } else if (ret == SATISFIABLE) {
    std::cout << "c ERROR: SHOULD NEVER OCCUR!" << std::endl;
  } else if (_settings.verbosity > 0) {
    printf("s UNKNOWN\n");
    printf("c Search is stopped by a limit (maybe time-limit)\n");
  }
  //  if (lcnt > 0) {
  //    //        if (_settings.GetPrintModel()) {
  //    //            printf("v ");
  //    //            for (int i = 0; i< _nbOfOrigVars; i++) {
  //    //                printf("%s%d ", mmodel[i]?"":"-", i+1);
  //    //                if ((i+1)%20 == 0 && i+1 < _nbOfOrigVars) printf("\nv
  //    ");
  //    //            }
  //    //            printf("\n");
  //    //        }
  //    std::cout << std::endl;
  //    if (_settings.GetPrintModel()) {
  //      std::cout << "c PrintModel is active!" << std::endl;
  //      //            PrintModel();
  //    }
  //  }

  if (lcnt > 0 && _settings.verbosity > 0)
    std::cout << "c Latest Answer = " << answer << "/" << localSumOfSoftWeights
              << "  loops " << lcnt << std::endl;

  if (_settings.verbosity > 5) std::cout << "END QMAXSAT" << std::endl;
  //    delete[] mmodel;
  return ret;
}

bool Pacose::AddEncoding(std::vector<SoftClause *> *tmpSoftClauses,
                         EncodingType *encodingType) {
  if (tmpSoftClauses == nullptr) {
    tmpSoftClauses = &_originalSoftClauses;
  }
  if (encodingType == nullptr) {
    encodingType = &_encoding;
  }
  if (*encodingType == DGPW18) {
    return true;
  } else if (*encodingType != HEURISTICQMAXSAT) {
    return true;
  } else {
    return false;
  }

  return true;
}

unsigned Pacose::SolveMaxSAT(std::vector<SoftClause *> *tmpSoftClauses,
                             EncodingType *encodingType) {
  if (tmpSoftClauses == nullptr) {
    tmpSoftClauses = &_originalSoftClauses;
  }
  if (encodingType == nullptr) {
    encodingType = &_encoding;
  }
  if (*encodingType == DGPW18) {
    //    return _cascCandidates[i - 1].dgpw->MaxSolveWeightedPartial(optimum);
  } else if (*encodingType != HEURISTICQMAXSAT) {
    SolveQMax();
  } else {
    return false;
  }

  return true;
}

unsigned Pacose::SolveProcedure(long long int nbVars) {
  _settings.formulaIsDivided = true;
  double timeStart;
  struct rusage resources;
  getrusage(RUSAGE_SELF, &resources);
  timeStart =
      resources.ru_utime.tv_sec + 1.e-6 * (double)resources.ru_utime.tv_usec;
  std::vector<unsigned> lastSatisfiableAssignment = {};

  DivideSCsIfPossible();
  DGPW::DGPW *dgpw = nullptr;

  if (nbVars == 0) {
    _nbVars = _satSolver->GetNumberOfVariables();
  } else {
    _nbVars = nbVars + 1;
  }

  std::vector<std::pair<uint64_t, uint32_t>> tares = {};
  std::vector<std::pair<uint64_t, uint32_t>> watchdogs = {};
  unsigned variablesOfEncoding = 0;
  unsigned clausesOfEncoding = 0;
  bool greedyPreproUsed = false;
  unsigned noUnsatSCs = 0;
  unsigned noApproxSolverCalls = 0;
  unsigned noSolverCalls = 0;
  double timeTrimming = 0;
  _settings.Print();

  if (_satSolver->Solve() == 20) {
    std::cout << "c Hard Clauses are not Satisfiable!" << std::endl;
    std::cout << "UNSATISFIABLE" << std::endl;
    return 20;
  }

  _encoding = _settings._encoding;
  if (_encoding == HEURISTIC20) {
    if (_settings.verbosity > 1) std::cout << "c HEURISTIC20" << std::endl;
    if (((_overallSoftWeights / _originalSoftClauses.size()) > 850) &&
        (_overallSoftWeights < 80000000000) &&
        (_originalSoftClauses.size() < 50000)) {
      _encoding = WARNERS;
      std::cout << "c Use Warners encoding!" << std::endl;
    } else {
      _encoding = DGPW18;
      std::cout << "c Use DGPW encoding!" << std::endl;
    }
  }

  for (unsigned i = static_cast<unsigned>(_sClauses.size()); i > 0; i--) {
    _settings.currentCascade.iteration = i - 1;
    _actualSoftClauses = &_sClauses[i - 1];

    if (_settings.simplify && _satSolver->Simplify() != 10) {
      std::cout << "c ERROR Instance not Satisfiable anymore, should never "
                   "happen - failure in encoding! - Return UNSAT!"
                << std::endl;
      return 20;
    }

    if (_encoding == DGPW18) {
      // for a kind of stratification only
      std::vector<unsigned> sClause;
      for (auto tare : tares) {
        sClause.push_back(tare.second);
        AddSoftClauseTo(_actualSoftClauses, sClause, tare.first);
        //            std::cout << "weight: " << tare.first << "  tare: " <<
        //            tare.second << std::endl;
        sClause.clear();
      }

      for (auto watchdog : watchdogs) {
        sClause.push_back(watchdog.second);
        AddSoftClauseTo(_actualSoftClauses, sClause, watchdog.first);
        //            std::cout << "weight: " << watchdog.first << "  watchdog:
        //            "
        //            << watchdog.second << std::endl;
        sClause.clear();
      }
    }

    // calc GCD
    // convert into MaxSAT if possible!
    Preprocess();

    unsigned long sumOfActualWeights = 0;
    for (auto sc : *_actualSoftClauses) {
      sumOfActualWeights += sc->weight;
    }
    uint64_t localSatWeight = CalculateLocalSATWeight();
    _localUnSatWeight = sumOfActualWeights - localSatWeight;
    if (_settings.verbosity > 0) {
      std::cout << "c local o: " << _localUnSatWeight << std::endl;
      if (_settings.verbosity > 3) {
        std::cout << "c number of SCs: " << _actualSoftClauses->size()
                  << std::endl;
        std::cout << "c sumOfActualWeights: " << sumOfActualWeights
                  << std::endl;
        std::cout << "c locSatWeight: " << localSatWeight << std::endl;
      }
    }

    wbSortAndFilter(_localUnSatWeight);

    int fixSCs = _settings.greedyPPFixSCs;
    int minSizeSCs = _settings.greedyMinSizeOfSet;

    // heuristic to choose when to use GreedyPreProcessing
    if (fixSCs == -1) {
      if (_encoding == DGPW18) {
        minSizeSCs = 650;
        fixSCs = 2;
      } else if (_encoding == WARNERS) {
        minSizeSCs = 12000;
        fixSCs = 2;
      } else {
        std::cout << "ERROR: no heuristic set yet!" << std::endl;
        assert(true);
        exit(1);
      }
      if (_settings.verbosity > 0)
        std::cout << "c MinSize for PP chose...: " << minSizeSCs << std::endl;
    }

    if (_settings.greedyPrepro != 0 &&
        _actualSoftClauses->size() > minSizeSCs) {
      greedyPreproUsed = true;
      double tmpTimeTrimming;
      struct rusage resources;
      getrusage(RUSAGE_SELF, &resources);
      tmpTimeTrimming = resources.ru_utime.tv_sec +
                        1.e-6 * (double)resources.ru_utime.tv_usec;

      GreedyPrepro *greedyPrePro = new GreedyPrepro(
          *_actualSoftClauses, &_settings, _satSolver, this, fixSCs);
      //      greedyPrePro->AddClausesToSATSolver(*_CLInterface.GetClauseDB(),
      //                                          _CLInterface.GetNoVars());
      _localUnSatWeight = greedyPrePro->StartPrepro();

      _satSolver->Solve();
      CalculateLocalSATWeight();
      //      if (greedyPrePro->GetUnsatWeight() == 0) {
      //        std::cout << "c UNSATWEIGHT == 0" << std::endl;
      //        for (auto sc : *_actualSoftClauses) {
      //          std::vector<unsigned> cl;

      //          for (auto lit : sc->clause) {
      //            unsigned literal = lit;
      //            cl.push_back(literal);
      //          }
      //          _CLInterface.AddClause(cl);
      //        }

      getrusage(RUSAGE_SELF, &resources);
      tmpTimeTrimming = resources.ru_utime.tv_sec +
                        1.e-6 * (double)resources.ru_utime.tv_usec -
                        tmpTimeTrimming;
      timeTrimming += tmpTimeTrimming;

      //      sumOfActualWeights = 0;
      //      for (auto sc : *_actualSoftClauses) {
      //        sumOfActualWeights += sc->weight;
      //      }

      //      sumOfActualWeights - CalculateLocalSATWeight();
      if (_settings.verbosity > 0)
        std::cout << "c local PrePro o: " << _localUnSatWeight << std::endl;
      //      wbSortAndFilter(_localUnSatWeight);
    }
    //    else if (i == 1 && !greedyPreproUsed) {
    //      std::cout << "c Never Used Greedy Prepro, kill yourself!" <<
    //      std::endl; exit(1);
    //    }

    if (_settings.simplify && _satSolver->Simplify() != 10) {
      std::cout << "c ERROR: Instance not Satisfiable anymore, return UNSAT!"
                << std::endl;
      exit(1);
    }

    if (_actualSoftClauses->size() == 0) {
      continue;
    } else if (_localUnSatWeight == 0) {
      // case all SCs are already SAT
      for (auto SC : *_actualSoftClauses) {
        unsigned relaxLit = SC->relaxationLit ^ 1;
        _satSolver->ResetClause();
        _satSolver->NewClause();
        _satSolver->AddLiteral(&relaxLit);
        _satSolver->CommitClause();
        _satSolver->ClearAssumption();
      }
      assert(_satSolver->Solve() == 10);
      continue;
    } else if (_actualSoftClauses->size() == 1) {
      // border case - only one SC
      // which is not yet SAT otherwise localUnsatWeight would be 0
      _satSolver->ClearAssumption();
      unsigned relaxLit = (*_actualSoftClauses)[0]->relaxationLit ^ 1;
      _satSolver->AddAssumption(&relaxLit);
      _satSolver->ResetClause();
      _satSolver->NewClause();
      if (_satSolver->Solve() == 10) {
        _satSolver->AddLiteral(&relaxLit);

      } else {
        relaxLit = relaxLit ^ 1;
        _satSolver->AddLiteral(&relaxLit);
      }
      _satSolver->CommitClause();
      _satSolver->ClearAssumption();
      continue;
    } /*else if (_localUnSatWeight == 1) {
      // there is only a single soft clause which is unsatisfiable
      _satSolver->ClearAssumption();
      // check if it is possible to solve all soft clauses
      _satSolver->AddAssumption();
      // if not add model
    }*/

    if (_encoding == DGPW18) {
      _cascCandidates[i - 1].dgpw = new DGPW::DGPW(this);

      if (i < static_cast<unsigned>(_sClauses.size()) && tares.empty() &&
          watchdogs.empty() && i >= 1 &&
          _cascCandidates[i - 1].weightsTillPoint >=
              _cascCandidates[i].ggtTillPoint &&
          _encoding == DGPW18) {
        lastSatisfiableAssignment = GetBestSCAssignment();
        //        assert(_cascCandidates[i].dgpw->Solve(lastSatisfiableAssignment)
        //        ==
        //               SATISFIABLE);
        _cascCandidates[i - 1].dgpw->SetInitialAssumptions(
            lastSatisfiableAssignment);
      }

      tares.clear();
      watchdogs.clear();

      _cascCandidates[i - 1].dgpw->SetSoftClauseVector(_actualSoftClauses);
      if (_settings.verbosity > 0)
        std::cout << "c greatest Common Divisor: " << _GCD << std::endl;
      _cascCandidates[i - 1].dgpw->SetGreatestCommonDivisor(
          static_cast<uint64_t>(_GCD));

      if (i > 1 && _cascCandidates[i - 2].weightsTillPoint >
                       _cascCandidates[i - 1].ggtTillPoint) {
        if (_settings.divideDGPW == DIVIDEALLSOLVEONLYWATCHDOGS && i > 1 &&
            _cascCandidates[i - 2].weightsTillPoint >=
                _cascCandidates[i - 1].ggtTillPoint) {
          _settings.currentCascade._solveTares = false;
        } else {
          _settings.currentCascade._solveTares = true;
        }
        _settings.currentCascade._onlyWithAssumptions = true;
      } else {
        _settings.currentCascade._solveTares = true;
        _settings.currentCascade._onlyWithAssumptions = false;
      }
    }
    //        assert(_satSolver->Solve() == 10);
    //        assert(_cascCandidates[i-1].dgpw->Solve() == 10);
    //        std::cout << _cascCandidates[i-1].dgpw->Solve() << std::endl;
    //    std::cout << "ASC: ";
    //    for (auto SC : *_actualSoftClauses) {
    //      std::cout << SC->weight << " ";
    //    }
    //    std::cout << std::endl;
    int64_t optimum = -1;
    if (_encoding == DGPW18 || _actualSoftClauses->size() == 1) {
      //        || i == 2) {
      unsigned long sumOfActualWeights = 0;
      for (auto sc : *_actualSoftClauses) {
        sumOfActualWeights += sc->weight;
      }
      //      _localUnSatWeight = sumOfActualWeights -
      //      CalculateLocalSATWeight();
      _cascCandidates[i - 1].dgpw->SetSatWeight(sumOfActualWeights -
                                                _localUnSatWeight);
      _cascCandidates[i - 1].dgpw->MaxSolveWeightedPartial(optimum);
      clausesOfEncoding += _cascCandidates[i - 1].dgpw->GetEncodingClauses();
      variablesOfEncoding +=
          _cascCandidates[i - 1].dgpw->GetEncodingVariables();
      noUnsatSCs += _cascCandidates[i - 1].dgpw->GetNoUnsatisfiableSCs();
      noApproxSolverCalls +=
          _cascCandidates[i - 1].dgpw->GetApproxSolverCalls();
      noSolverCalls += _cascCandidates[i - 1].dgpw->GetSolverCalls();
      //      if (_cascCandidates[i - 1].dgpw->GetGreedyPrepro() != 0)
      //        noGreedyPreproUsed++;

    } else if (_encoding != HEURISTICQMAXSAT) {
      //      std::cout << std::endl;
      //      if (i == 2) {
      //        std::cout << "ADD ALL RELAXLITS!" << std::endl;

      //        for (auto sc : *_actualSoftClauses) {
      //          _satSolver->ResetClause();
      //          _satSolver->NewClause();
      //          unsigned rlit = sc->relaxationLit ^ 1;
      //          //          std::cout << rlit << std::endl;
      //          _satSolver->AddLiteral(&rlit);
      //          _satSolver->CommitClause();
      //        }

      //        _satSolver->ClearAssumption();
      //        unsigned lastResult = _satSolver->Solve();
      //        std::cout << "SolveAfterAddingAllRelaxLits: " << lastResult
      //                  << std::endl;
      //        CalculateSATWeight();

      //      } else {
      _satSolver->ClearAssumption();
      SolveQMax(_actualSoftClauses, &_encoding);
      //      }
    }

    _satSolver->ClearAssumption();
    unsigned cr = _satSolver->Solve();
    if (cr != 10) {
      std::cout << "c ERROR: STRANGE SOLVING RESULT " << cr << " -> QUIT!"
                << std::endl;
      exit(0);
    }
    CalculateSATWeight();
    //    std::cout
    //        << "                                              Local SAT
    //        weight: "
    //        << std::endl;
    //    CalculateLocalSATWeight();

    //        std::cout << "c current Cascade iter...: " <<
    //        static_cast<unsigned>(_sClauses.size()) - i
    //                  << std::endl;
    unsigned long sumAllLowerWeights = 0;
    if (i > 1) {
      sumAllLowerWeights = _cascCandidates[i - 2].weightsTillPoint;
    }

    unsigned long weightOfTaresAndWatchdogs = 0;
    if (i != 1 && _cascCandidates[i - 2].weightsTillPoint >
                      _cascCandidates[i - 1].ggtTillPoint) {
      std::vector<unsigned> assumptions =
          _cascCandidates[i - 1].dgpw->GetLastAssumptions();

      unsigned long maxGreedyWeight = 0;
      if (_settings.useGreedyPreInBetween &&
          _settings.divideDGPW == DIVIDEALL) {
        // greedy Solving next cascade!
        dgpw = new DGPW::DGPW(this);
        std::vector<SoftClause *> tmpSCs = {};
        for (unsigned j = i - 2; j + 1 > 0; j--) {
          tmpSCs.insert(tmpSCs.end(), _sClauses[j].begin(), _sClauses[j].end());
        }
        //            std::cout << tmpSCs.size() << std::endl;

        //            *_actualSoftClauses = _sClauses[i - 2];
        dgpw->SetGreatestCommonDivisor(1);
        dgpw->SetSoftClauseVector(&tmpSCs);
        dgpw->SetFixedAssumptions(assumptions);

        _settings.currentCascade._solveTares = false;
        dgpw->MaxSolveWeightedPartial(optimum);
        clausesOfEncoding += dgpw->GetEncodingClauses();
        variablesOfEncoding += dgpw->GetEncodingVariables();
        noUnsatSCs += dgpw->GetNoUnsatisfiableSCs();
        dgpw->RemoveFixedAssumptions();
        maxGreedyWeight = dgpw->GetSATWeight();
        //                dgpw->~DGPW();
        if (_settings.verbosity > 0)
          std::cout << "c GreedyPrepro Max Weight: " << maxGreedyWeight
                    << std::endl;
      }

      unsigned long weightDiff = sumAllLowerWeights - maxGreedyWeight;

      if (_settings.verbosity > 1) {
        std::cout << "c sumAllLowerWeights: " << sumAllLowerWeights
                  << std::endl;
        std::cout << "c NEXT MAX CASC WEIGHT DIFF: " << weightDiff << std::endl;
      }

      tares = _cascCandidates[i - 1].dgpw->GetTareVector(weightDiff);
      watchdogs = _cascCandidates[i - 1].dgpw->GetWatchdogs(weightDiff);

      for (auto tare : tares) {
        weightOfTaresAndWatchdogs += tare.first;
        //                std::cout << "c weight: " << tare.first << "  tare: "
        //                << tare.second << std::endl;
      }

      for (auto watchdog : watchdogs) {
        weightOfTaresAndWatchdogs += watchdog.first;
        //                std::cout << "c weight: " << watchdog.first << "
        //                watchdog: " << watchdog.second << std::endl;
      }
      //      std::cout << std::endl;
    }
  }

  CalculateSATWeight();

  if (_settings.verbosity > 0) {
    std::cout << "c #SolverCalls approx....: " << noApproxSolverCalls
              << std::endl;
    std::cout << "c #SolverCalls...........: " << noSolverCalls << std::endl;
    std::cout << "c trimming used..........: " << greedyPreproUsed << std::endl;
    std::cout << "c #unsatisfiable SCs.....: " << noUnsatSCs << std::endl;
    std::cout << "c #clauses of encoding...: " << clausesOfEncoding
              << std::endl;
    std::cout << "c #variables of encoding.: " << variablesOfEncoding
              << std::endl;
    getrusage(RUSAGE_SELF, &resources);
    std::cout << "c time...................: "
              << (resources.ru_utime.tv_sec +
                  1.e-6 * (double)resources.ru_utime.tv_usec) -
                     timeStart
              << std::endl;
  }

  std::cout << "o " << _unSatWeight << std::endl;
  std::cout << "s OPTIMUM FOUND" << std::endl;

  PrintResult();
  return 10;
}

void Pacose::PrintResult() {
  if (!_settings._printModel) return;

  //  old print model version std::cout << "v ";
  //  for (unsigned i = 1; i < _nbVars; i++) {
  //    std::cout << (_satSolver->GetModel(i)) << " ";
  //  }
  //  std::cout << std::endl;

  std::cout << "v ";
  for (unsigned i = 1; i < _nbVars; i++) {
    std::cout << ((_satSolver->GetModel(i) ^ 1) % 2);
  }
  std::cout << std::endl;
}

void Pacose::CalculateOverallTimes() {
  if (_settings.verbosity > 4) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
  }
  DGPW::TimeVariables overallTime;
  int iter = 0;
  for (auto casc : _cascCandidates) {
    iter++;
    overallTime.AddTimeStruct(casc.dgpw->_timeVariables);
    casc.dgpw->_timeVariables->DumpVariables(iter);
  }
  overallTime.DumpVariables();
}

std::vector<unsigned> Pacose::GetBestSCAssignment() {
  std::vector<unsigned> SCModel;
  for (auto SC : _originalSoftClauses) {
    SCModel.push_back(SC->lastassignment);
  }
  return SCModel;
}

uint64_t Pacose::CalculateLocalSATWeight(
    std::vector<SoftClause *> *tmpSoftClauses) {
  if (_settings.verbosity > 4) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << "c Clauses.size: " << _satSolver->GetNumberOfClauses()
              << std::endl;
  }

  if (tmpSoftClauses == nullptr) {
    tmpSoftClauses = _actualSoftClauses;
  }

  long long int unSatWeight = 0;
  long long int satWeight = 0;

  // Process all soft clauses
  for (uint32_t i = 0; i != tmpSoftClauses->size(); ++i) {
    uint32_t relaxlit = (*tmpSoftClauses)[i]->relaxationLit;
    if (_satSolver->GetModel(relaxlit >> 1) == relaxlit) {
      std::vector<uint32_t> clause((*tmpSoftClauses)[i]->clause);
      uint32_t pos = 0;
      for (; pos != clause.size(); ++pos) {
        // clause satisfied without trigger?
        if (_satSolver->GetModel(clause[pos] >> 1) == clause[pos]) {
          satWeight += (*tmpSoftClauses)[i]->weight;
          (*tmpSoftClauses)[i]->lastassignment = relaxlit ^ 1;
          //                    std::cout << "++ " << i << ": " << i << "/"
          //                              <<
          //                              _originalSoftClauses[i]->originalWeight
          //                              << std::endl;
          break;
        }
      }
      if (pos == clause.size()) {
        (*tmpSoftClauses)[i]->lastassignment = relaxlit;
        unSatWeight += (*tmpSoftClauses)[i]->weight;
      }
    } else if (_satSolver->GetModel(relaxlit >> 1) != 0) {
      assert(_satSolver->GetModel(relaxlit >> 1) == (relaxlit ^ 1));
      (*tmpSoftClauses)[i]->lastassignment = relaxlit ^ 1;
      satWeight += (*tmpSoftClauses)[i]->weight;
      //            std::cout << " ++relax sat " << i << ": " << i << "/"
      //                      << _originalSoftClauses[i]->originalWeight <<
      //                      std::endl;
    }
  }

  if (_settings.verbosity > 0) {
    std::cout << "c new local max found: " << satWeight << std::endl;
    std::cout << "c new local o found: " << unSatWeight << std::endl;
  }

  return static_cast<uint64_t>(satWeight);
}

uint64_t Pacose::CalculateSATWeight() {
  if (_settings.verbosity > 4) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << "c SC.size: " << _originalSoftClauses.size() << std::endl;
    std::cout << "c Clauses.size: " << _satSolver->GetNumberOfClauses()
              << std::endl;
  }

  long long int unSatWeight = 0;
  long long int satWeight = 0;

  // Process all soft clauses
  for (uint32_t i = 0; i != _originalSoftClauses.size(); ++i) {
    uint32_t relaxlit = _originalSoftClauses[i]->relaxationLit;
    if (_satSolver->GetModel(relaxlit >> 1) == relaxlit) {
      std::vector<uint32_t> clause(_originalSoftClauses[i]->clause);
      uint32_t pos = 0;
      for (; pos != clause.size(); ++pos) {
        // clause satisfied without trigger?
        if (_satSolver->GetModel(clause[pos] >> 1) == clause[pos]) {
          satWeight += _originalSoftClauses[i]->originalWeight;
          _originalSoftClauses[i]->lastassignment = relaxlit ^ 1;
          //                    std::cout << "++ " << i << ": " << i << "/"
          //                              <<
          //                              _originalSoftClauses[i]->originalWeight
          //                              << std::endl;
          break;
        }
      }
      if (pos == clause.size()) {
        _originalSoftClauses[i]->lastassignment = relaxlit;
        unSatWeight += _originalSoftClauses[i]->originalWeight;
      }
    } else if (_satSolver->GetModel(relaxlit >> 1) != 0) {
      assert(_satSolver->GetModel(relaxlit >> 1) == (relaxlit ^ 1));
      _originalSoftClauses[i]->lastassignment = relaxlit ^ 1;
      satWeight += _originalSoftClauses[i]->originalWeight;
      //            std::cout << " ++relax sat " << i << ": " << i << "/"
      //                      << _originalSoftClauses[i]->originalWeight <<
      //                      std::endl;
    }
  }
  if (satWeight > _satWeight || _unSatWeight < unSatWeight) {
    _satWeight = satWeight;
    _unSatWeight = unSatWeight;
    if (_settings.verbosity > 0)
      std::cout << "c global SAT weight: " << _satWeight << std::endl;
    std::cout << "o " << _unSatWeight << std::endl;
  } else {
    if (_settings.verbosity > 0) {
      std::cout << "c new local max found: " << satWeight << std::endl;
      std::cout << "c new local o found: " << unSatWeight << std::endl;
    }
  }

  return static_cast<uint64_t>(_satWeight);
}

void Pacose::DumpCascCandidates() {
  if (_settings.verbosity == 0) return;
  //    int i = 0;
  //    for (auto cascade : _cascCandidates) {
  //        std::cout << "c Points[" << i << "]..............: " <<
  //        cascade.Points << std::endl; std::cout << "c weightsTillPoint[" << i
  //        << "]....: " << cascade.weightsTillPoint
  //                  << std::endl;
  //        std::cout << "c ggtTillPoint[" << i << "]........: " <<
  //        cascade.ggtTillPoint << std::endl; std::cout << "c
  //        allWeightsAreEqual[" << i << "]..: " << cascade.allWeightsAreEqual
  //                  << std::endl;

  //        i++;
  //    }

  int i = 0;
  for (auto cascade : _cascCandidates) {
    if (i == 0) {
      std::cout << "c Points[" << i << "]..............: " << cascade.Points
                << std::endl;
      std::cout << "c weightsTillPoint[" << i
                << "]....: " << cascade.weightsTillPoint << std::endl;
    } else {
      std::cout << "c Points[" << i << "]..............: "
                << _cascCandidates[i].Points - _cascCandidates[i - 1].Points
                << std::endl;
      std::cout << "c weightsTillPoint[" << i << "]....: "
                << _cascCandidates[i].weightsTillPoint -
                       _cascCandidates[i - 1].weightsTillPoint
                << std::endl;
    }
    std::cout << "c ggtTillPoint[" << i << "]........: " << cascade.ggtTillPoint
              << std::endl;
    std::cout << "c allWeightsAreEqual[" << i
              << "]..: " << cascade.allWeightsAreEqual << std::endl;

    i++;
  }
}

void Pacose::DivideSCs(std::vector<unsigned> &sortedSCs, int acceptedMode) {
  partitionInformation nextCandidate;
  if (_cascCandidates.size() == 0) {
    _cascCandidates.push_back(nextCandidate);
  } else {
    assert(_cascCandidates.size() == 1);
    _cascCandidates[0] = nextCandidate;
  }
  _cascCandidates.back().ggtTillPoint =
      _originalSoftClauses[sortedSCs[0]]->weight;

  unsigned j = 0;

  //    std::cout << "partitionPoints.size(): " << SCs.Points.size() <<
  //    std::endl;

  _overallSoftWeights = 0;
  //    unsigned long weightRange =
  //    _originalSoftClauses[sortedSCs.back()]->weight -
  //    _originalSoftClauses[sortedSCs[0]]->weight;

  //    unsigned long weightRange = 5;
  //    std::cout << "c weightRange: " << weightRange << std::endl;
  int mode = -1;

  for (unsigned i = 0; i < sortedSCs.size(); ++i) {
    j++;

    mode = -1;
    if (acceptedMode == 0 &&
        _overallSoftWeights <= _originalSoftClauses[sortedSCs[i]]->weight &&
        i > 0) {
      //      std::cout << "se" << _originalSoftClauses[sortedSCs[i]]->weight
      //                << std::endl;
      if (!(_settings.testIfDividable != 2 &&
            _overallSoftWeights ==
                _originalSoftClauses[sortedSCs[i]]->weight)) {
        //        std::cout << "mod0" << std::endl;
        mode = 0;
        // if it is only the first element, then undo the mode decision unless
        // it is a strict smaller!
        if (i == 1 &&
            _overallSoftWeights == _originalSoftClauses[sortedSCs[i]]->weight) {
          //          std::cout << "mod-1" << std::endl;
          mode = -1;
        }
      }
    }
    if (acceptedMode == 1 && j > _originalSoftClauses.size() / 15 &&
        i < _originalSoftClauses.size() - (_originalSoftClauses.size() / 15) &&
        (10 * _originalSoftClauses[sortedSCs[i - 1]]->weight) <
            _originalSoftClauses[sortedSCs[i]]->weight)
      mode = 1;
    if (acceptedMode == 2 && j > _originalSoftClauses.size() / 4 &&
        (_originalSoftClauses[sortedSCs[i]]->weight >
         _originalSoftClauses[sortedSCs[i - 1]]->weight * 3) &&
        i < _originalSoftClauses.size() - (_originalSoftClauses.size() / 4) &&
        (2 * _overallSoftWeights < _sumOfSoftWeights - _overallSoftWeights))
      mode = 2;
    if (acceptedMode == 3 && j > _originalSoftClauses.size() / 6 &&
        i < _originalSoftClauses.size() - (_originalSoftClauses.size() / 6) &&
        _originalSoftClauses[sortedSCs[i - 1]]->weight <
            _originalSoftClauses[sortedSCs[i]]->weight &&
        _originalSoftClauses[sortedSCs[i - 1]]->weight ==
            _cascCandidates.back().ggtTillPoint)
      mode = 3;
    if (acceptedMode == 4 && j > _originalSoftClauses.size() / 4 &&
        _originalSoftClauses.size() > 150 &&
        (_originalSoftClauses[sortedSCs[i]]->weight >
         _originalSoftClauses[sortedSCs[i - 1]]->weight * 1.5) &&
        i < _originalSoftClauses.size() - (_originalSoftClauses.size() / 4) &&
        (5 * _overallSoftWeights < _sumOfSoftWeights - _overallSoftWeights))
      mode = 4;
    if (acceptedMode == 5 && j > 20 &&
        (_originalSoftClauses[sortedSCs[i]]->weight >
         _originalSoftClauses[sortedSCs[i - 1]]->weight * 20) &&
        i < _originalSoftClauses.size() - 20 &&
        _overallSoftWeights < 2 * (_sumOfSoftWeights - _overallSoftWeights))
      mode = 5;
    if (acceptedMode == 6 && j > 20 &&
        _cascCandidates.back().allWeightsAreEqual &&
        (_originalSoftClauses[sortedSCs[i]]->weight >
         _originalSoftClauses[sortedSCs[i - 1]]->weight * 10) &&
        i < _originalSoftClauses.size() - 20 &&
        _overallSoftWeights < 5 * (_sumOfSoftWeights - _overallSoftWeights))
      mode = 6;

    // minimal size fo rone Softclause
    //    if (j < _settings.minSize) mode = -1;

    if (mode == acceptedMode) {
      j = 0;

      _cascCandidates.back().allWeightsAreEqual =
          _cascCandidates.back().ggtTillPoint ==
          _originalSoftClauses[sortedSCs[i - 1]]->weight;
      _cascCandidates.back().Points = i;
      _cascCandidates.back().weightsTillPoint = _overallSoftWeights;
      // already next cascade!
      _cascCandidates.push_back(nextCandidate);
      _cascCandidates.back().ggtTillPoint =
          _originalSoftClauses[sortedSCs[i]]->weight;
    }

    if (_cascCandidates.back().ggtTillPoint != 1) {
      _cascCandidates.back().ggtTillPoint = GreatestCommonDivisor(
          static_cast<long long int>(_cascCandidates.back().ggtTillPoint),
          static_cast<long long int>(
              _originalSoftClauses[sortedSCs[i]]->weight));
    }

    _overallSoftWeights += _originalSoftClauses[sortedSCs[i]]->weight;
  }

  _cascCandidates.back().allWeightsAreEqual =
      (_cascCandidates.back().ggtTillPoint ==
       _originalSoftClauses[sortedSCs.back()]->weight);

  _cascCandidates.back().Points = _originalSoftClauses.size();
  _cascCandidates.back().weightsTillPoint = _overallSoftWeights;

  //    std::cout << "_cascCandidates.size() " << _cascCandidates.size() <<
  //    std::endl;

  //    if ((_cascCandidates[0].Points < _originalSoftClauses.size() / 15)
  //        && (_cascCandidates[0].weightsTillPoint >
  //        _cascCandidates[1].ggtTillPoint)) { _cascCandidates.clear();
  //        _cascCandidates.push_back(nextCandidate);
  //    }
  //    for (unsigned i = 1; i < _cascCandidates.size() - 1; ++i) {
  //        DumpCascCandidates();
  //        if ((_cascCandidates[i].weightsTillPoint > _cascCandidates[i +
  //        1].ggtTillPoint)) {
  //            &&(_cascCandidates[i].Points - _cascCandidates[i - 1].Points
  //               < _originalSoftClauses.size() / 15) _cascCandidates.clear();
  //            _cascCandidates.push_back(nextCandidate);
  //        }
  //        std::cout << "" << std::endl;
  //    }
}

void Pacose::RemoveCascCand(unsigned i) {
  _cascCandidates[i - 1].Points = _cascCandidates[i].Points;
  _cascCandidates[i - 1].weightsTillPoint = _cascCandidates[i].weightsTillPoint;
  _cascCandidates[i - 1].ggtTillPoint = GreatestCommonDivisor(
      _cascCandidates[i].ggtTillPoint, _cascCandidates[i - 1].ggtTillPoint);
  _cascCandidates[i - 1].allWeightsAreEqual = false;
  _cascCandidates.erase(_cascCandidates.begin() + i);
}

bool Pacose::CheckMinWeightDist(std::vector<unsigned> &sortedSCs,
                                unsigned firstPoint, unsigned long biggerThan) {
  if (_settings.verbosity > 3) {
    //    std::cout << "c first element: "
    //              << _originalSoftClauses[sortedSCs[firstPoint]]->weight
    //              << std::endl;
    std::cout << "c biggerThan: " << biggerThan << std::endl;
  }

  // at first easy check simple distance between the weights is already smaller
  // than the needed distance
  for (unsigned i = firstPoint; i < sortedSCs.size() - 1; i++) {
    // test if original weights have at least a distance of biggerThan
    if (_originalSoftClauses[sortedSCs[i + 1]]->weight !=
            _originalSoftClauses[sortedSCs[i]]->weight &&
        _originalSoftClauses[sortedSCs[i + 1]]->weight - biggerThan <
            _originalSoftClauses[sortedSCs[i]]->weight) {
      //      std::cout << "i: " << i << " "
      //                << _originalSoftClauses[sortedSCs[i + 1]]->weight << " "
      //                << _originalSoftClauses[sortedSCs[i]]->weight <<
      //                std::endl;
      std::cout << "c valid weight distance..: false" << std::endl;
      return false;
    }
  }

  // check all possible weight combinations
  // ATTENTION: NP hard problem!!!
  std::set<unsigned long> allWeightCombinations = {
      0, static_cast<unsigned long>(-1)};
  std::set<unsigned long>::iterator it, it2;
  std::pair<std::set<unsigned long>::iterator, bool> ret;
  for (unsigned i = firstPoint; i < sortedSCs.size(); i++) {
    unsigned long actualWeight = _originalSoftClauses[sortedSCs[i]]->weight;
    std::set<unsigned long> allCurrentCombinations = allWeightCombinations;
    for (it = ++allCurrentCombinations.begin();
         it != --allCurrentCombinations.end(); ++it) {
      ret = allWeightCombinations.insert(actualWeight + *it);
      if (ret.second == false) {
        // element was already in, continue with next element
        continue;
      }

      std::set<unsigned long>::iterator before = ret.first;
      std::set<unsigned long>::iterator after = ret.first;

      if (ret.second && (!(*ret.first - *--before > biggerThan) ||
                         !(*++after - *ret.first > biggerThan))) {
        if (_settings.verbosity > 0)
          std::cout << "c valid weight distance..: false2" << std::endl;
        return false;
      }
    }

    ret = allWeightCombinations.insert(actualWeight);
    if (allWeightCombinations.size() > 100000) {
      if (_settings.verbosity > 0) {
        std::cout << "c too many weight combinations: "
                  << allWeightCombinations.size() << std::endl;
        std::cout << "c valid weight distance..: unknown" << std::endl;
      }
      return false;
    }
    if (ret.second == false) {
      // element was already in, continue with next element
      //      std::cout << "Element was already in!" << std::endl;
      continue;
    }
    // Check weight distance for actual weight
    std::set<unsigned long>::iterator before = ret.first;
    std::set<unsigned long>::iterator after = ret.first;
    if (ret.second && (!(*ret.first - *--before > biggerThan) ||
                       !(*++after - *ret.first > biggerThan))) {
      if (_settings.verbosity > 0)
        std::cout << "c valid weight distance..: false3" << std::endl;
      return false;
    }
  }

  if (_settings.verbosity > 0)
    std::cout << "c valid weight distance..: true" << std::endl;
  return true;
}

unsigned long Pacose::DivideSCsIfPossible() {
  if (_settings.testIfDividable > 0)
    _settings.divideDGPW = USEONLYGCD;
  else
    _settings.divideDGPW = NODIVISION;
  if (_settings.divideDGPW == NODIVISION) {
    _sClauses.push_back(_originalSoftClauses);
    partitionInformation nextCandidate;
    _cascCandidates.push_back(nextCandidate);
    return _sClauses.size();
  }

  _settings.minSize =
      std::floor(_settings.minSize * (0.01 * _originalSoftClauses.size()));
  if (_settings.minSize != 0)
    std::cout << "c min size of DGPW.......: " << _settings.minSize
              << std::endl;

  // get indices of sorted Bucket entries
  // generate Indice Vector to sort that vector!
  std::vector<unsigned> sortedSCIndices(_originalSoftClauses.size());
  std::size_t n(0);
  std::generate(std::begin(sortedSCIndices),
                std::begin(sortedSCIndices) +
                    static_cast<unsigned>(_originalSoftClauses.size()),
                [&] { return n++; });

  // stable sort - not changing order of SC's important for some of the
  // instances! especially spot5!
  std::stable_sort(std::begin(sortedSCIndices), std::end(sortedSCIndices),
                   [&](std::size_t i1, std::size_t i2) {
                     return (_originalSoftClauses[i2]->weight >
                             _originalSoftClauses[i1]->weight);
                   });

  if (_settings.divisionMode == -1) {
    for (int mode = 0; mode <= 6; mode++) {
      DivideSCs(sortedSCIndices, mode);
      if (_cascCandidates.size() > 1 ||
          (mode == 0 && _settings.divideDGPW == USEONLYGCD)) {
        if (_settings.verbosity > 0)
          std::cout << "c Accepted Mode..........: " << mode << std::endl;
        break;
      }
    }
  } else if (_settings.divisionMode >= 0 && _settings.divisionMode <= 6) {
    DivideSCs(sortedSCIndices, _settings.divisionMode);
  } else {
    std::cout << "c ERROR: no valid division!!" << std::endl;
    exit(0);
  }

  if (_cascCandidates.size() > 1) {
    //        std::cout << "c partitionFactor........: " <<
    //        _settings.partitionFactor << std::endl;

    // test if gcd of the bigger potential cascades i greater than the sum of
    // weights of the lower cascades
    for (unsigned i = static_cast<unsigned>(_cascCandidates.size() - 1); i > 0;
         --i) {
      unsigned long minWeightDistance = _cascCandidates[i - 1].weightsTillPoint;
      //      std::cout << "c wTP[" << i - 1 << "]!"
      //                << _cascCandidates[i - 1].weightsTillPoint << std::endl;

      //      std::cout << "ggtTillPoint: " << _cascCandidates[i].ggtTillPoint
      //                << std::endl;
      //      std::cout << "minSize: " << _settings.minSize << std::endl;
      //      std::cout << "Points: "
      //                << _cascCandidates[i].Points - _cascCandidates[i -
      //                1].Points
      //                << std::endl;

      if (_settings.testIfDividable == 2) {
        // test for bigger equal
        minWeightDistance = _cascCandidates[i - 1].weightsTillPoint - 1;
      }

      bool isBigger = false;
      if ((_cascCandidates[i].Points - _cascCandidates[i - 1].Points >=
           _settings.minSize) &&
          (i != 1 || _cascCandidates[0].Points >= _settings.minSize))
        isBigger = true;

      //      std::cout << "i " << i << " isbigger: " << isBigger << std::endl;
      //      std::cout << "i " << i << " minWeightDistance: " <<
      //      minWeightDistance
      //                << std::endl;
      //      std::cout << "i " << i << " _cascCandidates[i].ggtTillPoint: "
      //                << _cascCandidates[i].ggtTillPoint << std::endl;

      if ((_cascCandidates[i].ggtTillPoint > minWeightDistance) && isBigger) {
        if (_settings.verbosity > 0)
          std::cout << "c VALID subCascade[" << i << "]!" << std::endl;
      } else {
        if (!CheckMinWeightDist(sortedSCIndices, _cascCandidates[i - 1].Points,
                                minWeightDistance) ||
            !isBigger) {
          if (_settings.verbosity > 0)
            std::cout << "c INVALID subCascade[" << i << "]!" << std::endl;
          RemoveCascCand(i);
        }
      }
    }

    if (_cascCandidates.size() > 1) {
      if (_settings.verbosity > 0)
        std::cout << "c number of subproblems..: " << _cascCandidates.size()
                  << std::endl;
    }

    _settings.formulaIsDivided = true;
    unsigned nextPartition = 0;
    std::vector<SoftClause *> SC;
    _sClauses.push_back(SC);
    for (unsigned i = 0; i < sortedSCIndices.size(); ++i) {
      if (i == _cascCandidates[nextPartition].Points) {
        nextPartition++;
        _sClauses.push_back(SC);
      }
      _sClauses.back().push_back(_originalSoftClauses[sortedSCIndices[i]]);
    }
  }

  if (_cascCandidates.size() <= 1) {
    // COPY CONSTRUCTOR??
    if (_sClauses.size() == 0) {
      std::vector<SoftClause *> SC;
      _sClauses.push_back(SC);
    }
    _sClauses.back() = _originalSoftClauses;
  } else {
    DumpCascCandidates();
  }

  if (_settings.verbosity > 3) {
    std::cout << std::endl;
    for (auto i : _sClauses) {
      std::cout << "c ";
      for (auto j : i) {
        std::cout << j->weight << " ";
      }
      std::cout << std::endl << std::endl;
    }
  }

  if (_settings.divCheck) exit(0);

  return _sClauses.size();
}

void Pacose::Preprocess() {
  //    GreedyMaximizeInitialSATWeight();

  // not compatible with dividing SC's into subformulas!
  if (_settings.divideDGPW == NODIVISION) {
    AnalyzeSCsAndConvertIfPossible();
  }

  CalcGCDAndDivideIfPossible();
}

void Pacose::AnalyzeSCsAndConvertIfPossible() {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;
  if (!_settings.GetAnalyzeFormula()) return;

  if (_actualSoftClauses->size() == 0) {
    _settings.SetFormulaType(FormulaType::SAT);
    std::cout << "c is pure SAT formula......: true" << std::endl;
    return;
  }

  long long int minWeight = (*_actualSoftClauses)[0]->weight;
  long long int maxWeight = (*_actualSoftClauses)[0]->weight;
  long long int prevMaxWeight = (*_actualSoftClauses)[0]->weight;
  bool hasMoreThanTwoWeights = false;
  bool recalc = false;
  long long int sumOfMinWeights = 0;

  for (int i = 0; i < _actualSoftClauses->size(); i++) {
    //        std::cout << (*_actualSoftClauses)[i]->weight << " ";

    if ((*_actualSoftClauses)[i]->weight < minWeight) {
      minWeight = (*_actualSoftClauses)[i]->weight;
      recalc = true;
    } else if ((*_actualSoftClauses)[i]->weight > maxWeight) {
      prevMaxWeight = maxWeight;
      maxWeight = (*_actualSoftClauses)[i]->weight;
    }

    if (minWeight != maxWeight && minWeight != prevMaxWeight &&
        maxWeight != prevMaxWeight) {
      hasMoreThanTwoWeights = true;
    }

    if ((*_actualSoftClauses)[i]->weight == minWeight) {
      if (recalc) {
        sumOfMinWeights = 0;
        recalc = false;
      }
      sumOfMinWeights += minWeight;
    }
  }

  //  recognize maxSAT instances
  if (minWeight == maxWeight) {
    if (!_hasHardClauses) {
      for (int i = static_cast<int>(_actualSoftClauses->size() - 1); i < 0;
           i--) {
        if ((*_actualSoftClauses)[static_cast<size_t>(i)]->clause.size() == 1) {
          // wrongly interpreted as unit soft clause, has to be added as unit
          // hard clause
          _satSolver->NewClause();
          unsigned int ulit =
              (*_actualSoftClauses)[static_cast<size_t>(i)]->clause[0];
          _satSolver->AddLiteral(&ulit);
          _satSolver->CommitClause();
          _satSolver->ResetClause();

          //                    _satSolver->AddClause(_softClauses[i]->clause);
        } else {
          // deactivate trigger literal by adding it as negated unit clause.
          _satSolver->NewClause();
          unsigned int ulit =
              (*_actualSoftClauses)[static_cast<size_t>(i)]->relaxationLit ^ 1;
          _satSolver->AddLiteral(&ulit);
          _satSolver->CommitClause();
          _satSolver->ResetClause();
        }
      }
      _hasHardClauses = true;
      _settings.SetFormulaType(FormulaType::SAT);
      std::cout << "c converted to SAT.......: true" << std::endl;
      // CAN BE CONTINUED WITH SAT SOLVER CALL, NOT YET IMPLEMENTED!
    } else {
      _settings.SetFormulaType(FormulaType::MAXSAT);
      std::cout << "c is MaxSAT formula......: true" << std::endl;
    }
  }

  //    // has exactly two weights and no hard clauses - it can be converted to
  //    pure unweighted maxSAT std::cout << "!hasMoreThanTwoWeights: " <<
  //    !hasMoreThanTwoWeights
  //              << "  !_hasHardClauses: " << _hasHardClauses
  //              << "  minWeight: " << minWeight
  //              << "  maxWeight: " << maxWeight
  //              << "  sumOfMinWeights: " << sumOfMinWeights
  //              << "  maxWeight: " << maxWeight
  //              << std::endl;

  if (!hasMoreThanTwoWeights && !_hasHardClauses && minWeight != maxWeight &&
      sumOfMinWeights < maxWeight) {
    std::vector<SoftClause *> newSoftClauses;
    //        std::vector< SoftClause* > newSoftClauses;
    //        std::vector<std::vector<u_int32_t>> newLiterals;

    //        std::cout << "maxWeight " << maxWeight << std::endl;
    for (int ind = 0; ind < static_cast<int>(_actualSoftClauses->size());
         ++ind) {
      if ((*_actualSoftClauses)[static_cast<size_t>(ind)]->weight !=
          maxWeight) {
        newSoftClauses.push_back(
            (*_actualSoftClauses)[static_cast<size_t>(ind)]);
      }
      //            else if
      //            ((*_actualSoftClauses)[static_cast<size_t>(ind)]->clause.size()
      //            == 1)
      //            {
      //                std::cout << "in here" << std::endl;
      //                // wrongly interpreted as unit soft clause, has to be
      //                added as unit hard clause _satSolver->NewClause();
      //                unsigned int ulit =
      //                (*_actualSoftClauses)[static_cast<size_t>(ind)]->clause[0];
      //                _satSolver->AddLiteral(&ulit);
      //                _satSolver->CommitClause();
      //                _satSolver->ResetClause();
      //            }
      else {
        //                std::cout << (*_actualSoftClauses)[ind]->weight << " "
        //                << (*_actualSoftClauses)[ind]->relaxationLit << " ";
        // deactivate trigger literal by adding it as negated unit clause.
        _satSolver->NewClause();
        unsigned int ulit =
            (*_actualSoftClauses)[static_cast<size_t>(ind)]->relaxationLit ^ 1;
        //                std::cout << ulit << ", " << std::endl;
        _satSolver->AddLiteral(&ulit);
        _satSolver->CommitClause();
        _satSolver->ResetClause();
      }
    }
    _hasHardClauses = true;
    //        _sClauses[1] = newSoftClauses;
    (*_actualSoftClauses) = newSoftClauses;

    _settings.SetFormulaType(FormulaType::MAXSAT);
    std::cout << "c converted to MaxSAT....: true" << std::endl;
    std::cout << "c Remaining SoftClauses..: " << _actualSoftClauses->size()
              << std::endl;

    //        for( auto SC : *_actualSoftClauses ) {
    //            std::cout << SC->weight << " ";
    //        }
    //        std::cout << std::endl;

    _sumOfSoftWeights = sumOfMinWeights;
  }
}

long long int Pacose::GreatestCommonDivisor(long long int a, long long int b) {
  if (a == b) {
    return a;
  }
  //    std::cout << " a: " << a << "  b: " << b << std::endl;
  long long int temp;
  while (b > 0) {
    temp = b;
    b = a % b;
    a = temp;
  }
  //    std::cout << "xa: " << a << "  b: " << b << std::endl;
  return a;
}

void Pacose::CalcGCDAndDivideIfPossible() {
  if (!_settings.GetAnalyzeFormula()) return;

  _GCD = _minWeight;
  //    std::cout << "minweight: " << _minWeight << std::endl;

  for (int ind = 0; ind < _actualSoftClauses->size(); ++ind) {
    //        std::cout <<"SCI.weight: " << _softClauses[ind]->weight <<
    //        std::endl;
    if (_GCD == 1) break;
    _GCD = GreatestCommonDivisor(_GCD, (*_actualSoftClauses)[ind]->weight);
  }

  if (_GCD > 1) {
    //        std::cout << "c greatest common divisor: " << _GGT << std::endl;
    for (int ind = 0; ind < _actualSoftClauses->size(); ++ind) {
      (*_actualSoftClauses)[ind]->weight /= _GCD;
    }
  }
}

void Pacose::SetSumOfSoftWeights(unsigned long softWeights) {
  _sumOfSoftWeights = softWeights;
}
