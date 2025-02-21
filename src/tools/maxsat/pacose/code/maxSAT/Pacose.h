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

#ifndef PACOSE_H
#define PACOSE_H

#include "dgpw.h"
//#include "Encodings.h"
#include <iostream>
#include <vector>
#include "SATSolverProxy.h"
#include "Settings.h"

// class OutOfTimeException{};
// class DGPW;

struct SoftClause;
// class SATSolverProxy;
// class SATSolverType;
class Encodings;
// class DGPW;

// Pacose is an extension to the glucose solver
// containing
/**
 * @brief The Pacose class
 *      1.  (Partial) (Weighted) MaxSAT Solver
 *          Initially derived from QMaxSAT (2017 Competition Version)
 */
class Pacose {
 public:
  Pacose();

  ~Pacose();

  bool parseWcnfFile(std::string wcnfFile);
  bool parseWcnfFileOld(std::string wcnfFile);
  unsigned SolveProcedure(long long int nbVars = 0);
  unsigned SolveQMax(std::vector<SoftClause *> *tmpSoftClauses = nullptr,
                     EncodingType *encodingType = nullptr);
  void Preprocess();

  // Variables
  Settings _settings;
  //    Encodings _encodings;

  SATSolverProxy *_satSolver;
  unsigned _nbVars;
  unsigned _nbClauses;
  uint64_t _top;  // hard clause weight
  std::vector<SoftClause *> _originalSoftClauses;
  std::vector<SoftClause *> *_actualSoftClauses;
  std::vector<std::vector<SoftClause *>> _sClauses;

  uint64_t CalculateSATWeight();
  long long int GetSATWeight() { return _satWeight; }

  void InitSatSolver(SATSolverType solverType);
  void AddSoftClause(std::vector<uint32_t> &clause, uint64_t weight = 1);

  bool _hasHardClauses;
  void SetSumOfSoftWeights(unsigned long softWeights);

 private:
  struct partitionInformation {
    DGPW::DGPW *dgpw;
    unsigned Points;
    unsigned long weightsTillPoint;
    unsigned long ggtTillPoint;
    bool allWeightsAreEqual;
  };

  std::vector<partitionInformation> _cascCandidates;
  EncodingType _encoding;

  // Variables
  int _cpuLimit;
  int _memLimit;
  int _nbOfOrigVars;
  unsigned long _sumOfSoftWeights;
  unsigned long _overallSoftWeights;
  // fulfilled softclauses
  long long int _satWeight;
  // o-value
  long long int _unSatWeight;
  long long int _localUnSatWeight;
  long long int _minWeight;
  long long int _maxWeight;
  long long int _GCD;

  Encodings *_encodings;

  // relaxation literals or unit clause values!
  std::vector<unsigned> _blockings;
  unsigned _negRelaxLit;
  // weights of relaxation literals
  std::vector<long long int> _weights;

  unsigned SignedToUnsignedLit(int literal);
  void CalcGCDAndDivideIfPossible();
  long long GreatestCommonDivisor(long long a, long long b);
  void AnalyzeSCsAndConvertIfPossible();
  unsigned long DivideSCsIfPossible();
  void genCardinals(long long tmpUnSATWeight, long long &divisor,
                    std::vector<unsigned> &lits,
                    std::vector<unsigned> &linkingVar,
                    std::vector<long long> &linkingWeight,
                    std::vector<long long> &divisors,
                    std::vector<std::vector<unsigned>> &linkingVars,
                    std::vector<std::vector<long long>> &linkingWeights,
                    int compression);
  void HeuristicQMaxSAT(long long sum, long long k);
  void wbSortAndFilter(long long UnSATWeight);
  void DumpCascCandidates();
  void AddSoftClauseTo(std::vector<SoftClause *> *softClauseVector,
                       std::vector<unsigned> &clause, uint64_t weight);
  void GreedyMaximizeInitialSATWeight(unsigned maxTime = 10,
                                      unsigned maxSolves = 1000);
  std::vector<unsigned> GetBestSCAssignment();
  void DivideSCs(std::vector<unsigned> &sortedSCs, int acceptedMode = 0);
  void CalculateOverallTimes();
  void RemoveCascCand(unsigned i);
  bool CheckMinWeightDist(std::vector<unsigned> &sortedSCs, unsigned firstPoint,
                          unsigned long biggerThan);

  bool AddEncoding(std::vector<SoftClause *> *tmpSoftClauses = nullptr,
                   EncodingType *encodingType = nullptr);

  /**
   * @brief SolveMaxSAT straight MaxSAT solving of SoftClauseVector plus
   * encoding type without any preprocessing
   * @param tmpSoftClauses
   * @param encodingType
   */
  unsigned SolveMaxSAT(std::vector<SoftClause *> *tmpSoftClauses = nullptr,
                       EncodingType *encodingType = nullptr);
  uint64_t CalculateLocalSATWeight(
      std::vector<SoftClause *> *tmpSatClauses = nullptr);
  void PrintResult();
};

#endif  // PACOSE_H
