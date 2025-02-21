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

#ifndef SETTINGS_H
#define SETTINGS_H

#include <stdint.h>
#include <iostream>
#include <string>

enum EncodingType {
  WARNERS,   // adder warners 1996, [Warners 1998?]
  BAILLEUX,  // totalizer [Bailleux & Boufkhad 2003]
  ASIN,  // [Asin et. al 2011] Robert Asin, Robert Nieuwenhuis, Albert Oliveras,
         // Enric Rodriguez-Carbonell "Cardinality Networks: a theoretical and
         // empirical study"
  OGAWA,       // modulo Totalizer [Ogawa et. al 2013]
  BAILLEUXW2,  // BailleuxW2
  WMTO,        // Weighted MaxSAT Totalizer
  MRWTO,       // Mixed Radix Weighted Totalizer
  MRWTO2,      // Mixed Radix Weighted Totalizer
  DGPW18,      // Dynamic Global Polynomial Watchdog [Paxian & Reimer 2018]
  HEURISTICQMAXSAT,  // auto heuristic [Koshi, 2014]
                     // selects between "warn, bail, ogaw"
  HEURISTIC20  // HEURISTIC 20 Competition, choosing between warners and dgpw
};

enum SolverType { SOLVER, SIMPSOLVER, PARALLELSOLVER, MAXSATSOLVER };

enum FormulaType { SAT, MAXSAT, WEIGHTEDMAXSAT };

enum EncodeStrategy { ENCODEALL, ENCODEONLYIFNEEDED };

enum PartitionStrategy {
  NOPARTITION,
  GROUPBYWEIGHTADDATLAST,
  GROUPBYWEIGHT,
  GROUPBYBIGGESTREPEATINGENTRY,
  GROUPBYCOLUMNS
};

enum MultipleCascadeSolvingState {
  SINGLECASCADE,
  MAINCASCADE,
  TARECASCADE,
  TARETARES
};

enum MultipleCascadeDivideStrategy {
  SOLVEINNORMALCASCADEMODE,
  SORTEDNUMBEROFSOFTCLAUSES,
  RANDOMNUMBEROFSOFTCLAUSES,
  SOFTCLAUSESINORDER,
  SORTEDGOODDIVIDEPOINTS
};

enum InterimResult {
  NOINTERIMRESULT,
  CUTATTOP,
  CUTATBOTTOM,
  CUTBOTH,
};

enum StructureInfo {
  ISSAT,
  CONVERTTOMAXSAT,
  ISMAXSAT,
  DIVWEIGHTSBYDIVISOR,
  ISWEIGHTEDMAXSAT,
};

enum SorterType {
  BITONIC,
  ODDEVEN,
  TOTALIZER,
};

enum DivideDGPWStrategy {
  NODIVISION,
  USEONLYGCD,
  DIVIDEALL,
  DIVIDEALLSOLVEONLYWATCHDOGS,
  DIVIDEALLSOLVEGREEDY,
};

struct Settings {
 public:
  Settings()
      : _encoding(HEURISTIC20),
        _solverType(SOLVER),
        _compression(-1),
        _printModel(true),
        _analyzeFormula(true) {
    ResetDGPW();
    ResetCore();
    // std::cout << "CARD=" << card << std::endl;
    verbosity = 0;
  }

  /**
   * @brief ResetCore
   *          sets all settings to standard values.
   */
  void ResetCore() {
    _encoding = HEURISTIC20;
    _solverType = SOLVER;
    _formulaType = WEIGHTEDMAXSAT;
    _compression = -1;
    _printModel = true;
    _analyzeFormula = true;
  }

  struct CurrentCascade {
    bool _onlyWithAssumptions = false;
    bool _solveTares = true;
    unsigned iteration = 0;
  } currentCascade;

  void ResetDGPW(void) {
    reconf = 99;
    // MaxSAT Settings
    _encoding = HEURISTIC20;
    _compression = -1;
    greedyPPFixSCs = -1;
    simplify = false;
    networkType = TOTALIZER;
    preciseTarget = false;
    targetOpt = -1;
    greedyPrepro = 1;
    divCheck = false;
    greedyMinSizeOfSet = 1;
    useGreedyPreInBetween = false;
    atLeastnEqualWeights = 1;
    greedyPPTimeLimit = 600;
    greedyPPSATPercentage = 65;
    greedyPPUNSATPercentage = 65;
    greedyPPPropagationPerSecond = 4000000;
    greedyPPTimeoutFactor = 18;
    greedyPPSSCSwitchFactor = 75;

    // Weighted MaxSAT Settings
    encode01 = false;
    lastPos1 = true;
    base = 2;
    groupHeuristic = 1;
    percentOff = 100;
    percentOffReinsert = true;
    equalWeight = 0;
    partitionStrategy = NOPARTITION;
    solveAtFirst = true;
    encodeStrategy = ENCODEONLYIFNEEDED;
    createGraphFile = "";
    onlyByTares = false;
    solveAsWeighted = true;
    formulaIsDivided = false;
    partitionFactor = 1;
    testIfDividable = 2;
    checkIfSolutionIsUnique = false;

    // multiple cascade mode
    mcDivideStrategy = SOLVEINNORMALCASCADEMODE;
    cascadeDivider = 0;
    maxBucketSize = 0;
    nOfCasc = 0;
    tareCascadeOnlyByTares = false;
    sepHiWeight = false;
    weightPlusOne = true;
    minSize = 0;

    // interim results
    interimResult = NOINTERIMRESULT;

    featureTest = 0;
    analyze = true;

    // Paper Options
    adderCaching = true;
    coneOfInfluence = true;
    exactBounding = true;
    plainVariant = false;
    dGPW = true;

    SetPaperOptions();
  }

  void SetPaperOptions() {
    if (plainVariant || adderCaching || coneOfInfluence || exactBounding) {
      encodeStrategy = ENCODEALL;
      partitionStrategy = NOPARTITION;
      lastPos1 = false;
      weightPlusOne = false;
      encode01 = false;
      solveAtFirst = true;
      solveAsWeighted = true;
      analyze = true;
      if (adderCaching && coneOfInfluence && exactBounding) {
        dGPW = true;
        plainVariant = false;
      }
    }
    if (adderCaching) {
      partitionStrategy = GROUPBYWEIGHT;
      //              partitionStrategy = GROUPBYBIGGESTREPEATINGENTRY;
      //              partitionStrategy = GROUPBYCOLUMNS;
      //              createGraphFile = "graph";
    }
    if (coneOfInfluence) {
      encodeStrategy = ENCODEONLYIFNEEDED;
      lastPos1 = true;
    }
    if (exactBounding) {
      weightPlusOne = true;
      // had wrong value in competition version
      // slightly different from what we explained in the paper
      //              onlyByTares = true;
    }
  }

  void SetSolverType(SolverType solverType) { _solverType = solverType; }

  void SetEncoding(EncodingType encoding) { _encoding = encoding; }

  void SetFormulaType(FormulaType formulaType) { _formulaType = formulaType; }

  void SetCompression(int compression) { _compression = compression; }

  void SetPrintModel(bool printModel) { _printModel = printModel; }

  bool SetAnalyzeFormula() { return _analyzeFormula; }

  SolverType GetSolverType() { return _solverType; }

  EncodingType GetEncoding() { return _encoding; }

  FormulaType GetFormulaType() { return _formulaType; }

  int GetCompression() { return _compression; }

  bool GetPrintModel() { return _printModel; }

  bool GetAnalyzeFormula() { return _analyzeFormula; }

  void Print() const {
    if (verbosity == 0) return;
    std::cout << "c Plain Variant..........: "
              << (plainVariant ? "true" : "false") << std::endl;
    std::cout << "c Adder Caching..........: "
              << (adderCaching ? "true" : "false") << std::endl;
    std::cout << "c Cone Of Influence......: "
              << (coneOfInfluence ? "true" : "false") << std::endl;
    std::cout << "c Exact Bounding.........: "
              << (exactBounding ? "true" : "false") << std::endl;
    std::cout << "c DGPW combined..........: " << (dGPW ? "true" : "false")
              << std::endl;
    std::cout << "c verbosity..............: " << verbosity << std::endl;
    std::cout << "c ------------------------" << std::endl;
    std::cout << "c MaxSAT options:" << std::endl;
    std::cout << "c weight plus 1..........: " << weightPlusOne << std::endl;
    std::cout << "c network type...........: ";
    switch (networkType) {
      case BITONIC:
        std::cout << "Bitonic sorter" << std::endl;
        break;
      case ODDEVEN:
        std::cout << "Odd-Even sorter" << std::endl;
        break;
      case TOTALIZER:
        std::cout << "Totalizer" << std::endl;
        break;
    }
    //          std::cout << "c decision strategies....: ";
    std::cout << "c encode 01..............: " << (encode01 ? "true" : "false")
              << std::endl;
    //          std::cout << "c sort soft clauses......: ";
    std::cout << "c base...................: " << base << std::endl;
    std::cout << "c analyze................: " << (analyze ? "true" : "false")
              << std::endl;
    std::cout << "c partitionStrategy......: " << partitionStrategy
              << std::endl;
    std::cout << "c   heuristic............: " << groupHeuristic << std::endl;

    std::cout << "c defined target optimum.: " << targetOpt << std::endl;
  }

  //  private:
  EncodingType _encoding;
  SolverType _solverType;
  FormulaType _formulaType;
  // compression of encoding due to comp in QMaxSAT;
  // 0 highest compression
  int _compression;
  bool _printModel;
  bool _analyzeFormula;

  // General Settings
  uint32_t verbosity;
  double cpuLimit;
  uint32_t memLimit;
  uint32_t maxWidth;
  SorterType networkType;
  bool preciseTarget;
  int32_t targetOpt;

  // Weighted MaxSAT Settings
  bool analyze;
  bool encode01;
  bool lastPos1;
  uint32_t base;
  uint32_t groupHeuristic;
  uint32_t percentOff;
  bool percentOffReinsert;
  uint32_t equalWeight;
  PartitionStrategy partitionStrategy;
  bool solveAtFirst;
  EncodeStrategy encodeStrategy;
  std::string createGraphFile;
  bool onlyByTares;
  bool solveAsWeighted;
  bool formulaIsDivided;

  // Paper Options
  bool adderCaching;
  bool coneOfInfluence;
  bool exactBounding;
  bool plainVariant;
  bool dGPW;

  // multiple cascade mode
  MultipleCascadeDivideStrategy mcDivideStrategy;
  uint32_t cascadeDivider;
  uint32_t maxBucketSize;
  uint32_t nOfCasc;
  bool tareCascadeOnlyByTares;
  bool sepHiWeight;
  bool weightPlusOne;

  // interim results
  InterimResult interimResult;
  DivideDGPWStrategy divideDGPW;

  unsigned minSize;
  unsigned partitionFactor;
  //  bool useGreedyPrepro;
  int greedyPrepro;
  int greedyPPTimeLimit;
  int greedyPPPropagationPerSecond;
  int greedyPPTimeoutFactor;
  int greedyPPSSCSwitchFactor;
  int greedyPPFixSCs;
  int greedyPPSATPercentage;
  int greedyPPUNSATPercentage;
  int greedyMinSizeOfSet;

  bool useGreedyPreInBetween;
  unsigned atLeastnEqualWeights;
  int divisionMode;
  int testIfDividable;
  bool divCheck;
  bool checkIfSolutionIsUnique;
  // CMS configuration Option
  unsigned reconf;
  bool simplify;

  int featureTest;
};

#endif  // SETTINGS_H
