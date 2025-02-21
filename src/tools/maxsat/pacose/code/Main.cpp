/********************************************************************************************
SATSolverProxy.cpp -- Copyright (c) 2020, Tobias Paxian
    Parts are taken from the QMaxSAT 2017 MaxSAT evaluation version.

    Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal in
    the Software without restriction, including without limitation the rights to
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

#include <iostream>
#include "CLI11.hpp"  // Parser
#include "MaxSATDimacs.h"
#include "Pacose.h"  // MaxSAT solver
#include "Settings.h"
//#include "externalPrepro/externalPrepro.h"

SATSolverType returnSolverType(std::string solver) {
  if (solver == "antom") {
    return SATSolverType::ANTOM;
  }
  if (solver == "glucose3") {
    return SATSolverType::GLUCOSE3;
  }
  if (solver == "glucose421") {
    return SATSolverType::GLUCOSE421;
  }
  if (solver == "cryptominisat") {
    return SATSolverType::CRYPTOMINISAT;
  }
  if (solver == "cadical") {
    return SATSolverType::CADICAL;
  }
  if (solver == "lingeling") {
    return SATSolverType::LINGELING;
  }

  return SATSolverType::MAPLEGLUCOSE;
}

//=================================================================================================
int main(int argc, char **argv) {
  // Option Parsing with CLI11
  CLI::App app{"\nP a c o s e, an iterative MaxSAT Solver\n"};

  // Define options
  std::string encoding = "auto";
  std::string solver = "glucose421";
  std::string maxCnfFile;

  Pacose *pacose = new Pacose();
  Settings *settings = &pacose->_settings;
  std::string encode01Parser = "false";
  int sortSoftClauses = 0;
  int partitionStrategyParser;
  int encodeStrategyParser = 1;    // ??
  int mcDivideStrategyParser = 0;  // ??
  int interimResultParser = 0;     // ??
  int divideDGPWParser = 0;
  // uint32_t verbosity;
  // double cpuLimit;
  // uint32_t memLimit;

  app.add_set("-e, --encoding", encoding,
              {"dgpw", "warn", "bail", "asin", "ogaw", "bailw2", "wmt", "mrwto",
               "mrwto2", "qauto", "auto"},
              " \t Different MaxSAT encoding to choose from.", true);

  app.add_option("-p, --printModel", settings->_printModel,
                 "\t Print solution model. (default: true)", true);

  app.add_set(
      "--setQMSATCompression", settings->_compression,
      {-1, 0, 1, 2, 10, 11, 20, 21, 30, 99},
      "\t Wcnf is dividable (division mode 1!)\n"
      "\t\t\t -1: heuristic to choose between 1 and 99 -- only WARNERS\n"
      "\t\t\t 0: highest compression\n"
      "\t\t\t 1,2,10,11: other compression modes, not available for "
      "all encodings\n");

  app.add_option("-f, --file, 1", maxCnfFile, "\t The used max CNF file")
      ->required()
      ->check(CLI::ExistingFile);
  app.add_set("-s, --solver", solver,
              {"antom", "glucose3", "glucose421", "cryptominisat", "cadical",
               "lingeling", "mapleglucose"},
              "The chosen solver.", true);

  app.add_flag(
      "--simplify", settings->simplify,
      "\t simplifies instance (at the moment only glucose4).");  // if argument
                                                                 // adderCaching
                                                                 // exists then
  // true, else false;
  app.add_set("--reconf", settings->reconf, {0, 3, 4, 6, 7, 12, 13, 14, 15, 16},
              "\t sort soft clause parts for partialMode (default: 0) \n");

  app.add_option("-v, --verbose", settings->verbosity,
                 "\t Increase verbosity \n "
                 "\t \t \t \t 0 reset verbosity",
                 true);

  app.add_option("-c, --cpuLimit", settings->cpuLimit,
                 "\t CPU limit in seconds \n "
                 "\t \t \t \t 0.0: no CPU limit at all (default)");
  app.add_option("-m, --memLimit", settings->memLimit,
                 "\t Memory limit in MB \n "
                 "\t \t \t \t  0: no memory limit at all (default)");

  //  if (encoding == "dgpw") {
  app.add_set("--encode01", encode01Parser, {"true", "false", "lastPos1"},
              "encode complete 01-comparators or just half (default:true)");
  app.add_set("--sortsoft", sortSoftClauses, {0, 1, 2, 3},
              "\t sort soft clause parts for partialMode (default: 0) \n"
              "\t \t \t \t 0: no sorting \n"
              "\t \t \t \t 1: sort soft clauses due to the larger number of "
              "conflicting soft "
              "clauses \n"
              "\t \t \t \t 2: sort soft clauses due to the smaller number of "
              "conflicting "
              "soft clauses \n"
              "\t \t \t \t 3: random sort ");
  app.add_set("--partitionStrategy", partitionStrategyParser, {0, 1, 2, 3, 4},
              "\t sets the strategy how to fill the buckets (default: 3) \n"
              "\t \t \t \t 0: Standard, fill buckers with ungrouped values \n"
              "\t \t \t \t 1: Grouping weights - depth problematic \n"
              "\t \t \t \t 2: Grouping weights - balanced depth\n"
              "\t \t \t \t 3: Group weights then by biggest repeating entry "
              "due to a heuristic");
  // brauche ich diese if?
  // if (settings->partitionStrategy == GROUPBYBIGGESTREPEATINGENTRY) {
  // if this case a heuristic needs to be chosen,
  // max difference in sizes to merge SoftClauseNodes
  app.add_set(
      "--heuristic", settings->groupHeuristic, {0, 1, 2, 3, 4, 5, 6, 7, 8},
      " \t 0: Standard combination of other heuristics \n"
      "\t \t \t \t 1: number of reduced (ternary clauses * 2 + binary "
      "clauses)\n"
      "\t \t \t \t 2: sum of bucket sizes the SC occurs\n"
      "\t \t \t \t 3: clisest size in percentage\n"
      "\t \t \t \t 4: the one with least occurences in other buckets\n"
      "\t \t \t \t 5: how many merges are possible after this merge\n"
      "\t \t \t \t 6: greatest depth of submerges (till depth 3, then adds "
      "possible merges)\n"
      "\t \t \t \t ?: number of reduceed ternary clauses claculating all "
      "possible subMerges");
  app.add_option(
      "--percentOff", settings->percentOff,
      "\t The max percentage two weights differ with partitionStrategy=3 "
      "(default: 100), [0, 100]",
      true);
  app.add_option("--percentOffReinsert", settings->percentOffReinsert,
                 "\t reinsert if size is reached again (default: true)",
                 true);  // kommt das auch nur vor wenn partitionstrategyParser
                         // == 3 ?, anstatt option geht auch flag.
  app.add_option("--equalWeight", settings->equalWeight,
                 "\t merges equal weights every value rounds (default: 0 <- "
                 "not checking)",
                 true);
  app.add_option("--analyze", settings->analyze,
                 " \t converts formulas to ((Partial) Max) SAT if possible "
                 "(default: true)\n"
                 "\t \t \t \t recognizes if formula has a common divisor - "
                 "and divides all "
                 "weights by it",
                 true);  // auch noch im if? und folgende??
  app.add_option(
      "--solveAtFirst", settings->solveAtFirst,
      "\t starts with a solver call at the beginning. (default: true)", true);
  app.add_set("--encodeStrategy", encodeStrategyParser, {0, 1},
              "\t sets strategy of how to encode the buckets (default: 1)\n"
              "\t \t \t \t 0: Encode all.\n"
              "\t \t \t \t 1: Encode only if needed");
  app.add_option("--createGraph", settings->createGraphFile,
                 "\t MISSING COMMENT");  // was ist mit zeile 223 - 231
  app.add_option("--cascDiv", settings->cascadeDivider,
                 "\t divider of the SCs or max bucket size\n"
                 "\t \t \t \t if also sets the limit for the max bucket "
                 "size, if cascade is reunioned\n"
                 "\t \t \t \t if not defined extra",
                 true);
  app.add_option("--maxBucketSize", settings->maxBucketSize,
                 "\t Sets the max bucket size, if cascade is reunioned", true);
  app.add_option("--nOfCasc", settings->nOfCasc,
                 "\t Sets the number of cascades in multiple mode, Stronger "
                 "thatn cascDiv.",
                 true);
  app.add_option("--onlyByTares", settings->onlyByTares,
                 "\t solve whole cascade only by tares, adding therefore as "
                 "many buckets as necessary.\n"
                 "\t \t \t \t Works only with encodeStrategy = 1 (default)",
                 true);  // flag? - default??? und condition habe ich nicht
                         // gefundne im parser? anderer Code?
  app.add_option("--tcOnlyByTares", settings->tareCascadeOnlyByTares, " ",
                 true);  // flag?
  app.add_set("--mcDivideStrategy", mcDivideStrategyParser, {0, 1, 2, 3, 4},
              "\t sets the strategy for dividing the Softclauses into "
              "multiple cascades.\n"
              "\t \t \t \t 0: solves the cascade in normal mode. (default) \n"
              "\t \t \t \t 1: sorts the SoftClauses, max SCs in one sub "
              "cascade <= cascadeDiv. \n"
              "\t \t \t \t 2: SCs are randomly distributed, max SCs in one "
              "sub cascade <= "
              "cascadeDiv.\n"
              "\t \t \t \t 3: the order of the SCs is not changed!.\n"
              "\t \t \t \t 4: sorts the SoftClauses, max bucket size in one "
              "cascade <= cascadeDiv.\n"
              "\t \t \t \t 5: SCs are randomly distributed, max bucket size "
              "in one cascade <= "
              "cascadeDiv.\n");
  app.add_set("--interimResult", interimResultParser, {0, 1, 2, 3},
              "\t Interim results are used while processing Cascade.\n"
              "\t \t \t \t 0: No Interim Result is used. (standard)\n"
              "\t \t \t \t 1: Use results to cut at Top.\n"
              "\t \t \t \t 2: Use results to cut at Bottom (best with sorted "
              "weights).\n"
              "\t \t \t \t 3: Use results to cut at both.\n");

  app.add_option("--sepHiWeight", settings->sepHiWeight,
                 "\t if highest weight is bigger than sum of all other SCs, "
                 "process it seperately.\n"
                 "\t \t \t \t same if difference between two weights is "
                 "bigger than 10x.");
  app.add_option("--weightPlus1", settings->weightPlusOne, "\t MISSING COMMENT",
                 true);
  app.add_option("--solveAsWeighted", settings->solveAsWeighted,
                 "\t MISSING COMMENT", true);

  // paper options
  app.add_flag("--adderCaching", settings->adderCaching,
               "\t MISSING COMMENT");  // if argument adderCaching exists then
                                       // true, else false;
  app.add_flag("--coneOfInfluence", settings->coneOfInfluence,
               "\t MISSING COMMENT");
  app.add_option("--exactBounding", settings->exactBounding,
                 "\t Try to solve the last solved weight plus 1.", true);
  app.add_flag("--plain", settings->plainVariant, "\t MISSING COMMENT");
  app.add_option("--DGPW", settings->dGPW, "\t MISSING COMMENT", true);

  // at the moment for divideByColum features
  app.add_option("-t, --featureTest", settings->featureTest,
                 "\t MISSING COMMENT", true);

  app.add_set(
      "--divideDGPW", divideDGPWParser, {0, 1, 2, 3},
      "\t Divide the DGPW into smaller subproblems if possible!\n"
      "\t \t \t \t 0: Do not divide DGPW. (standard)\n"
      "\t \t \t \t 1: Divide only if sum of weights is strictly smaller "
      "than the next set.\n"
      "\t \t \t \t 2: Divide at any point, use interim results to "
      "reduce size of next DGPW.\n"
      "\t \t \t \t 3: Same as two, but only solve Watchdogs from "
      "higher cascades.\n"
      "\t \t \t \t (4: Solving cascades greedy at first then use "
      "additional cascade to "
      "increase weight further if possible.)\n");

  app.add_option(
      "--minSize", settings->minSize,
      "\t Percentage of the number of SCs for which a SC can be splittet.", 0);

  app.add_set(
      "--divisionMode", settings->divisionMode, {-1, 0, 1, 2, 3, 4, 5, 6, 10},
      "\t Divide the DGPW into smaller subproblems if possible!\n"
      "\t \t \t \t -1: Test all modes, choose first suitable. (standard)\n"
      "\t \t \t \t 0: StrictlyBigger. (standard)\n"
      "\t \t \t \t 10: BiggerEqual, smallest bigger equal subset > 1. "
      "(standard)\n"
      "\t \t \t \t 0-6: Different division techniques.\n",
      10);

  app.add_option("--greedyPrepro", settings->greedyPrepro,
                 "\t Activate greedy preprocessor!\n"
                 "\t\t\t 0: without greedy preprocessor\n"
                 "\t\t\t 1: binary search prepro\n"
                 "\t\t\t 2: old greedy variant",
                 1);

  //  app.add_option("--greedyPPTimeLimit", settings->greedyPPTimeLimit,
  //                 "\t Time Limit for the greedy Prepro", 600);
  app.add_option("--greedyPPFixSCs", settings->greedyPPFixSCs,
                 "\t Fix Softclauses either --greedyPrepro=1:\n"
                 "\t\t\t -1: Automatic selection if used, due to number of SCs."
                 "\t\t\t 0: never, only binary search of the never sat SCs\n"
                 "\t\t\t 1: only at the first run\n"
                 "\t\t\t 2: always fixed",
                 0);
  app.add_option("--greedyPPSATPercentage", settings->greedyPPSATPercentage,
                 "\t Should be between 1 and 100", 60);
  app.add_option("--greedyPPUNSATPercentage", settings->greedyPPUNSATPercentage,
                 "\t Should be between 1 and 100", 65);
  app.add_option(
      "--greedyPPPropagationPerSecond", settings->greedyPPPropagationPerSecond,
      "\t How many propagations * timelimit for one PP solver call!", 2500000);
  app.add_option("--greedyPPTimeoutFactor", settings->greedyPPTimeoutFactor,
                 "\t Timeout due to soft clauses X TimeoutFactor!", 12);
  app.add_option(
      "--greedyPPSSCSwitchFactor", settings->greedyPPSSCSwitchFactor,
      "\t After which percentage does it switch from FixSC=2 to FixSC=0 mode!",
      66);

  app.add_option("--greedyMinSizeSCs", settings->greedyMinSizeOfSet,
                 "\t Minimal number of SCs to use GreedyPP.", 100);

  app.add_flag("--useGreedyPreInBetween", settings->useGreedyPreInBetween,
               "\t Activate preprocessor for in between greedy solutions!");

  app.add_set("--testIfDividable", settings->testIfDividable, {0, 1, 2},
              "\t Wcnf is dividable (division mode 1!)\n"
              "\t\t\t 0: no test\n"
              "\t\t\t 1: test if the weight difference is strictly bigger\n"
              "\t\t\t 2: test if the weight difference is bigger equal "
              "(smallest subset for bigger equal has size one)\n");

  app.add_flag("--checkIfSolutionIsUnique", settings->checkIfSolutionIsUnique,
               "\t Test if solution is unique. More Solver calls, chance for "
               "less encoding!");

  app.add_option("--nEqualWeights", settings->atLeastnEqualWeights,
                 "\t At least 8 equal weights - to", 1);

  app.add_option("--partitionFactor", settings->partitionFactor,
                 "\t Factor of difference between two weights to divide "
                 "DGPW, only with divideDGPW>1",
                 10);

  app.add_flag("--DivTest", settings->divCheck, "\t Test if dividable!");

  // Nadja - add here commandline options

  //  }

  CLI11_PARSE(app, argc, argv)

  if (encoding == "dgpw") {
    settings->_encoding = DGPW18;
  } else if (encoding == "warn") {
    settings->_encoding = WARNERS;
  } else if (encoding == "bail") {
    settings->_encoding = BAILLEUX;
  } else if (encoding == "asin") {
    settings->_encoding = ASIN;
  } else if (encoding == "ogaw") {
    settings->_encoding = OGAWA;
  } else if (encoding == "bailw2") {
    settings->_encoding = BAILLEUXW2;
  } else if (encoding == "wmt") {
    settings->_encoding = WMTO;
  } else if (encoding == "mrwto") {
    settings->_encoding = MRWTO;
  } else if (encoding == "mrwto2") {
    settings->_encoding = MRWTO2;
  } else if (encoding == "qauto") {
    settings->_encoding = HEURISTICQMAXSAT;
  } else if (encoding == "auto") {
    settings->_encoding = HEURISTIC20;
  }

  settings->percentOff =
      (settings->percentOff > 100)
          ? 100
          : settings->percentOff;  // man könnte auch ein check(CLI::Range(min,
                                   // max) machen.

  if (partitionStrategyParser == 0) {
    settings->partitionStrategy = NOPARTITION;
  } else if (partitionStrategyParser == 1) {
    settings->partitionStrategy = GROUPBYWEIGHTADDATLAST;
  } else if (partitionStrategyParser == 2) {
    settings->partitionStrategy = GROUPBYWEIGHT;
  } else if (partitionStrategyParser == 3) {
    settings->partitionStrategy = GROUPBYBIGGESTREPEATINGENTRY;
  } else if (partitionStrategyParser == 4) {
    settings->partitionStrategy = GROUPBYCOLUMNS;
  }

  if (encode01Parser == "true") {
    settings->encode01 = true;
    settings->lastPos1 = false;
  }
  if (encode01Parser == "false") {
    settings->encode01 = false;
    settings->lastPos1 = false;
  }
  if (encode01Parser == "lastPos1") {
    settings->encode01 = false;
    settings->lastPos1 = true;
  }

  if (encodeStrategyParser == 0) {
    settings->encodeStrategy = ENCODEALL;
  }

  if (encodeStrategyParser == 0) {
    settings->encodeStrategy = ENCODEALL;
  } else if (encodeStrategyParser == 1) {
    settings->encodeStrategy = ENCODEONLYIFNEEDED;
  }

  if (mcDivideStrategyParser == 0) {
    settings->mcDivideStrategy = SOLVEINNORMALCASCADEMODE;
  } else if (mcDivideStrategyParser == 1) {
    settings->mcDivideStrategy = SORTEDNUMBEROFSOFTCLAUSES;
  } else if (mcDivideStrategyParser == 2) {
    settings->mcDivideStrategy = RANDOMNUMBEROFSOFTCLAUSES;
  } else if (mcDivideStrategyParser == 3) {
    settings->mcDivideStrategy = SOFTCLAUSESINORDER;
  } else if (mcDivideStrategyParser == 4) {
    settings->mcDivideStrategy = SORTEDGOODDIVIDEPOINTS;
  }

  if (divideDGPWParser == 0) {
    settings->divideDGPW = NODIVISION;
  } else if (divideDGPWParser == 1) {
    settings->divideDGPW = USEONLYGCD;
  } else if (divideDGPWParser == 2) {
    settings->divideDGPW = DIVIDEALL;
  } else if (divideDGPWParser == 3) {
    settings->divideDGPW = DIVIDEALLSOLVEONLYWATCHDOGS;
  } else if (divideDGPWParser == 4) {
    settings->divideDGPW = DIVIDEALLSOLVEGREEDY;
  }

  if (interimResultParser == 0) {
    settings->interimResult = NOINTERIMRESULT;
  } else if (interimResultParser == 1) {
    settings->interimResult = CUTATTOP;
  } else if (interimResultParser == 2) {
    settings->interimResult = CUTATBOTTOM;
  } else if (interimResultParser == 3) {
    settings->interimResult = CUTBOTH;
  }

  if (settings->plainVariant) {
    settings->adderCaching = false;
    settings->coneOfInfluence = false;
    settings->exactBounding = false;
    settings->dGPW = false;
    settings->SetPaperOptions();
  }

  if (settings->adderCaching || settings->coneOfInfluence ||
      settings->exactBounding) {
    settings->SetPaperOptions();
  }

  if (settings->dGPW) {
    settings->adderCaching = true;
    settings->coneOfInfluence = true;
    settings->exactBounding = true;
    settings->plainVariant = false;
    settings->SetPaperOptions();
  }

  std::cout << "c This is Pacose 2020" << std::endl;
  std::cout << "c Based on QMAXSAT 2018, GLUCOSE 4.2.1 and MiniSAT"
            << std::endl;
  if (settings->verbosity > 0)
    std::cout << "c file...................: " << maxCnfFile << std::endl;
  pacose->InitSatSolver(returnSolverType(solver));

  //  externalPrepro *externalPrepro =
  //      new class externalPrepro();  // hier evtl. pointer auf Pacose
  //      übergeben

  //  // Code von Nadja
  //  externalPrepro->ParseWcnfFile();

  //  externalPrepro->StartPreprocessing();

  //  if (externalPrepro->SendPreprocessedFileToPacose()) {

  // own pacose Parser
  //   if (pacose->parseWcnfFile(maxCnfFile)) {

  long long int noVars = parse_DIMACS(&maxCnfFile, pacose, settings->verbosity);
  //    parser from QMaxSAT,
  if (noVars != 0) {
    //        std::cout << "c successfully parsed!" << std::endl;
    int returnValue = pacose->SolveProcedure(noVars);
    //  externalPrepro->RestoreSolution();
    //  externalPrepro->PrintSolution();
    return returnValue;

  } else {
    std::cout << "UNSAT" << std::endl;
    std::cout << "c NOT successfully parsed!" << std::endl;
    return 20;
  }
}
