
/********************************************************************************************
 * antom.h -- Copyright (c) 2014-2016, Sven Reimer
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this
 * software and associated documentation files (the "Software"), to deal in the Software
 * without restriction, including without limitation the rights to use, copy, modify, merge,
 * publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
 * to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or
 * substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 * BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 ********************************************************************************************/

#ifndef ANTOM_ANTOM_H_
#define ANTOM_ANTOM_H_

// Include standard headers.
#include "antombase.h"

namespace antom
{
  enum SolverType
  {
    ANTOMSOLVER,
    CRYPTOMINISATSOLVER
  };

  class SoftClause;
  class Sorter;

  // The "MaxAntom" class.
  class Antom : public AntomBase
  {
  public:

    // Constructor.
    Antom(uint32_t threads = 1);

    // Destructor.
    ~Antom(void);

    void SetControl(void);

    /* MaxSAT interface */

    // Returns the last used index. Does not necessary meet with "variables()"
    uint32_t LastIndex(void) const;

    // Add a clause to the soft clause database
    bool AddSoftClause(std::vector<uint32_t>& clause, uint32_t weight = 1);

    // Set predefined bounds
    void SetMaxBounds(uint32_t low, uint32_t high);
    void SetLowerBound(uint32_t low);
    void SetHigherBound(uint32_t high);

    // Set a target for optimal solution.
    // Declares minimal number of unsatisfied clauses
    // If target is reached, stop maxSAT-computation.
    void SetOptTarget(int32_t target);
    // States whether the target should be reached as precise as possible
    void SetPreciseTarget(bool val);

    // Search modes for maxSAT
    void SetSearchMode(SearchMode val);

    // Partial modes for maxSAT
    void SetPartialMode(PartialMode val);

    // Solve partial maxSAT Problem
    uint32_t MaxSolve(int32_t& optimum);
    uint32_t MaxSolve(const std::vector< uint32_t >& externalAssumptions, int32_t& optimum);

    // Solve weighted partial maxSAT Problem
    uint32_t MaxSolveWeightedPartial(int32_t& optimum);
    uint32_t MaxSolveWeightedPartial(const std::vector< uint32_t >& externalAssumptions, int32_t& optimum);

    #if 0
    // Sets back-end SAT-Solver for MaxSAT
    void SetSolverProxy(SolverType type);
    SolverProxy* GetSolverProxy(void) const;
    #endif

    /* Some maxSAT related interface functions */
    void SetMaxWidth(uint32_t width);
    void SetSplittedWidth(uint32_t width);
    void SetMaxParts(uint32_t parts);
    void SetSortSoftClauses(uint32_t val);
    void CalcSplitWidth(void);
    void SetSkip(bool val);
    void SetCSC(uint32_t val);
    void SetRelaxationLits(bool val);
    void SetIncompleteMode(bool val);
    void SetHorizontalWidth(uint32_t val);
    void SetGridMode(uint32_t val);
    void SetNetworktype(SorterType val);
    void SetMaxInprocessing(bool val);

    // val = 0 --> switch off incremental mode for maxSAT (default)
    // val = 1 --> use maxSAT part of SAT incrementally (soft clauses are deleted after each call).
    // val = 2 --> use maxSAT part of SAT incrementally (soft clauses are kept after each call).
    void SetIncrementalMode(uint32_t val);

    // De-/activates constant detection with SAT in preprocessing
    void SetSatConst(uint32_t val);

    // Resets MaxSAT related data
    void Reset(void);
    void InstanceReset(void);

    void SetStrategyforSoftVariables(bool pos);

    void SetEncode01Mode(bool val);
    void SetDecStratMode(uint32_t val);
    void SetBypassesMode(uint32_t val);
    void SetBaseMode(uint32_t val);

    void enableBBintegration();


    // settings for BB integetration
    void setStopAfterFirstConflict(bool flag);

    // tells the solver that it will receive its initial clause DB from the BB solver
    bool receiveInitialClausesFromBB();

    double getTotalRuntimeInSeconds();

    void waitForStart();

    friend class Sorter;
  private:

    void DataReset(void);

    // Solves the current (partial) MaxSAT formula. Returns SAT/UNSAT and modifies "optimum", representing the minimum number of
    // unsatisfied clauses.
    uint32_t MaxSolveLocal(int32_t& optimum);
    uint32_t MaxSolveLocal(const std::vector<uint32_t>& externalAssumptions, int32_t& optimum);

    // Check for constants in the candidates list by applying the SAT solver
    // Attention! Use with care, process might be costly
    // If quickcheck is set, only assumptions are deduced, instead of solving the instance
    // -> less powerful, but much faster
    bool FindConstantsWithSat(std::vector< uint32_t >& candidates, bool quickcheck);

    uint32_t SatBasedSearch(Sorter* sorter, uint32_t& pos);
    uint32_t UnsatBasedSearch(Sorter* sorter, uint32_t& pos);
    uint32_t BinaryBasedSearch(Sorter* sorter, uint32_t& pos);

    uint32_t GetPreciseTarget(Sorter* sorter, uint32_t& pos, int32_t gap);

    bool AddBound(uint32_t lit);

    uint32_t SetLowerPartialParts(std::vector< uint32_t >& localassumptions) const;
    int32_t UpdateSorterOptimum(void);

    uint32_t CountSatisfiedSoftClauses(Sorter* sorter, const std::vector<uint32_t>& model);

    //void mergePartTrigger ( Sorter& trigger1, Sorter& trigger2 );
    void CheckAllConflictingSoftclauses(void);

    // Store already found assignments of relaxation lits to "assumptions"
    void CollectRelaxationLits(std::vector< uint32_t >& assumptions);

    Sorter* CreateNextSorter(uint32_t depth = 0, uint32_t targetdepth = 1);
    Sorter* MergeNextSorter(void);
    Sorter* GetNextSorter(void);

    // Deactives soft clasues, deletes sorter clauses and variables
    void InvalidateSoftClauses(void);

    //Weighted MaxSAT stuff

    // Cout the model of the buckets and which position has the first set bit.
    void DumpBucketModel(const std::vector<uint32_t>& model);

    // Get rough solution by proceeding last bucket
    uint32_t ProceedLastBucket(int32_t& optimum);
    // Get precise solution by adjusting the bucket's tares
    uint32_t AdjustTareValues(int32_t& optimum);

    void SetTareBounds(void);

    void ExitTimeout(void);

    void SetDecisionStrategiesForMaxSAT(void);


    // end settings for BB integetration

    #ifndef NDEBUG
    void CheckGates(void) const;
    #endif

    std::vector< SoftClause* > _softClauses;

    std::vector< std::vector< Sorter* > > _sorterTree;

    std::vector< uint32_t > _bestModel;

    uint32_t _maxSorterDepth;
    SearchMode _searchMode;

    uint32_t _maxsatResult;

    // partial status variables
    PartialMode _partialMode;
    bool _lastpartialsorter;
    uint32_t _splittedWidth;
    uint32_t _maxWidth;
    uint32_t _maxParts;
    uint32_t _horizontalwidth;
    uint32_t _bypassGrid;
    uint32_t _horizontalbypasses;
    uint32_t _verticalbypasses;
    uint32_t _sortsoftclauses;
    uint32_t _comparator;
    uint32_t _skipped;
    uint32_t _setCSC;
    bool _doskipping;
    bool _setRelaxationLits;
    bool _incompleteMode;
    bool _maxInprocess;

    uint32_t _satconst;

    // bitonic sorter
    // odd even sorter
    SorterType _networkType;
    bool _encode01Mode;

    uint32_t _triggervar;

    // last index of original problem (without sorter network variables)
    uint32_t _lastIndex;

    // Bounds for maximization mode
    uint32_t _low;
    uint32_t _high;
    int32_t _target;
    bool _preciseTarget;

    //Weighted MaxSAT stuff
    uint32_t _currentBucketForTare;

    uint32_t _baseMode;
    uint64_t _sumOfSoftWeights;
    uint64_t _satWeight;

    // partial status variables
    uint32_t _highestBucketMultiplicator;
    uint32_t _currentBucketMultiplicator;
    // the value calculated highestBucketMultiplicator * bucketInfo[1];
    uint64_t _estimatedSatWeight;
    uint64_t _diffToSatWeight;
    std::vector<uint32_t> _collectedAssumptions;

    // testing new functions
    int _decStratMode;

  };
}

#endif
