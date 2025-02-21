/***********************************************************************[Hardware_reuse_sorters.cc]
Copyright (c) 2021, Marek Piotrów

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

/* An input sequence to the sorter consists of element of type Formula with possible repetitions.
   Each sequence is sorted and the repetitions are counted. If the count is <= 7, it is stored on 
   the 3 least significant bits of the element (they are all zeroes in the input values); otherwise,
   the counter is added as the next element to an vector with unique elements only. Previous input 
   sequences are stored in this way in the vector prev_elem. A greedy algorithm is used to find previous 
   sequences that are subsequences of the currect one. If found, they are replaced in the current one 
   with a pairs of elements: a constant Subseq and an index of a subsequence in prev_seq. Those special 
   pairs are put at the beginning of compacted current sequence. 
*/

struct OutSequence {   // items of output sequences are stored in the vector out_fmls 
    unsigned index;    // contains index of the next-to-last item of the corresponding output sequence;
    int      ineq;     // the corresponding values of ineq in the encoded constraint
    int      cost;     // contains the cost of the corresponding sorter encoding in FEnv
    OutSequence(unsigned _index, int _ineq, int _cost): index(_index), ineq(_ineq), cost(_cost) {}
} ;

class ReuseSorters {
  private:
    static const unsigned maxCount = 46; // the default value of opt_base_max - 1
    static const int minStoredSeqSize = 16; // the minimum length of a sequence to be stored and reused
    static const int minAnalyzedSeqSize = 24; // the minimum length of a sequence to be analyzed for subsequences
    vec<unsigned> prev_elem; // contains elements of compacted previous sequences;
    vec<int> prev_seq,       // contains indices of the next-to-last elements of prev_elem;
             prev_seq_size;  // contains lengthes of input sequences stored in prev_seq and prev_elem;
    Map<Pair<unsigned, unsigned>, vec<int>* > fnmap; // (input, count) -> { indices in prev_seq, where the pair appered } 
    Map<unsigned, vec<int>* > fnnmap[maxCount]; // as above, but for count c in {1, 2, ..., maxCount} and in fnnmap[c-1]
    vec<Formula> out_fmls;   // if opt_shared_fmls is false, out_fmls contains elements of outputs of previous sorters;
    vec<OutSequence> out_seq;   // if opt_shared_fmls is false, out_seq contains info about stored sorter encodings   
    int tmp_seq_cnt;         // The number of last sequences that can be kept or removed. 
                             // They are removed if the corresponding encoding is replaced by a smaller BDD/Adder one. 

    void updateCoverIndices(unsigned nr, vec<int>& cover, vec<int>& rev_cover, vec<Pair<unsigned,unsigned> >& usedfs, 
                                                          Map<unsigned, unsigned>& unusedmap, bool root_level);
    int  encodeWithReuse(vec<unsigned>&outfs, int from, int to, vec<Formula>& outvars, int k, int ineq, bool soft_constr);
    void insertSeqElemsIntoMaps(unsigned nr, vec<Pair<unsigned,unsigned> >& usedfs, unsigned& reused_size, bool root_level);
  public:
    ReuseSorters(): tmp_seq_cnt(0) {}
    int  encodeBySorter(vec<Formula>& fs, int max_sel, int ineq, bool soft_constr);
    void keepLastSequences(void);
    void removeLastSequences(void);
} ;


