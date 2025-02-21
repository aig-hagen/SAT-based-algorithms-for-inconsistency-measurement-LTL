#ifndef ONEWATCH_H
#define ONEWATCH_H

/***********[onewatch.h]
Copyright (c) 2022, Fahiem Bacchus

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***********/

#include <algorithm>
#include <ostream>

#ifdef GLUCOSE
#include "glucose/core/SolverTypes.h"
#else
#include "minisat/core/SolverTypes.h"
#endif

#include "maxhs/ds/packed.h"
#include "maxhs/utils/io.h"

struct Watch {
  size_t cls_size{0};
  int cls_idx{0};
  Minisat::Lit blocker{Minisat::lit_Undef};
};

inline std::ostream& operator<<(std::ostream& os, const Watch& w) {
  os << "[cls_size = " << w.cls_size << " cls_idx = " << w.cls_idx
     << " blocker = " << w.blocker << "] ";
  return os;
}

class OneWatch {
 public:
  OneWatch(const MaxHS::Packed_vecs<Minisat::Lit> pv) {
    for (size_t i{0}; i < pv.size(); ++i) {
      if (!pv[i].empty()) {
        auto lt = pv[i][0];
        ensure_size(lt);
        o_lists[Minisat::toInt(lt)].push_back(
            {pv[i].size(), static_cast<int>(i),
             (pv[i].size() > 1 ? pv[i][1] : pv[i][0])});
      }
    }
    for (auto& ol : o_lists)
      sort(ol.begin(), ol.end(), [](const Watch& x, const Watch& y) {
        return x.cls_size < y.cls_size;
      });
  }
  const vector<Watch>& operator[](Minisat::Lit lt) const {
    if (static_cast<size_t>(toInt(lt)) < o_lists.size())
      return o_lists[Minisat::toInt(lt)];
    else
      return no_watches;
  }

 private:
  vector<vector<Watch>> o_lists;
  void ensure_size(Minisat::Lit lt) {
    size_t idx = Minisat::toInt(lt);
    if (idx >= o_lists.size()) o_lists.resize(idx + 1);
  }
  vector<Watch> no_watches;
};

#endif
