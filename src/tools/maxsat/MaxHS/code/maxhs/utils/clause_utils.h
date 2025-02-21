/***********[clause_utils.h]
Copyright (c) 2022 Fahiem Bacchus

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

#ifndef CLAUSE_UTILS_H
#define CLAUSE_UTILS_H
#include <vector>

#ifdef GLUCOSE
#include "glucose/core/SolverTypes.h"
#else
#include "minisat/core/SolverTypes.h"
#endif

#ifdef GLUCOSE
namespace Minisat = Glucose;
#endif

namespace MaxHS {

inline bool rm_dups_and_check_tautology(std::vector<Lit> cls) {
  // remove duplicates from clause. If the clause is a tautology
  // return false (here is no clause).
  if (cls.size() > 1) {
    std::sort(cls.begin(), cls.end());
    size_t cur_size{1}, examine{1};
    for (; examine < cls.size(); ++examine) {
      if (cls[cur_size - 1] == cls[examine])
        continue;
      else {
        if (cls[cur_size - 1] == ~cls[examine])
          return false;
        else
          cls[cur_size++] = cls[examine];
      }
    }
    cls.resize(cur_size);
  }
  return true;
}

}  // namespace MaxHS
#endif
