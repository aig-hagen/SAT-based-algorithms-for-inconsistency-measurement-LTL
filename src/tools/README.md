This directory contains various tools needed to compute beyond-NP measures (MI, MIC, MC, P, MV, NC).

To download the git submodules, issue `git submodule update --init --recursive`.

- cadical: SAT solver used by RIME. Install according to README and update the path to cadical in rime/Makefile before compiling RIME.
- ganak: exact model counter used by MSSCounting and exactMUScounter. Install according to README and copy binary to MSSCounting/tools and exactMUScounter/tools.
- rime: MSS/MCS enumerator. Install according to README and copy binary to exactMUScounter/tools.
- mustool: MUS enumerator. Install according to README.
- exactMUScounter: MUS counter. Install according to README. Remember to copy binaries for ganak and rime to exactMUScounter/tools.
- MSSCounting: MSS/MCS counter. Install according to README. Remember to copy binary for ganak to MSSCounting/tools.
- MARCO: MUS enumerator. Install according to README.
- forqes: Smallest MUS extractor. Binary for x86-64 Linux provided.
- mcscache: MSS/MCS enumerator. Binary for x86-64 Linux provided.
- umuser: Computes the union of MUSes. Binary for x86-64 Linux provided.

Newer GCC versions may complain about `printf` not being defined while compiling MCSMUS (used by RIME and MUST). In this case add `#include <cstdio>` to `control.cc` as suggested by the compiler.
