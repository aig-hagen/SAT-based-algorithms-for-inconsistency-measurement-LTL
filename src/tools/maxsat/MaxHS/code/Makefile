###################################################################################################
#  make			 "for release version"
#  make d		 "debug version"
#  make p		 "profiling version"
#  make clean   	 "remove object files and executable"
###################################################################################################
.PHONY:	r d p rlns dlns plns lr ld lp all install install-bin clean distclean
all: r rlns

## Configurable options ###########################################################################

## Cplex library location (configure these variables)
LINUX_CPLEXLIBDIR   ?= /w/63/fbacchus/CPLEX_Studio201/cplex/lib/x86-64_linux/static_pic
LINUX_CPLEXINCDIR   ?= /w/63/fbacchus/CPLEX_Studio201/cplex/include

## Cplex libary on mac vm
LINUXVM_CPLEXLIBDIR   ?= /home/fbacchus/CPLEX_Studio201/cplex/lib/x86-64_linux/static_pic
LINUXVM_CPLEXINCDIR   ?= /home/fbacchus/CPLEX_Studio201/cplex/include
#
#If you want to build on macos
#DARWIN_CPLEXLIBDIR   ?= /Users/fbacchus/Applications/IBM/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_osx/static_pic/
#DARWIN_CPLEXINCDIR   ?= /Users/fbacchus/Applications/IBM/ILOG/CPLEX_Studio1210/cplex/include
DARWIN_CPLEXLIBDIR   ?= /Users/fbacchus/Applications/IBM/ILOG/CPLEX_Studio201/cplex/lib/x86-64_osx/static_pic/
DARWIN_CPLEXINCDIR   ?= /Users/fbacchus/Applications/IBM/ILOG/CPLEX_Studio201/cplex/include


ifeq "$(shell uname)" "Linux"
CPLEXLIBDIR   =$(LINUX_CPLEXLIBDIR)
CPLEXINCDIR   =$(LINUX_CPLEXINCDIR)
endif
ifeq "$(shell uname)" "Darwin"
CPLEXLIBDIR   =$(DARWIN_CPLEXLIBDIR)
CPLEXINCDIR   =$(DARWIN_CPLEXINCDIR)
endif
ifeq "$(shell hostname)" "ubuntu-vm"
CPLEXLIBDIR   =$(LINUXVM_CPLEXLIBDIR)
CPLEXINCDIR   =$(LINUXVM_CPLEXINCDIR)
endif

# Directory to store object files, libraries, executables, and dependencies:
BUILD_DIR      ?= build

# Include debug-symbols in release builds?
#MAXHS_RELSYM ?= -g

MAXHS_REL    ?= -O3 -DNDEBUG -DNCONTRACTS -DNTRACING
# for cadical if needed
#MAXHS_REL    += -DLOGGING
#MAXHS_REL   += -DQUITE
#MAXHS_DEB    ?= -O0 -DDEBUG -D_GLIBCXX_DEBUG -ggdb
MAXHS_DEB    ?= -O0 -DDEBUG -ggdb
MAXHS_PRF    ?= -O3 -DNDEBUG


ifneq "$(CXX)" "clang++"
MAXHS_REL    += -flto 
endif
# Target file names
MAXHS      = maxhs#          Name of Maxhs main executable.
MAXHS_SLIB = lib$(MAXHS).a#  Name of Maxhs static library.
MAXHSLNS   = maxhs_lns#      Name of Maxhs large neighbour search incomplete solver

#-DIL_STD is a IBM/CPLEX issue

MAXHS_CXXFLAGS = -DIL_STD -I. -I$(CPLEXINCDIR)
MAXHS_CXXFLAGS += -Wall -Wno-parentheses -Wextra -Wno-deprecated
MAXHS_CXXFLAGS += -std=c++17

MAXHS_LDFLAGS  = -Wall -lz -L$(CPLEXLIBDIR) -lcplex -lpthread -ldl -flto

ECHO=@echo

ifeq ($(VERB),)
VERB=@
else
VERB=
endif

#files
ALL_SRCS = $(wildcard maxhs/sys/*.cc) $(wildcard maxhs/core/*.cc) $(wildcard maxhs/ifaces/*.cc) $(wildcard maxhs/utils/*.cc) minisat/utils/System.cc
SATSOLVER_SRCS = $(wildcard cadical/src/*.cpp)
MAXPRE_SRCS = $(wildcard maxpre/src/*.cpp)

ALL_SRCS += $(SATSOLVER_SRCS)
ALL_SRCS += $(MAXPRE_SRCS)

#accomodate cadical .cpp files
TMP = $(ALL_SRCS:.cpp=.cc)
ALL_OBJS = $(TMP:.cc=.o)

MAXHS_SRCS = $(filter-out %maxhslns.cc %lslnssolver.cc, $(ALL_SRCS))
MAXHS_OBJS = $(filter-out %maxhs.o %maxhslns.o %lslnssolver.o, $(ALL_OBJS))

MAXHSLNS_SRCS = $(filter-out %maxhs.cc, $(ALL_SRCS))
MAXHSLNS_OBJS = $(filter-out %maxhs.o %maxhslns.o, $(ALL_OBJS))

r:	$(BUILD_DIR)/release/bin/$(MAXHS)
d:	$(BUILD_DIR)/debug/bin/$(MAXHS)
p:	$(BUILD_DIR)/profile/bin/$(MAXHS)

rlns:	$(BUILD_DIR)/release/bin/$(MAXHSLNS)
dlns:	$(BUILD_DIR)/debug/bin/$(MAXHSLNS)
plns:	$(BUILD_DIR)/profile/bin/$(MAXHSLNS)


lr:	$(BUILD_DIR)/release/lib/$(MAXHS_SLIB)
ld:	$(BUILD_DIR)/debug/lib/$(MAXHS_SLIB)
lp:	$(BUILD_DIR)/profile/lib/$(MAXHS_SLIB)

## Build-type Compile-flags:
$(BUILD_DIR)/release/%.o:           MAXHS_CXXFLAGS +=$(MAXHS_REL) $(MAXHS_RELSYM)
$(BUILD_DIR)/debug/%.o:				MAXHS_CXXFLAGS +=$(MAXHS_DEB)
$(BUILD_DIR)/profile/%.o:			MAXHS_CXXFLAGS +=$(MAXHS_PRF) -pg

## Build-type Link-flags:
$(BUILD_DIR)/profile/bin/$(MAXHS):		MAXHS_LDFLAGS += -pg
ifeq "$(shell uname)" "Linux"
$(BUILD_DIR)/release/bin/$(MAXHS):		MAXHS_LDFLAGS += --static -z muldefs
#$(BUILD_DIR)/release/bin/$(MAXHS):		MAXHS_LDFLAGS += -z muldefs
endif
$(BUILD_DIR)/release/bin/$(MAXHS):		MAXHS_LDFLAGS += $(MAXHS_RELSYM)

$(BUILD_DIR)/profile/bin/$(MAXHSLNS):		MAXHS_LDFLAGS += -pg
ifeq "$(shell uname)" "Linux"
$(BUILD_DIR)/release/bin/$(MAXHSLNS):		MAXHS_LDFLAGS += --static -z muldefs
endif
$(BUILD_DIR)/release/bin/$(MAXHSLNS):		MAXHS_LDFLAGS += $(MAXHS_RELSYM)

## Executable dependencies
$(BUILD_DIR)/release/bin/$(MAXHS):	 	$(BUILD_DIR)/release/maxhs/core/maxhs.o $(foreach o,$(MAXHS_OBJS),$(BUILD_DIR)/release/$(o))
$(BUILD_DIR)/debug/bin/$(MAXHS):	 	$(BUILD_DIR)/debug/maxhs/core/maxhs.o $(foreach o,$(MAXHS_OBJS),$(BUILD_DIR)/debug/$(o))
$(BUILD_DIR)/profile/bin/$(MAXHS):	 	$(BUILD_DIR)/profile/maxhs/core/maxhs.o $(foreach o,$(MAXHS_OBJS),$(BUILD_DIR)/profile/$(o))
#$(BUILD_DIR)/debug/bin/$(MAXHS):	 	$(BUILD_DIR)/debug/maxhs/core/maxhs.o $(BUILD_DIR)/debug/lib/$(MAXHS_SLIB)
#$(BUILD_DIR)/profile/bin/$(MAXHS):	 	$(BUILD_DIR)/profile/maxhs/core/maxhs.o $(BUILD_DIR)/profile/lib/$(MAXHS_SLIB)

$(BUILD_DIR)/release/bin/$(MAXHSLNS):	 	$(BUILD_DIR)/release/maxhs/core/maxhslns.o $(foreach o,$(MAXHSLNS_OBJS),$(BUILD_DIR)/release/$(o))
$(BUILD_DIR)/debug/bin/$(MAXHSLNS):	 	$(BUILD_DIR)/debug/maxhs/core/maxhslns.o $(foreach o,$(MAXHSLNS_OBJS),$(BUILD_DIR)/debug/$(o))
$(BUILD_DIR)/profile/bin/$(MAXHSLNS):	 	$(BUILD_DIR)/profile/maxhs/core/maxhslns.o $(foreach o,$(MAXHSLNS_OBJS),$(BUILD_DIR)/profile/$(o))

## Library dependencies
#$(BUILD_DIR)/release/lib/$(MAXHS_SLIB):	$(foreach o,$(MAXHS_OBJS),$(BUILD_DIR)/release/$(o))
#$(BUILD_DIR)/debug/lib/$(MAXHS_SLIB):		$(foreach o,$(MAXHS_OBJS),$(BUILD_DIR)/debug/$(o))
#$(BUILD_DIR)/profile/lib/$(MAXHS_SLIB):	$(foreach o,$(MAXHS_OBJS),$(BUILD_DIR)/profile/$(o))

## Compile rules 
$(BUILD_DIR)/release/%.o:	%.cc
	$(ECHO) Compiling: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $(MAXHS_CXXFLAGS) $(CXXFLAGS) -c -o $@ $< -MMD -MF $(BUILD_DIR)/release/$*.d
$(BUILD_DIR)/release/%.o:	%.cpp
	$(ECHO) Compiling: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $(MAXHS_CXXFLAGS) $(CXXFLAGS) -c -o $@ $< -MMD -MF $(BUILD_DIR)/release/$*.d

$(BUILD_DIR)/profile/%.o:	%.cc
	$(ECHO) Compiling: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $(MAXHS_CXXFLAGS) $(CXXFLAGS) -c -o $@ $< -MMD -MF $(BUILD_DIR)/profile/$*.d
$(BUILD_DIR)/profile/%.o:	%.cpp
	$(ECHO) Compiling: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $(MAXHS_CXXFLAGS) $(CXXFLAGS) -c -o $@ $< -MMD -MF $(BUILD_DIR)/release/$*.d

$(BUILD_DIR)/debug/%.o:	%.cc
	$(ECHO) Compiling: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $(MAXHS_CXXFLAGS) $(CXXFLAGS) -c -o $@ $< -MMD -MF $(BUILD_DIR)/debug/$*.d
$(BUILD_DIR)/debug/%.o:	%.cpp
	$(ECHO) Compiling: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $(MAXHS_CXXFLAGS) $(CXXFLAGS) -c -o $@ $< -MMD -MF $(BUILD_DIR)/release/$*.d

## Linking rule
$(BUILD_DIR)/release/bin/$(MAXHS) $(BUILD_DIR)/debug/bin/$(MAXHS) $(BUILD_DIR)/profile/bin/$(MAXHS):
	$(ECHO) Linking Binary: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $^ $(MAXHS_LDFLAGS) $(LDFLAGS) -o $@

$(BUILD_DIR)/release/bin/$(MAXHSLNS) $(BUILD_DIR)/debug/bin/$(MAXHSLNS) $(BUILD_DIR)/profile/bin/$(MAXHSLNS):
	$(ECHO) Linking Binary: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(CXX) $^ $(MAXHS_LDFLAGS) $(LDFLAGS) -o $@

## Static Library rule
%/lib/$(MAXHS_SLIB):
	$(ECHO) Linking Static Library: $@
	$(VERB) mkdir -p $(dir $@)
	$(VERB) $(AR) -rcs $@ $^

clean:
	rm -f $(foreach t, release debug profile, $(foreach o, $(ALL_OBJS), $(BUILD_DIR)/$t/$o)) \
      $(foreach t, release debug profile, $(foreach d, $(ALL_OBJS:.o=.d), $(BUILD_DIR)/$t/$d)) \
	  $(foreach t, release debug profile, $(BUILD_DIR)/$t/bin/$(MAXHS))\


## Include generated dependencies
DEPS = $(foreach b, release, $(foreach d, $(ALL_OBJS:.o=.d), $(BUILD_DIR)/$b/$d))

-include $(DEPS)
