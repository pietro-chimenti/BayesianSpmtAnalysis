CXXSRCS = $(wildcard *.cxx) 

PRGSRCS = $(wildcard *.cpp)

# compiler and flags
CXX       = g++
CXXFLAGS  = -g -O2 -Wall -fPIC -Wno-deprecated 
LD        = /usr/bin/ld -m elf_x86_64
LDFLAGS   = -g -O2 

CXXFLAGS += $(shell bat-config --cflags)
LIBS := $(shell bat-config --libs)

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS = $(addsuffix .o,$(basename $(CXXSRCS)))
MYPROGS = $(basename $(PRGSRCS))
PRGOBJS = $(addsuffix .o,$(basename $(PRGSRCS)))

GARBAGE = $(CXXOBJS) $(PRGOBJS) link.d $(MYPROGS) $(wildcard *.pdf) $(wildcard *.root) $(wildcard *.jpg) log.txt

# targets
all : $(MYPROGS)

.PHONY : all clean print

link.d : $(addsuffix .h,$(basename $(CXXSRCS))) $(CXXSRCS) $(PRGSRCS)
	$(CXX) -MM $(CXXFLAGS) $(filter-out %.h,$^) > link.d;
	@$(foreach prog,$(MYPROGS), echo $(prog) : $(prog).o >> link.d;)

-include link.d

$(CXXOBJS) $(PRGOBJS) :
	$(CXX) $(CXXFLAGS) -c $(filter $(basename $@).%,$(filter-out %.h,$^)) -o $@

$(MYPROGS) : $(CXXOBJS)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

clean :
	rm -f $(GARBAGE)

print :
	@echo compiler  : $(CXX)
	@echo c++ srcs  : $(CXXSRCS) $(PRGSRCS)
	@echo c++ objs  : $(CXXOBJS) $(PRGOBJS)
	@echo c++ flags : $(CXXFLAGS)
	@echo ld flags  : $(LDFLAGS)
	@echo libs      : $(LIBS)
