# Define C++ compiler
CXX = g++

# Define linker
LINKER = g++

# Debug flag: either nothing or -ggdb
DEBUGFLAG = -ggdb -g3

# Opitmization flag
OPTFLAG = -O2

# Path to StPicoEvent (either absolute or relative)
PICOPATH = /home/lirikk/StPicoEvent

# ROOT flag: -D__ROOT__
ROOTFLAG =

# Define C++ flags (-D_VANILLA_ROOT_ is needed to avoid StMessMgr confusion)
# -std=c++11
CXXFLAGS = $(shell root-config --cflags) -I. -I$(PICOPATH) $(OPTFLAG) -fPIC -Wall -pipe
CXXFLAGS += -D_VANILLA_ROOT_ $(DEBUGFLAG)

# Define libraries to be added
LIBS = $(shell root-config --libs) 
LIBS += -L$(PICOPATH) -lStPicoDst

# Define path to include directories
INCS = $(shell root-config --incdir)

# Linker flags //  -D__ROOT__
LDFLAGS = $(shell root-config --cflags) -I. -I$(PICOPATH) $(OPTFLAG) -fPIC -Wall -pipe -D_VANILLA_ROOT_

all: picoAnalyzerStandalone

# Compile C++ files
picoAnalyzerStandalone: picoAnalyzerStandalone.o
	$(LINKER) $< $(LDFLAGS) $(LIBS) $(GLIBS) -o $@

# Compile C++ files
%.o: %.C
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	@echo 'Compiling: $@ $<'

# PHONY options
.PHONY: clean distclean

# Remove object files and dictionary
clean:
	rm -vf *.o

# Remove object files, dictionary and library
distclean:
	rm -vf *.o picoAnalyzerStandalone
