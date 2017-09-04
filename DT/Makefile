# Makefile
# CMS DT Pattern Recognition - PR
# Author: 20170515 Sara Vanini

# ROOT
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

# BOOST
BOOSTLIBS = /usr/lib/x86_64-linux-gnu/

# Linux
CXX           = g++ -m64
CXXFLAGS      = -w -Wno-deprecated -fPIC
CXXFLAGS     += $(ROOTCFLAGS)
# default
CXXFLAGS += -O0
# debugging option
#CXXFLAGS += -DDEBUG -gdwarf-3 -g3 -O2
# include root
CXXFLAGS += -I$(ROOTSYS)/include
# libs
LIBS           = $(ROOTLIBS) -L$(BOOSTLIBS) -lboost_program_options  -lm -ldl -rdynamic
GLIBS        = $(ROOTGLIBS)

# files
TOP := $(realpath ../)
CPP_FILES :=  $(wildcard $(TOP)/DT/src/*.cpp)
CPP_FILES +=  $(wildcard $(TOP)/Common/src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

#---- rules --------------------------------------------------------------------------
all: runPR

runPR: $(OBJ_FILES) main.o
	${CXX}  -o runPR ${CXXFLAGS} main.o $(OBJ_FILES) $(LIBS) $(GLIBS)

obj/%.o: $(TOP)/*/src/%.cpp
	${CXX} ${CXXFLAGS} -c -o $@ $<

obj/%.o: $(TOP)/src/%.cpp
	${CXX} ${CXXFLAGS} -c -o $@ $<

main.o:
	${CXX} ${CXXFLAGS} -c main.c

clean:
	rm $(OBJ_FILES) main.o runPR
