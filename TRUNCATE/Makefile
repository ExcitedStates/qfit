
ROOT = /Users/zundert/Workspace/qfit-dependencies
INCLUDE_DIR = $(ROOT)/include
LIB_DIR = $(ROOT)/lib
LIBRARY = mmdb

#CPPFLAGS = -xW $(addprefix -I,$(INCLUDE_DIR)) 
CPPFLAGS = $(addprefix -I,$(INCLUDE_DIR)) #-std=c++98 -fpermissive
LDFLAGS = $(addprefix -L,$(LIB_DIR)) $(addprefix -l,$(LIBRARY))
CXX = g++-6
CC = gcc

EXE = $(patsubst %.cpp,%,$(wildcard *.cpp))

all: 	$(EXE)

$(EXE):	
	$(CXX) $(CPPFLAGS) $(addprefix -L,$(LIB_DIR)) $@.cpp $(addprefix -l,$(LIBRARY)) -o $@
