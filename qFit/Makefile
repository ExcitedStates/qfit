ROOT = /Users/zundert/Workspace/qfit-dependencies
CELLAR = /usr/local/Cellar

# Header file directories
LOOPTK_INCLUDE = $(ROOT)/looptk/src/core \
		 $(ROOT)/looptk/src/utils
CGAL_INCLUDE = $(CELLAR)/cgal/4.9/include
CLIPPER_INCLUDE = $(ROOT)/include
CPLEX_INCLUDE = $(ROOT)/cplex-12.6/cplex/include \
		$(ROOT)/cplex-12.6/concert/include
XML2_INCLUDE = /usr/include/libxml2

INCLUDE_DIR = $(CGAL_INCLUDE) $(LOOPTK_INCLUDE) $(CLIPPER_INCLUDE) \
	      $(CPLEX_INCLUDE) $(XML2_INCLUDE)

# Library directories
CLIPPER_LIB = $(ROOT)/lib
LOOPTK_LIB = $(ROOT)/looptk/lib
GSL_LIB = $(ROOT)/lib
CPLEX_LIB = $(ROOT)/cplex-12.6/cplex/lib/x86-64_osx/static_pic \
	    $(ROOT)/cplex-12.6/concert/lib/x86-64_osx/static_pic
CGAL_LIB = $(CELLAR)/cgal/4.9/lib

LIBDIR = $(LOOPTK_LIB) $(CGAL_LIB) $(CLIPPER_LIB) $(CPLEX_LIB)

# Libraries
CPLEX_LIBRARY = ilocplex concert cplex
LOOPTK_LIBRARY = looptk
CGAL_LIBRARY = CGAL
CLIPPER_LIBRARY = clipper-ccp4 clipper-mmdb clipper-contrib clipper-core mmdb
XML2_LIBRARY = xml2
GSL_LIBRARY = gsl

LIB = $(CGAL_LIBRARY) $(CPLEX_LIBRARY) $(LOOPTK_LIBRARY) \
      $(CLIPPER_LIBRARY) $(XML2_LIBRARY) gsl

# Variables determining optimization level and presence of
# debugging information in the output.  The default is to
# activate full optimization (-O3) and omit debugging symbols.
# For best debugging (e.g. with gdb or valgrind), override with
# DEBUG_VAR = -g and OPT_VAR = -O0.
#DEBUG_VAR = -g 

OPT_VAR = -O3 -DMATH_DOUBLE -DIL_STD -DLOOPTK_LIGHT -std=c++98

# The primary flags for the compiler and linker.
CPPFLAGS = $(DEBUG_VAR) $(OPT_VAR) $(addprefix -I,$(INCLUDE_DIR))
CXX = g++-6

LDFLAGS = -pthread $(addprefix -L, $(LIBDIR)) $(addprefix -l,$(LIB))

OBJECTS = $(patsubst src/%.cpp,src/%.o,$(wildcard src/*.cpp)) 
EXECUTABLE = qFit

# Build targets.

all: $(EXECUTABLE)

qFit : $(OBJECTS) 
	$(CXX) -o bin/$@ $(OBJECTS) $(LDFLAGS)

reDokkum_frag : $(OBJECTS)
	$(CXX) -o bin/$@ $(OBJECTS) $(LDFLAGS)

docs:
	./scandoc -i ./default.pl *.h
	mv *.html $(DOCDIR)

default: $(EXECUTABLE)
 
# Cleaning targets.

clean :
	rm -f src/*.o

realclean : clean
	rm -f $(DOCDIR)/*.html

