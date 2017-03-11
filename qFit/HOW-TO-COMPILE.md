How to compile with GCC
=======================

- Change source code of header files where CSolver.h is included. The include
  needs to be put on top (probably due to *string* versus *cstring* include).
- Compile clipper with latest gcc/g++ and use -std=c++98 by setting CXX and
  CXXFLAGS env-variables.
- See Makefile for all required libraries and compiler flags.


Notes
-----
I copied the static library looptk.a to the main lib/ folder and renamed it
liblooptk.a for a more consistent Makefile.

