/*
    qFit: Multiconformer modeling by constrained fitting of rotamer occupancies
    Henry van den Bedem, Ankur Dhanik, Jean-Claude Latombe, Ashley Deacon. Acta Cryst. D65:1107–1117 (2009)
    e-mail: vdbedem@slac.stanford.edu

        Copyright (C) 2009-2012 Stanford University

	Permission is hereby granted, free of charge, to any person obtaining a copy of
	this software and associated documentation files (the "Software"), to deal in
	the Software without restriction, including without limitation the rights to
	use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
	of the Software, and to permit persons to whom the Software is furnished to do
	so, subject to the following conditions: 

	This entire text, including the above copyright notice and this permission notice
	shall be included in all copies or substantial portions of the Software. 

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
	OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
	IN THE SOFTWARE.

    
*/


#ifndef CFRAGMENTSELECTOR_H_
#define CFRAGMENTSELECTOR_H_

#ifndef CSOLVER_H_
#include "CSolver.h"
#endif

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <numeric>
#include "PBasic.h"
#include "PChain.h"
#include "PExtension.h"
#include "PResources.h"
#include "PTools.h"
#include "PMyatomsf.h"
#include "PIKAlgorithms.h"
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-mmdb.h>
#include <mmdb/mmdb_manager.h>
#include "POptimize.h"
#ifndef __RAMPROB_H_
#include "ramprob.h"
#endif
#include "myPDBIO.h"
#include "PScaleXMap.h"

#ifndef CSIDECHAINSAMPLER_H_
#include "CSideChainSampler.h"
#endif

// Types to hold vector-of-ints (Vi) and vector-of-vector-of-ints (Vvi)
typedef std::vector<int> Vi;
typedef std::vector<Vi> Vvi;


class CFragmentSelector
{
private:
  CMMDBManager myCMMDBManager;
  CSolver mCSolver;
  std::vector<std::string> atm_ids;
  int resnum;
  std::vector<std::vector<char> > orderings;
  std::vector<bool> isMSE;
  
  void generate_paths_recursive ( const std::vector<std::vector<char> > &a, std::vector<char> &s, const int D );
  void GetAtomIDs ( const PCChain pChn, const int lowseq, const int highseq );
  Atom_list GetAtoms ( const PCChain pChn, const int lowseq, const int highseq, const std::vector<char> path );
  Atom_list GetMCAtoms ( const PCChain pChn, const int lowseq, const int highseq, const std::vector<char> path );
  void SetAtoms ( const PCChain pChn, const int lowseq, const int highseq, const Atom_list& atm_lst, const char altLoc, const double occupancy );
  void GetAltLocs ( const int lowseq, const int highseq, const PCChain pChn, std::vector<std::vector<char> >& altChars );

public:
	CFragmentSelector();
	virtual ~CFragmentSelector();
	void Init ( std::string inpdb, std::string inmtz, const std::string ChnID, const int lSeqNum, const int hSeqNum );
	void DoFindConfs ( std::string inpdb, const std::string ChnID, const int lSeqNum, const int hSeqNum );
    void cart_product(
        Vvi& rvvi,  // final result
        Vi&  rvi,   // current result
        Vvi::const_iterator me, // current input
        Vvi::const_iterator end); // final input
	
};

#endif /*CFRAGMENTSELECTOR_H_*/
