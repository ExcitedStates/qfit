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

//#include "PMain.h"

#ifndef CFRAGMENTSELECTOR_H_
#include "CFragmentSelector.h"
#endif

#ifndef CSIDECHAINSAMPLER_H_
#include "CSideChainSampler.h"
#endif

// NOTE TO WHOEVER COMPILES THIS NEXT:
// This same file is used to compile both the qFit and reDokkum_frag executables,
// depending on whether the REDOKKUM_FRAG variable is defined.
// If it is NOT defined, do "make qFit"          for the qFit executable.
// If it IS     defined, do "make reDokkum_frag" for the reDokkum_frag executable.
//#define REDOKKUM_FRAG

int main(int argc, char *argv[])
{
    WarningIndicator wi = SUPPRESS_WARNINGS;

    bool flip_peptide = false;
    string pdb, mtz, chn, pdb2;
    int idx, hidx;

    string s_flip_peptide = argv[1];
    if ( s_flip_peptide == "--flip-peptide" )
    {
        flip_peptide = true;
        pdb = argv[2];
        mtz = argv[3];
        chn = argv[4];
        idx = atoi(argv[5]);
        #ifdef REDOKKUM_FRAG
            hidx = atoi(argv[6]);
            pdb2 = argv[7];
        #endif
    }
    else if ( s_flip_peptide.find("pdb") == std::string::npos )
    {
        cout << "First argument should be PDB file or --flip-peptide\n";
        exit(0);
    }
    else
    {
        flip_peptide = false;
        pdb = argv[1];
        mtz = argv[2];
        chn = argv[3];
        idx = atoi(argv[4]);
        #ifdef REDOKKUM_FRAG
            hidx = atoi(argv[5]);
            pdb2 = argv[6];
        #endif
    }

    std::string resourcedir = getenv ( "LOOPTK_RESOURCEDIR" );
    LoopTK::Initialize(wi, resourcedir );
	
    if(argc<=1) {
        cerr << "ERROR: Must specify parameter file" << endl;
        exit(1);
    }

    string dbgfile;

#ifndef REDOKKUM_FRAG
    CSideChainSampler mCSideChainSampler;
    mCSideChainSampler.Init(pdb, mtz, chn, idx, flip_peptide);
    mCSideChainSampler.DoFindConfs(chn, idx, idx);
#else
   CFragmentSelector mCFragmentSelector;
   mCFragmentSelector.Init(pdb, mtz, chn, idx, hidx);
   mCFragmentSelector.DoFindConfs(pdb2, chn, idx, hidx);
#endif
    
   cerr << "Initialization done" << endl;
   PResources::FreeResources();
   return 0; 
}
