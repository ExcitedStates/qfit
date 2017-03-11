#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits>
#include <cmath>
#include <vector>

#include "mmdb/mmdb_manager.h"

bool isAA ( const PCResidue pRes )
{
  std::string resname = pRes->GetResName ( );
//std::cout <<  ( pRes->isAminoacid() || resname == "MSE" ) << std::endl;
  return ( pRes->isAminoacid() || resname == "MSE" );
}

bool isAAGEQ7Chain ( const PCChain pC, int& seq_start, int& seq_stop )
{
  PPCResidue pRsdTbl;
  int nRsds;
  bool success = false;
  
  pC->GetResidueTable ( pRsdTbl, nRsds );
  if ( nRsds < 7 ) return success;
  
  seq_start = -100;
  seq_stop = -100;
  
  int i;
  for ( i=0; i<nRsds; ++i )
  {
    if ( isAA ( pRsdTbl[i] ) )
    {
       if ( seq_start == -100 )
         seq_start = pRsdTbl[i]->GetSeqNum();
    }
    else
      if ( seq_start != -100 )
        break;
  }

  if ( seq_stop == -100 && seq_start != -100 )
    seq_stop = pRsdTbl[i-1]->GetSeqNum();
  
  if ( seq_stop - seq_start >= 6 ) 
    // sequence difference of 6 means 7 residues
    success = true;
  return success;
}

int main ( int argc, char* argv[] )
{
  int start, stop;
  if ( argc != 2 )
  {
    std::cout << "Usage: chains pdbin" << std::endl;
    exit(1);
  }
  
  PPCChain pChnTbl;
  int nChns;

  CMMDBManager myCMMDBManager;

  std::string inpdb = argv[1];
  int RC = myCMMDBManager.ReadPDBASCII ( const_cast<char*>( inpdb.c_str() ) );

  myCMMDBManager.GetChainTable ( 1, pChnTbl, nChns );
  myCMMDBManager.FinishStructEdit ( );

  for ( int i=0; i<nChns; ++i )
  {
    if ( isAAGEQ7Chain ( pChnTbl[i], start, stop ) )
      std::cout << pChnTbl[i]->GetChainID() << " " << start << " " << stop << std::endl;
  }

  return 0;
}

