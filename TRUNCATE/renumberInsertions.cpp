#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits>
#include <cmath>
#include <vector>
#include <mmdb/mmdb_manager.h>

int main ( int argc, char* argv[] )
{
  if ( argc != 2 )
  {
    std::cout << "Usage: insertionRenumber pdbin" << std::endl;
    exit(1);
  }
  
  PPCChain pChnTbl;
  int nChns;

  CMMDBManager myCMMDBManager;

  std::string inpdb = argv[1];
  int RC = myCMMDBManager.ReadPDBASCII ( const_cast<char*>( inpdb.c_str() ) );

  myCMMDBManager.GetChainTable ( 1, pChnTbl, nChns );
  myCMMDBManager.FinishStructEdit ( );

  PPCResidue pRsdTbl;
  int nRsds;
 
  for ( int i=0; i<nChns; ++i )
  {
    int diff = 0;
    pChnTbl[i]->GetResidueTable ( pRsdTbl, nRsds );
    for ( int j=0; j<nRsds; ++j )
    {
      if ( (std::string)(pRsdTbl[j]->GetInsCode()) != "" )
      {
        diff++;
       *(pRsdTbl[j])->insCode=' ';
      }
      pRsdTbl[j]->seqNum+=diff;
    }
  }

  myCMMDBManager.FinishStructEdit ( );
  RC = myCMMDBManager.WritePDBASCII ( const_cast<char*>( inpdb.c_str() ) );

  return 0;
}

