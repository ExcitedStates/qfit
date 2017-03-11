#include <iostream>
#include <stdlib.h>
#include <mmdb/mmdb_manager.h>

int main ( int argc, char *argv[] )
{
  if ( argc != 4 )
  {
    std::cout << "Usage: flipornot pdbin chainID resnum" << std::endl;
    std::cout << "(Make sure pdbin has only 1 conformation!)" << std::endl;
    exit(1);
  }

  PCChain pChn;
  PPCResidue pRsdTbl;
  int nRsds;
  
  CMMDBManager myCMMDBManager;

  std::string inpdb = argv[1];
  int RC = myCMMDBManager.ReadPDBASCII ( const_cast<char*>( inpdb.c_str() ) );

  myCMMDBManager.GetModel( 1 )->CalcSecStructure ( false );

  pChn = myCMMDBManager.GetChain ( 1, argv[2] );

  if ( ! pChn )
  {
    std::cout << "No Chain " << argv[2] << " in file." << std::endl;
    exit ( 2 );
  }

  pChn->GetResidueTable ( pRsdTbl, nRsds );
  
  for ( int i = 0; i < nRsds; i++ )
  {
    if ( ! pRsdTbl[i]->isAminoacid() ) continue;
    if ( pRsdTbl[i]->GetSeqNum() != atoi(argv[3]) ) continue;

    if ( pRsdTbl[i]->SSE != SSE_Strand && pRsdTbl[i]->SSE != SSE_Helix )
    {
      std::cout << "--flip-peptide" << std::endl;
    }
  }

  return 0;
}
