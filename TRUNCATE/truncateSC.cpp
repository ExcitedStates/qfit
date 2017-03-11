#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <mmdb/mmdb_manager.h>

bool isMC ( const PCAtom pAtm )
{
  return ( pAtm->CheckID ( "CA", NULL, NULL ) || pAtm->CheckID ( "N", NULL, NULL ) || pAtm->CheckID ( "C", NULL, NULL ) || pAtm->CheckID ( "O", NULL, NULL ) || pAtm->CheckID ( "CB", NULL, NULL ) );
}

int main ( int argc, char* argv[] )
{
  if ( argc != 5 )
  {
    std::cout << "Usage: truncateSC pdbin pdbout chainID seqNum." << std::endl;
    exit(1);
  }
  
  PCChain pChn;
  PCResidue pRsd;
  PPCAtom pAtmTbl;
  int nAtms;
  CMMDBManager myCMMDBManager;

  std::string inpdb = argv[1];
  int RC = myCMMDBManager.ReadPDBASCII ( const_cast<char*>( inpdb.c_str() ) );

  pChn = myCMMDBManager.GetChain ( 1, argv[3] );
  myCMMDBManager.PDBCleanup ( PDBCLEAN_ATNAME | PDBCLEAN_INDEX | PDBCLEAN_ALTCODE_STRONG );
  myCMMDBManager.FinishStructEdit ( );

  if ( !pChn )
  {
    std::cout << "No Chain " << argv[3] << " in file." << std::endl;
    exit ( 2 );
  }

  pRsd = pChn->GetResidue ( atoi(argv[4]), "" );

  if ( !pRsd )
  {
    std::cout << "No Residue " << argv[4] << " in chain " << argv[3] << "." << std::endl;
    exit ( 3 );
  }

  pRsd->GetAtomTable ( pAtmTbl, nAtms );

  for ( int i=0; i<nAtms; i++ )
  {
    if ( !isMC ( pAtmTbl[i] )  || ( pAtmTbl[i]->altLoc[0] > 'A' ) )
      pRsd->DeleteAtom ( i );
    else
    {
      pAtmTbl[i]->altLoc[0] = char(0);
      pAtmTbl[i]->occupancy = 1.0;
    }
  }
  
  myCMMDBManager.PDBCleanup ( PDBCLEAN_ATNAME | PDBCLEAN_INDEX | PDBCLEAN_ALTCODE_STRONG );
  myCMMDBManager.FinishStructEdit ( );

  myCMMDBManager.WritePDBASCII ( argv[2] );

  return 0;
}

