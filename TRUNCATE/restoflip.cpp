#include <iostream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <string>
#include <mmdb/mmdb_manager.h>

int main ( int argc, char *argv[] )
{
  if ( argc != 6 )
  {
    std::cout << "Usage: restoflip [flip|noflip] pdbin chainID firstResNum lastResNum" << std::endl;
    std::cout << "(Make sure pdbin has only 1 conformation!)" << std::endl;
    exit(1);
  }

  PCChain pChn;
  PPCResidue pRsdTbl;
  int nRsds;
  
  CMMDBManager myCMMDBManager;

  std::string flip_or_noflip = argv[1];

  std::string inpdb = argv[2];
  int RC = myCMMDBManager.ReadPDBASCII ( const_cast<char*>( inpdb.c_str() ) );

  myCMMDBManager.GetModel( 1 )->CalcSecStructure ( false );

  pChn = myCMMDBManager.GetChain ( 1, argv[3] );

  if ( ! pChn )
  {
    std::cout << "No Chain " << argv[3] << " in file." << std::endl;
    exit ( 2 );
  }

  int firstResNum = atoi( argv[4] );
  int lastResNum  = atoi( argv[5] );

  std::vector<int> helixOrSheet;
  std::vector<int> other;

  pChn->GetResidueTable ( pRsdTbl, nRsds );

  for ( int i = 0; i < nRsds; i++ )
  {
    if ( ! pRsdTbl[i]->isAminoacid() ) continue;
    if ( pRsdTbl[i]->GetSeqNum() < firstResNum || pRsdTbl[i]->GetSeqNum() > lastResNum ) continue;
    
    if ( pRsdTbl[i]->SSE == SSE_Strand || pRsdTbl[i]->SSE == SSE_Helix )
    {
      helixOrSheet.push_back( pRsdTbl[i]->GetSeqNum() );
    }
    else
    {
      other.push_back( pRsdTbl[i]->GetSeqNum() );
    }
  }

  // Output

  std::vector<int> v;
  if ( flip_or_noflip == "flip" )
    v = other;
  else
    v = helixOrSheet;

  std::string resListStr = "";
  for ( std::vector<int>::iterator it = v.begin(); it != v.end(); ++it )
  {
    if ( it != v.begin() )
    {
      resListStr += ",";
    }
    //resListStr += std::to_string( *it );
    std::ostringstream os;
    os << *it;
    resListStr += os.str();
  }

  std::cout << resListStr << std::endl;

  return 0;
}
