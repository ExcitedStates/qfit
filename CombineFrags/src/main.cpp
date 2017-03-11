#ifndef CRANKFRAGS_H_
#include "CRankFrags.h"
#endif

int readModelFile ( PCMMDBManager pMMDBManager, const std::string infile )
{
  int RC;
  int kin = isMMDBBIN ( const_cast<char*>( infile.c_str() ) );
  
  if ( kin<0 )  
  {
    printf ( " ***** INPUT FILE\n\n %s\n\nDOES NOT EXIST.", const_cast<char*>( infile.c_str() ) );
    return 2;
  }
  if ( kin==0 )  
    kin = MMDB_FILE_Binary;
  else if ( isPDB(const_cast<char*>( infile.c_str() ))==0 )  
    kin = MMDB_FILE_PDB;
  else if (isCIF(const_cast<char*>( infile.c_str() ))==0)  
    kin = MMDB_FILE_CIF;
  else  
  {
    printf ( " ***** UNRECOGNIZABLE FORMAT OF INPUT FILE\n       %s\n", const_cast<char*>( infile.c_str() ) );
        return 3;
    }
  
  switch ( kin )  
  {
    default               :
    case MMDB_FILE_PDB    : 
      //printf ( " ...reading PDB file %s\n", const_cast<char*>( infile.c_str() ) );
      RC = pMMDBManager->ReadPDBASCII ( const_cast<char*>( infile.c_str() ) );
      break;
    case MMDB_FILE_CIF    :
      printf ( " ...reading CIF file %s\n", const_cast<char*>( infile.c_str() ) );
      RC = pMMDBManager->ReadCIFASCII ( const_cast<char*>( infile.c_str() ) );
      break;
    case MMDB_FILE_Binary :
      printf ( " ...reading MMDB binary file %s\n", const_cast<char*>( infile.c_str() ) );
      RC = pMMDBManager->ReadMMDBF ( const_cast<char*>( infile.c_str() ) );
   }
  return RC;
}

bool ParseGapInput ( const std::string gs, const std::string ge, std::string& CID, 
  int& strt, int& end )
{
  bool success = false;

  if ( isdigit ( gs[0] ) && isdigit ( ge[0] ) )
  {
    strt = atoi ( gs.c_str ( ) );
    end = atoi ( ge.c_str ( ) );
    success = true;
  }
  else if ( ( isalpha ( gs[0] ) && isalpha ( ge[0] ) ) &&
      ( gs.substr ( 0,1 ) == ge.substr ( 0,1 ) ) )
  {
    CID = gs.substr ( 0,1 );
    strt = atoi ( gs.substr (1, gs.size()-1 ).c_str ( ) );
    end = atoi ( ge.substr (1, ge.size()-1 ).c_str ( ) );
    success = true;
  }

  return success;
}

bool GetSubChain ( PCMMDBManager pMMDBManager, const std::string ChnID, const int startSeqNum, const int endSeqNum )
{
  int selHnd = pMMDBManager->NewSelection();
  
  pMMDBManager->Select ( selHnd, STYPE_ATOM, 0, "*", ANY_RES, "*", 
              ANY_RES, "*", "*","*","*", "*", SKEY_NEW );
  
  pMMDBManager->Select ( selHnd, STYPE_ATOM, 0, const_cast<char*>(ChnID.c_str ( ) ), startSeqNum, "*",
              endSeqNum, "*", "*", "*", "*", "*", SKEY_XOR );
  
  pMMDBManager->DeleteSelObjects ( selHnd );
  pMMDBManager->DeleteSelection (  selHnd );
  pMMDBManager->FinishStructEdit ( );
  return true;
}

bool SetTSubChain ( PCMMDBManager pMMDBManager, const std::string ChnID, const int startSeqNum, const int endSeqNum )
{
  int selHnd = pMMDBManager->NewSelection();
  PPCAtom psel;
  int nsel;
  
  pMMDBManager->Select ( selHnd, STYPE_ATOM, 0, "*", ANY_RES, "*", 
              ANY_RES, "*", "*","*","*", "*", SKEY_NEW );
  
  pMMDBManager->Select ( selHnd, STYPE_ATOM, 0, const_cast<char*>(ChnID.c_str ( ) ), startSeqNum+1, "*",
              endSeqNum-1, "*", "*", "*", "*", "*", SKEY_XOR );
  
  pMMDBManager->GetSelIndex( selHnd, psel, nsel );

  for ( int i=0;i<nsel;i++ )
    psel[i]->SetCoordinates ( psel[i]->x, psel[i]->y, psel[i]->z, 1.0, 20.0 );
  
  pMMDBManager->DeleteSelection (  selHnd );
  pMMDBManager->FinishStructEdit ( );
  return true;
}

int main ( int argc, char* argv[] )
{
  int gapStartSeqNum, gapEndSeqNum;
  CMMDBManager mModelCMMDBManager, mFragCMMDBManager, mReadCMMDBManager;
  PPCModel pMdlTbl;
  int nModels;
  std::string ChnID;
  
  if ( argc < 4 )
  {
    std::cout << "Usage: CombineFrags model.pdb frag1.pdb [frag2.pdb] ... [fragN.pdb]\n" 
    << "CombineFrags --scorebest ResID1 ResID2 model1.pdb [model2.pdb] ... [modelN.pdb]\n"
    << "or CombineFrags --scorepart ResID1 ResID2 model1.pdb [model2.pdb] ... [modelN.pdb]\n";
    return 0;
  }

  CRankFrags mCRankFrags;
  
  if ( std::string ( argv[1] ) == "--scorebest" || std::string ( argv[1] ) == "--scorepart")
  {
  	if ( !ParseGapInput ( argv[2], argv[3], ChnID, gapStartSeqNum, gapEndSeqNum ) )
    {
      std::cout << "Error: Ambiguity in specification of flanking residues on command line." << std::endl;
      return 1;
    }
    readModelFile ( &mModelCMMDBManager, argv[4] );
    for ( int i=5; i<argc; i++ )
    { 
      readModelFile ( &mReadCMMDBManager, argv[i] );
      mFragCMMDBManager.AddModel ( mReadCMMDBManager.GetModel ( 1 ) );
      std::string t ( argv[i] );
      mCRankFrags.CalcMap ( t.erase ( t.find ( ".pdb") , 4 ) + ".mtz");
      SetTSubChain ( &mFragCMMDBManager, ChnID, gapStartSeqNum, gapEndSeqNum );
      mCRankFrags.calcEDMCC ( &mFragCMMDBManager, i-4);
  	  GetSubChain ( &mFragCMMDBManager, ChnID, gapStartSeqNum, gapEndSeqNum );
  	  //mFragCMMDBManager.AddModel ( mReadCMMDBManager.GetModel ( 1 ) );
  	  //mCRankFrags.calcEDMCC ( &mFragCMMDBManager, i-3);
    }
    if ( std::string ( argv[1] ) == "--scorepart" )
    {
      for ( int l=gapStartSeqNum+1; l<gapEndSeqNum; l++ )
        mModelCMMDBManager.GetModel(1)->GetChain ( ChnID.c_str() )->DeleteResidue ( l, "");
      mModelCMMDBManager.PDBCleanup ( PDBCLEAN_INDEX );
      mModelCMMDBManager.FinishStructEdit ( );
      mCRankFrags.InsertBestTerminal ( &mModelCMMDBManager, &mFragCMMDBManager, true );
      mCRankFrags.InsertBestTerminal ( &mModelCMMDBManager, &mFragCMMDBManager, false );
      std::string s ( argv[4] );
      s.insert ( s.find( ".pdb"), "_best" );
      mModelCMMDBManager.WritePDBASCII ( pstr(s.c_str()) ); 
    }
    else
    {
      int bestfrag =  mCRankFrags.Rankum ( &mFragCMMDBManager );
      PCChain pBestFrag  = mFragCMMDBManager.GetModel ( bestfrag )->GetChain ( 0 );
      mFragCMMDBManager.WritePDBASCII ( "tst.pdb" );
      mCRankFrags.InsertFragment ( mModelCMMDBManager.GetModel ( 1 )->GetChain ( pBestFrag->GetChainID ( ) ), pBestFrag );
      mModelCMMDBManager.PDBCleanup ( PDBCLEAN_INDEX );
      mModelCMMDBManager.FinishStructEdit();
      std::string s ( argv[4] );
      mModelCMMDBManager.WritePDBASCII ( pstr(s.c_str()) );
    }
  }
  else
  {
  	readModelFile ( &mModelCMMDBManager, argv[1] ); 
    for ( int i=2; i<argc; i++ )
    {
      readModelFile ( &mReadCMMDBManager, argv[i] );
      mFragCMMDBManager.AddModel ( mReadCMMDBManager.GetModel ( 1 ) );
    }
  
    mFragCMMDBManager.GetModelTable ( pMdlTbl, nModels );
  
    for ( int i=0; i<nModels; i++ )
    {
      CMMDBManager outCMMDBManager;
      outCMMDBManager.Copy ( &mModelCMMDBManager, MMDBFCM_All );
      mCRankFrags.InsertFragment ( outCMMDBManager.GetModel ( 1 )->GetChain ( pMdlTbl[i]->GetChain ( 0 )->GetChainID ( ) ), pMdlTbl[i]->GetChain ( 0 ) );
      std::ostringstream oss;
      oss << i;
      outCMMDBManager.PDBCleanup ( PDBCLEAN_INDEX );
      outCMMDBManager.FinishStructEdit();
      std::string s ( argv[1] );
      std::string t ( argv[2+i] );
      s.insert ( s.find( ".pdb"), "_" + t.erase ( t.find ( ".pdb") , 4 ) );
      outCMMDBManager.WritePDBASCII ( (pstr)s.c_str() );
    }
  }
  return 0;
}