#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits>
#include <cmath>
#include <vector>
#include <mmdb/mmdb_manager.h>

template<class T, class U> struct ComparePair1stDsc 
{
  bool operator()( const std::pair<T,U>& a, const std::pair<T,U>& b ) const 
  {
    return a.first > b.first; 
  }
};

bool isMC ( const PCAtom pAtm )
{
  return ( pAtm->CheckID ( "CA", NULL, NULL ) || pAtm->CheckID ( "N", NULL, NULL ) || pAtm->CheckID ( "C", NULL, NULL ) || pAtm->CheckID ( "O", NULL, NULL ) );
}

bool isOK ( const bool ok )
{
  return ok;
}


void RotateFinalChi ( CMMDBManager &mCMMDBManager, const int seqnum, AltLoc aloc )
{
  PCAtom head, tail;
  int selHnd, natms;
  PPCAtom atms;
  PCResidue pres = mCMMDBManager.GetChain ( 1, "A" )->GetResidue ( seqnum, "" );
  std::string name = pres->GetResName ( );
  selHnd = mCMMDBManager.NewSelection();
  
  if ( name == "ASN" || name == "ASP" || name == "TYR" || name == "HIS" || name == "PHE" )
  {
    tail = pres->GetAtom ( "CB", NULL, aloc );
    head = pres->GetAtom ( "CG", NULL, aloc );
    
    mCMMDBManager.Select (
      selHnd,   // selection handle
      STYPE_ATOM,  // select atoms
      0,           // any model
      "A",         // chains "A" and "B" only
      seqnum,"*", // this residue
      seqnum,"*", //   any insertion code
      "*",         // any residue name
      "!N,CA,C,O,CB,CG", // C-alphas only
      "*",         // carbons only
      aloc,        // any alternative location indicator
      SKEY_NEW     // NEW-selection
    );
  }
  else if ( name == "GLN" || name == "GLU" )
  {
    tail = pres->GetAtom ( "CG", NULL, aloc );
    head = pres->GetAtom ( "CD", NULL, aloc );
    
    mCMMDBManager.Select (
      selHnd,   // selection handle
      STYPE_ATOM,  // select atoms
      0,           // any model
      "A",         // chains "A" and "B" only
      seqnum,"*", // this residue
      seqnum,"*", //   any insertion code
      "*",         // any residue name
      "!N,CA,C,O,CB,CG,CD", // C-alphas only
      "*",         // carbons only
      aloc,        // any alternative location indicator
      SKEY_NEW     // NEW-selection
    );  
  }
  else if ( name == "ARG" )
  {
    tail = pres->GetAtom ( "NE", NULL, aloc );
    head = pres->GetAtom ( "CZ", NULL, aloc );
    
    mCMMDBManager.Select (
      selHnd,   // selection handle
      STYPE_ATOM,  // select atoms
      0,           // any model
      "A",         // chains "A" and "B" only
      seqnum,"*", // this residue
      seqnum,"*", //   any insertion code
      "*",         // any residue name
      "!N,CA,C,O,CB,CG,CD,NE,CZ", // C-alphas only
      "*",         // carbons only
      aloc,        // any alternative location indicator
      SKEY_NEW     // NEW-selection
    );
  }
  mCMMDBManager.GetSelIndex ( selHnd, atms, natms );

  VectorRotation ( atms, natms, M_PI, head->x - tail->x,
                                      head->y - tail->y,
                                      head->z - tail->z,
                                      tail->x,
                                      tail->y,
                                      tail->z
                                      );
  mCMMDBManager.DeleteSelection ( selHnd );
  mCMMDBManager.FinishStructEdit ( );
}

double GetRMSD ( PPCAtom atoms, PPCAtom atoms_ref, const int sz )
{
  double d  = 0;
  int idx = -1;
  for ( int i=0; i<sz; ++i )
  {
    for ( int ii=0; ii<sz; ++ii )
      if ( atoms[ii]->CheckID ( atoms_ref[i]->name, NULL, NULL ) )
      {
        idx = ii;
	ii=sz;
      }
    
    if ( idx < 0 )
    {
      std::cerr << "Atom with same name not found." << std::endl;
      return 1001;
    }
    d += ( atoms[idx]->x - atoms_ref[i]->x )*( atoms[idx]->x - atoms_ref[i]->x ) +
         ( atoms[idx]->y - atoms_ref[i]->y )*( atoms[idx]->y - atoms_ref[i]->y ) +
         ( atoms[idx]->z - atoms_ref[i]->z )*( atoms[idx]->z - atoms_ref[i]->z );
  }
  return sqrt ( d/sz);
}
    

double Edist ( double frst[3], double scnd[3] )
{
  double d = 0;
  for ( int i=0; i<3; i++ )
    d += ( frst[i] - scnd[i] )*( frst[i] - scnd[i] );
  return sqrt ( d );
}

int main ( int argc, char* argv[] )
{
  if ( argc != 4 )
  {
    std::cout << "Usage: postrmsd pdbin1 pdbin2 chainID." << std::endl;
    exit(1);
  }
      
  PCChain pChnTbl[2];
  CMMDBManager myCMMDBManager[2];

  for ( int i=0; i<2; i++ )
  {
    int RC = myCMMDBManager[i].ReadPDBASCII ( const_cast<char*>( argv[i+1] ) );
    pChnTbl[i] = myCMMDBManager[i].GetChain ( 1, argv[3] );
    if ( !pChnTbl[i] )
    {
      std::cout << "No Chain " << argv[3] << " in file." << std::endl;
      exit ( 2 );
    }
    myCMMDBManager[i].PDBCleanup ( PDBCLEAN_ATNAME | PDBCLEAN_INDEX | PDBCLEAN_ALTCODE_STRONG );
    myCMMDBManager[i].FinishStructEdit ( );

  }
  
  for ( int i=0; i<pChnTbl[0]->GetNumberOfResidues(); i++ )
  {
    int nAltLocs_ref, alflag_ref, nAltLocs, alflag;
    PAltLoc pALoc, pALoc_ref;
    rvector  occupancy, occupancy_ref;
    occupancy = occupancy_ref = NULL;
    pALoc = pALoc_ref = NULL;    
       
    PCResidue pRsd_ref = pChnTbl[0]->GetResidue(i);
    
    std::string name = pRsd_ref->GetResName ( );
    const bool is2Rotate = ( name == "ASN" || name == "ASP" || name == "TYR" || name == "HIS" || name == "PHE" || 
                       name == "GLN" || name == "GLU" || name == "ARG" );
	
    pRsd_ref->GetAltLocations ( nAltLocs_ref,pALoc_ref,occupancy_ref,alflag_ref );
    
    PCResidue pRsd = pChnTbl[1]->GetResidue(pRsd_ref->GetSeqNum(),"");
/*    
    char S[128];
    if ( pRsdTbl[i]->GetAtom ( "CA", "*", char(0) ) )
      std::cout << pRsdTbl[i]->GetAtom ( "CA", "*", char(0) )->GetAtomID ( S ) << std::endl;
    
    AltLoc aL;
    aL[0] = ' ';
    aL[1] = char(0);
    if ( !pRsd->GetAtom ( "CA", "*", aL ) )
      continue;     
    
    std::vector<std::pair<float,int> > sort_occ_idx;
    
    for ( int ii=0; ii<nAltLocs; ++ii )
      sort_occ_idx.push_back ( std::make_pair ( occupancy[ii], ii) );
*/     
    std::cout << pRsd_ref->GetSeqNum() << pRsd_ref->GetResName () << " ";
    for ( int ii=0; ii<nAltLocs_ref; ++ii )
    {
      PPCAtom      SelAtom_ref, SelAtom;
      int          selHnd_ref, nSelAtoms_ref, selHnd,nSelAtoms;
      double rmsd, rmsd_min;
      int cntr_min;
      rmsd_min = 1000;
      
      AltLoc aLoc;
      aLoc[0] = pALoc_ref[ii][0];
      aLoc[1] = char(0);
	
      selHnd_ref = myCMMDBManager[0].NewSelection();
      
      myCMMDBManager[0].Select (
        selHnd_ref,   // selection handle
        STYPE_ATOM,  // select atoms
        0,           // any model
        "A",         // chains "A" and "B" only
        pRsd_ref->GetSeqNum (),"*", // this residue
        pRsd_ref->GetSeqNum (),"*", //   any insertion code
        "*",         // any residue name
        "!N,CA,C,O", // C-alphas only
        "*",         // carbons only
        aLoc,        // any alternative location indicator
        SKEY_NEW     // NEW-selection
      );

     myCMMDBManager[0].GetSelIndex ( selHnd_ref,SelAtom_ref,nSelAtoms_ref );
     
     if ( nSelAtoms_ref == 0 )
       continue;

     pRsd->GetAltLocations ( nAltLocs,pALoc,occupancy,alflag );
     
     for ( int iii=0; iii<nAltLocs; ++iii )
     {
       AltLoc aLoc2;
       aLoc2[0] = pALoc[iii][0];
       aLoc2[1] = char(0);
       
       if ( aLoc2[0] == ' ' && nAltLocs > 1 )
         continue;
       
       const int N = is2Rotate ? 2 : 1;
       
       for ( int v=0; v<N; v++ )
       {
         selHnd = myCMMDBManager[1].NewSelection();
         myCMMDBManager[1].Select (
           selHnd,   // selection handle
           STYPE_ATOM,  // select atoms
           0,           // any model
           "A",         // chains "A" and "B" only
           pRsd->GetSeqNum (),"*", // this residue
           pRsd->GetSeqNum (),"*", //   any insertion code
           "*",         // any residue name
           "!N,CA,C,O", // C-alphas only
           "*",         // carbons only
           aLoc2,        // any alternative location indicator
           SKEY_NEW     // NEW-selection
         );
       
         myCMMDBManager[1].GetSelIndex ( selHnd,SelAtom,nSelAtoms );
       
         if ( nSelAtoms != nSelAtoms_ref )
         {
           std::cout << "WARNING: Number of Atoms differ. ";
	   break;
         }
	 
         int min_size = nSelAtoms <= nSelAtoms_ref ? nSelAtoms : nSelAtoms_ref;
       
         if ( min_size > 1 || (min_size == 1 && name == "ALA" ))
           rmsd = GetRMSD ( SelAtom, SelAtom_ref, min_size );
         else
	 {
           std::cout << "WARNING: Number of Atoms equals 0. ";
	   continue;
	 }
       
         if ( rmsd < rmsd_min )
         {
           rmsd_min = rmsd;
	   cntr_min = iii;
         }
         myCMMDBManager[1].DeleteSelection ( selHnd );
	 if ( is2Rotate )
	   RotateFinalChi ( myCMMDBManager[1], pRsd->GetSeqNum (), aLoc2 );
       }
    }
     std::cout << aLoc[0] << "<->" << pALoc[cntr_min][0] << " " << rmsd_min << " ";
    
  }
  std::cout << std::endl;
 }
  return 0;
}

