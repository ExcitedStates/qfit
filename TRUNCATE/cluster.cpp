#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits>
#include <cmath>
#include <vector>
#include <algorithm>
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

double Edist ( double frst[3], double scnd[3] )
{
  double d = 0;
  for ( int i=0; i<3; i++ )
    d += ( frst[i] - scnd[i] )*( frst[i] - scnd[i] );
  return sqrt ( d );
}

int main ( int argc, char* argv[] )
{
  if ( argc != 5 )
  {
    std::cout << "Usage: cluster pdbin pdbout chainID distance." << std::endl;
    exit(1);
  }
  
   
  std::ifstream frags("fragids.txt", std::ios::in);
  std::vector<int> fragids;

  if ( !frags.is_open() )
    std::cout << "Unable to open file fragids.txt." << std::endl;
  
  int num;
  while ( frags >> num )
    fragids.push_back ( num );
    
  PCChain pChn;
  PPCAtom pAtmTbl;
  int nAtms;
  PPCResidue pRsdTbl;
  int nRsds;
  const double TOL = 0.01;
  std::vector<int> multi;

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

  pChn->GetResidueTable ( pRsdTbl, nRsds );
  
  int tot_atms_deleted = 0;
  
  for ( int i=0; i<nRsds; i++ )
  {
    int nAltLocs, alflag;
    PAltLoc pALoc;
    rvector  occupancy;
    occupancy = NULL;
    pALoc = NULL;
    double* CApos;
    CApos= NULL;
    bool single_pos = true;
    double d = 0.0;
    double tot_occ = 0.0;
    
    if  ( std::find ( fragids.begin(), fragids.end(), pRsdTbl[i]->GetSeqNum ( ) ) != fragids.end() )
    {
      std::cout << "Residue " << pRsdTbl[i]->GetSeqNum ( ) << "is part of of a frag. Continuing." << std::endl;
      continue;
    }
    //std::cout << " Processing Sequence Number :" << pRsdTbl[i]->GetSeqNum ( ) << std::endl;
       
    pRsdTbl[i]->GetAltLocations ( nAltLocs,pALoc,occupancy,alflag );
/*    
    char S[128];
    if ( pRsdTbl[i]->GetAtom ( "CA", "*", char(0) ) )
      std::cout << pRsdTbl[i]->GetAtom ( "CA", "*", char(0) )->GetAtomID ( S ) << std::endl;
*/    
    AltLoc aL;
    aL[0] = ' ';
    aL[1] = char(0);
    if ( !pRsdTbl[i]->GetAtom ( "CA", "*", aL ) )
      continue;     
    
    std::vector<std::pair<float,int> > sort_occ_idx;
    
    for ( int ii=0; ii<nAltLocs; ++ii )
      sort_occ_idx.push_back ( std::make_pair ( occupancy[ii], ii) );
     
    for ( int ii=0; ii<nAltLocs; ++ii )
    {
      PPCAtom      refSelAtom, SelAtom;
      int          refselHnd, refnSelAtoms, selHnd,nSelAtoms;
      
      AltLoc aLoc;
      aLoc[0] = pALoc[ii][0];
      aLoc[1] = char(0);
      
//      std::cout << aLoc[0] << std::endl;

      refselHnd = myCMMDBManager.NewSelection();
      myCMMDBManager.Select (
        refselHnd,      // selection handle
        STYPE_ATOM,  // select atoms
        0,           // any model
        argv[3],         // chains "A" and "B" only
        pRsdTbl[i]->GetSeqNum (),"*", // this residue
        pRsdTbl[i]->GetSeqNum (),"*", //   any insertion code
        "*",         // any residue name
        "!N,CA,C,O",        // C-alphas only
        "*",         // carbons only
        aLoc,         // any alternative location indicator
        SKEY_NEW     // NEW-selection
              );
	 
     myCMMDBManager.GetSelIndex ( refselHnd,refSelAtom,refnSelAtoms );
       if ( refnSelAtoms == 0 )
         continue;

/*     
     char S[128];
     for ( int ll=0; ll<refnSelAtoms; ++ll )
	std::cout << refSelAtom[ll]->GetAtomID ( S ) << std::endl;
*/
     const int nSC = refnSelAtoms;
  
     selHnd = myCMMDBManager.NewSelection();
     myCMMDBManager.SelectNeighbours (
       selHnd,        // selection handle
       STYPE_ATOM, // select residues
       refSelAtom,       // array of already selected atoms
       refnSelAtoms,     // number of already selected atoms
       0.0,atof(argv[4]),      // distance range 0 to 10 angstroms
       SKEY_OR        // OR-selection
             );
     
     myCMMDBManager.GetSelIndex ( selHnd,SelAtom,nSelAtoms );
     
/*
     std::cout << "Neighbors:" << std::endl;
     for ( int ll=0; ll<nSelAtoms; ++ll )
       std::cout << SelAtom[ll]->GetAtomID ( S ) << std::endl;
*/
     for ( int iv=ii+1; iv<nAltLocs; ++iv )
     {
       AltLoc aLoc2;
       aLoc2[0] = pALoc[iv][0];
       aLoc2[1] = char(0);

       int n = 0;
       for ( int iii=0; iii<nSelAtoms; iii++ )
         if ( SelAtom[iii]->altLoc[0] == aLoc2[0] )
	   n++;
       
       if ( n == nSC != 0)
       {
         double occ;
         for ( int ll=0; ll<nSC; ++ll )
	 {
	   PCAtom patm = pRsdTbl[i]->GetAtom ( refSelAtom[ll]->name, NULL, aLoc2 );
	   occ = patm->occupancy;
	   refSelAtom[ll]->occupancy += patm->occupancy;
	 }
	 int atms_deleted = pRsdTbl[i]->DeleteAtom ( NULL, NULL, aLoc2 );
	 tot_atms_deleted += atms_deleted;
         std::cout << pRsdTbl[i]->GetSeqNum () << pRsdTbl[i]->GetResName() << " " << aLoc2[0];
	 std::cout << " DELETED: " << atms_deleted << " atoms." << std::endl;
	 std::cout << "Occupancy " << occ << " added to " << refSelAtom[0]->altLoc[0] << std::endl;
	 std::cout << "***" << std::endl;
       }
     }
   }
    
   myCMMDBManager.FinishStructEdit ( );
   pRsdTbl[i]->GetAltLocations ( nAltLocs,pALoc,occupancy,alflag );
//   std::cout << "New nAltLocs: " << nAltLocs << std::endl;
   if ( nAltLocs == 2 && pALoc[0][0] == ' ' )
   {
      PPCAtom atmTbl;
      int natms;
      pRsdTbl[i]->GetAtomTable ( atmTbl, natms );
       for ( int iii=0; iii<natms; iii++ )
        atmTbl[iii]->altLoc[0] = char(0);
   }
 }
 
 myCMMDBManager.PDBCleanup ( PDBCLEAN_ATNAME | PDBCLEAN_INDEX | PDBCLEAN_ALTCODE_STRONG );
 myCMMDBManager.FinishStructEdit ( );
 std::cout << "Total atoms deleted: " << tot_atms_deleted << std::endl;
 myCMMDBManager.WritePDBASCII ( argv[2] );

  return 0;
}

