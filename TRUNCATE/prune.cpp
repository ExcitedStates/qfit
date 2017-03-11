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
  return ( pAtm->CheckID ( "CA", NULL, NULL ) || pAtm->CheckID ( "N", NULL, NULL ) || pAtm->CheckID ( "C", NULL, NULL ) || pAtm->CheckID ( "O", NULL, NULL ) 
        || pAtm->CheckID ( "CB", NULL, NULL ));
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

std::string resname, seqnum;
float cc_all, md_all, md_main, md_side;
int main ( int argc, char* argv[] )
{
  if ( argc != 4 )
  {
    std::cout << "Usage: prune pdbin pdbout chainID." << std::endl;
    exit(1);
  }
  
   
  std::ifstream cc("cc.log", std::ios::in);
  std::vector<int> fragids;

  if ( !cc.is_open() )
    std::cout << "Unable to open file cc.log." << std::endl;
  
  std::string cc_content;
  while ( cc >> cc_content && cc_content != "SIDE" );

  cc.ignore ( 256, '\n');
  
  std::vector<std::pair<std::string,std::vector<float> > > seqnum_md_side;
  
  while ( cc.peek() != '\n' )
  {
    cc >> resname >> seqnum >> cc_all >> md_all >> md_main >> md_side;
    std::string::iterator sit;
    sit = seqnum.end()-1;
    seqnum.erase ( sit );
    if ( resname == "GLY" || resname == "ALA" || resname == "PRO" )
    {
      cc.ignore ( 256, '\n');
      continue;
    }
    if ( seqnum_md_side.size() == 0 || seqnum_md_side.back().first != seqnum )
    {
      std::vector<float> vcc;
      seqnum_md_side.push_back ( std::make_pair ( seqnum, vcc ) );
      seqnum_md_side.back().second.push_back ( md_side );
    }
    else
      seqnum_md_side.back().second.push_back ( md_side );
    cc.ignore ( 256, '\n');
  }
  
  for ( int i=0; i<seqnum_md_side.size(); ++i )
  {
    std::cout << seqnum_md_side[i].first << " ";
    for ( int j=0; j<seqnum_md_side[i].second.size(); ++j )
       std::cout << seqnum_md_side[i].second[j] << " ";
    std::cout <<std::endl;
  }
   
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
//  myCMMDBManager.PDBCleanup ( PDBCLEAN_ATNAME | PDBCLEAN_INDEX | PDBCLEAN_ALTCODE_STRONG );
//  myCMMDBManager.FinishStructEdit ( );

  if ( !pChn )
  {
    std::cout << "No Chain " << argv[3] << " in file." << std::endl;
    exit ( 2 );
  }
  char S[128];
  for ( int i=0; i<seqnum_md_side.size(); ++i )
  {
   //if ( atoi ( seqnum_md_side[i].first.c_str() ) != 314 )continue;
   PPCAtom atmtbl;
   int natms;
   PCResidue pRsd = pChn->GetResidue ( atoi ( seqnum_md_side[i].first.c_str() ), "" );
   pRsd->GetAtomTable ( atmtbl, natms );

   std::vector<int> todel;
   int ndel_side = 0;
    for ( int j=0; j<seqnum_md_side[i].second.size(); ++j )
      if ( seqnum_md_side[i].second[j] < 0.75 )
      { 
        ndel_side++;
        if ( seqnum_md_side[i].second.size() == 1 )
          for ( int k=0; k<natms; k++ )
	  {
	    if ( !isMC ( atmtbl[k] ) )
	    {
	      std::cout << "Atom deleted: " << atmtbl[k]->GetAtomID ( S ) << std::endl;
	      todel.push_back ( k );
	      //pRsd->DeleteAtom ( k );
	    }
	  }
        else
          for ( int k=0; k<natms; k++ )
          {
            if ( !isMC ( atmtbl[k] ) && atmtbl[k]->altLoc[0] == 'A'+j )
	    {
	      std::cout << "Atom deleted: " << atmtbl[k]->GetAtomID ( S ) << std::endl;
	      todel.push_back ( k );
	    }
	  }
       }
     std::vector<char> vaLoc;
     for ( int k=0; k<todel.size(); ++k )
     {
       vaLoc.push_back ( pRsd->GetAtom(todel[k])->altLoc[0] );
       if ( todel[k] )
         pRsd->DeleteAtom ( todel[k] );
     }
     std::vector<char>::iterator new_end=unique(vaLoc.begin(), vaLoc.end());
     vaLoc.erase(new_end, vaLoc.end());
     std::cout << "size: " << vaLoc.size() << " ndel " << ndel_side << " seqnum_md_side[i].second.size() " << seqnum_md_side[i].second.size() << std::endl;
     if ( pRsd->GetAtom ( "CA", NULL, "" ) /*just one mc conf*/ )
     {
       if ( seqnum_md_side[i].second.size() - ndel_side == 0 /*all sc deleted*/ )
       {
         for ( int j=0; j<vaLoc.size()-1; ++j )
           pRsd->DeleteAtom ( "CB", NULL, (AltLoc){vaLoc[j],char(0)} );
	 PCAtom pAtm = pRsd->GetAtom ( "CB", NULL, (AltLoc){vaLoc.back(),char(0)} );
	 if ( !pAtm )
	   std::cout << "No Atom with altLoc " << vaLoc.back() << " in residue " << seqnum_md_side[i].second.size() << std::endl;
	 pAtm->occupancy = 1.00;
	 pAtm->altLoc[0] = char(0);
       }
       else
       {
         for ( int j=0; j<vaLoc.size(); ++j )
           pRsd->DeleteAtom ( "CB", NULL, (AltLoc){vaLoc[j],char(0)} );
	 myCMMDBManager.FinishStructEdit ( );
	 if ( seqnum_md_side[i].second.size() - ndel_side == 1 ) /*just one left*/
	 {
	   PPCAtom atbl;
	   int na;
	   pRsd->GetAtomTable ( atbl, na );
	   for ( int m=0; m<na; m++ )
	     atbl[m]->altLoc[0] = char(0);
	 }
       }
     }
   }
 
 myCMMDBManager.PDBCleanup ( PDBCLEAN_ATNAME | PDBCLEAN_INDEX | PDBCLEAN_ALTCODE_STRONG );
 myCMMDBManager.FinishStructEdit ( );
 //std::cout << "Total atoms deleted: " << tot_atms_deleted << std::endl;
 myCMMDBManager.WritePDBASCII ( argv[2] );

  return 0;
}

