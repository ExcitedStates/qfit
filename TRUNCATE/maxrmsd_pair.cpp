#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits>
#include <cmath>
#include <vector>
#include <map>
#include <mmdb/mmdb_manager.h>

template<class T, class U> struct ComparePair1stDsc 
{
  bool operator()( const std::pair<T,U>& a, const std::pair<T,U>& b ) const 
  {
    return a.first > b.first; 
  }
};

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
  
  if ( seq_stop - seq_start > 7 ) 
    success = true;
  return success;
}

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

int main ( int argc, char* argv[] )
{
  if ( argc != 3 )
  {
    std::cout << "Usage: maxrmsd pdbin1 pdbin2." << std::endl;
    exit(1);
  }
  
  PCChain pChn;
  PPCChain pChnTbl;
  PPCAtom pAtmTbl;
  int nAtms;
  PPCResidue pRsdTbl;
  int nRsds, nChns;
  const double TOL = 0.01;
  std::vector<int> multi;

  CMMDBManager myCMMDBManager[2];

  //std::string pdbin1 = argv[1+npdb];
  int RC = myCMMDBManager[0].ReadPDBASCII ( const_cast<char*>( argv[1] ) );
  RC = myCMMDBManager[1].ReadPDBASCII ( const_cast<char*>( argv[2] ) );

  std::string pdbID = std::string(argv[1]).substr(0,4);

//  std::cout << pdbID << std::endl;
  std::map<int,int> buried;
  std::ifstream fin;
  fin.open("henry_buried.txt");
  char pdbid[4];
  char buf[100];
  fin.getline(buf,99);
  int resid, bur;
  while(!fin.eof()){
    int r = sscanf(buf, "%s %d %d", pdbid, &resid, &bur);
    if ( (std::string)pdbid == pdbID )
      buried[resid] = bur;
    //std::cout <<  pdbid << " " << resid << " " << bur << std::endl;
    fin.getline(buf,99);
  }
  
  myCMMDBManager[0].GetChainTable ( 1, pChnTbl, nChns );

  std::map<int,float> vmaxrmsd[2];
  int nCA[2];
  nCA[0] = nCA[1] = 0;

  for ( int iChn=0; iChn<1; iChn++ )
  {
    int s_start, s_stop;
    std::string resid;
    int start = 10000;
    int stop = -10000;
    if ( !isAAGEQ7Chain ( pChnTbl[iChn], s_start, s_stop ) )
      continue;
 
    start = (s_start < start ) ? s_start : start;
    stop = (s_stop > stop ) ? s_stop : stop;

    for ( int npdb=0; npdb<2; npdb++ )
    {
      if ( npdb == 0 )
        pChn = pChnTbl[iChn];
      else
        pChn = myCMMDBManager[1].GetChain (1, pChnTbl[iChn]->GetChainID() );

      std::cerr << pChnTbl[iChn]->GetChainID() << std::endl;

      if ( !pChn )
      {
        std::cout << "No Chain " << argv[2] << " in file." << std::endl;
        exit ( 2 );
      }

      pChn->GetResidueTable ( pRsdTbl, nRsds );
  
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
        double rmsd, maxrmsd;
        maxrmsd = 0.0;
    
 /*   if  ( std::find ( fragids.begin(), fragids.end(), pRsdTbl[i]->GetSeqNum ( ) ) != fragids.end() )
    {
      std::cout << "Residue " << pRsdTbl[i]->GetSeqNum ( ) << "is part of of a frag. Continuing." << std::endl;
      continue;
    }
    */
    //std::cout << " Processing Sequence Number :" << pRsdTbl[i]->GetSeqNum ( ) << std::endl;
    
	/* MAIN CHAIN ONLY */
    
    if ( !pRsdTbl[i]->GetAtom ( "CA", "C", "" ) )
    {
//      std::cout << pRsdTbl[i]->GetSeqNum() << std::endl;
      nCA[npdb]++;
    }

    if ( npdb == 0 )
      if ( buried.find (  pRsdTbl[i]->GetSeqNum() ) != buried.end() )
        if ( buried[ pRsdTbl[i]->GetSeqNum()] > 24 )
          continue;
        
   if ( npdb == 1 ) 
   {
     PCResidue pres = pChnTbl[iChn]->GetResidue ( pRsdTbl[i]->GetSeqNum(), "" );
     
     if ( pres )
       if ( (std::string)pres->GetResName() != (std::string)pRsdTbl[i]->GetResName() )
       {
	std::cerr << pres->GetSeqNum() << pres->GetResName() << " " << (std::string)pRsdTbl[i]->GetResName() << std::endl;
	continue;
       }
   }
       /* if (npdb == 0 )
          resid = pRsdTbl[i]->GetResName();
	else
          if ( resid != (std::string)pRsdTbl[i]->GetResName() )
	  {
	    std::cerr << resid << " " << (std::string)pRsdTbl[i]->GetResName() << std::endl;
            continue;
	  }*/
       
	if ( (std::string)pRsdTbl[i]->GetResName() == "ALA" || (std::string)pRsdTbl[i]->GetResName() == "GLY" || (std::string)pRsdTbl[i]->GetResName() == "PRO" || !pRsdTbl[i]->isAminoacid())
	  continue;
        
	pRsdTbl[i]->GetAltLocations ( nAltLocs,pALoc,occupancy,alflag );

/*    
    char S[128];
    if ( pRsdTbl[i]->GetAtom ( "CA", "*", char(0) ) )
      std::cout << pRsdTbl[i]->GetAtom ( "CA", "*", char(0) )->GetAtomID ( S ) << std::endl;
*/    
        AltLoc aL;
        aL[0] = ' ';
        aL[1] = char(0);
 
        std::vector<std::pair<float,int> > sort_occ_idx;
    
        //std::cout <<  pRsdTbl[i]->GetSeqNum () << " ";
        if ( nAltLocs <= 1 )
        {
          //std::cout << "0" << std::endl;
          vmaxrmsd[npdb][pRsdTbl[i]->GetSeqNum ()] = maxrmsd;
          continue;
        }
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

         refselHnd = myCMMDBManager[npdb].NewSelection();
         myCMMDBManager[npdb].Select (
           refselHnd,      // selection handle
           STYPE_ATOM,  // select atoms
           0,           // any model
           pChn->GetChainID(),         // chains "A" and "B" only
           pRsdTbl[i]->GetSeqNum (),"*", // this residue
           pRsdTbl[i]->GetSeqNum (),"*", //   any insertion code
           "*",         // any residue name
           "!N,CA,C,O",        // C-alphas only
           "*",         // carbons only
           aLoc,         // any alternative location indicator
           SKEY_NEW     // NEW-selection
              );
	 
           myCMMDBManager[npdb].GetSelIndex ( refselHnd,refSelAtom,refnSelAtoms );
           if ( refnSelAtoms == 0 )
             continue;

         for ( int iv=ii+1; iv<nAltLocs; ++iv )
         {
           AltLoc aLoc2;
           aLoc2[0] = pALoc[iv][0];
           aLoc2[1] = char(0);

           const int nSC = refnSelAtoms;
           selHnd = myCMMDBManager[npdb].NewSelection();
           myCMMDBManager[npdb].Select (
             selHnd,   // selection handle
             STYPE_ATOM,  // select atoms
             0,           // any model
             pChn->GetChainID(),         // chains "A" and "B" only
             pRsdTbl[i]->GetSeqNum (),"*", // this residue
             pRsdTbl[i]->GetSeqNum (),"*", //   any insertion code
             "*",         // any residue name
             "!N,CA,C,O", // C-alphas only
             "*",         // carbons only
             aLoc2,        // any alternative location indicator
             SKEY_NEW     // NEW-selection
           );
       
          myCMMDBManager[npdb].GetSelIndex ( selHnd,SelAtom,nSelAtoms );
          int n = 0;
          for ( int iii=0; iii<nSelAtoms; iii++ )
            if ( SelAtom[iii]->altLoc[0] == aLoc2[0] )
	      n++;
     
          if ( n == nSC != 0)
          {
            rmsd = GetRMSD ( refSelAtom, SelAtom, n );
            //std::cout << aLoc << aLoc2 << " " << rmsd << " ";
	    maxrmsd = (rmsd > maxrmsd ) ? rmsd : maxrmsd;
          }
        }
      }
      //std::cout << maxrmsd << std::endl;
      vmaxrmsd[npdb][pRsdTbl[i]->GetSeqNum ()] = maxrmsd;
    }
   }
    std::map<int,float>::iterator it1, it2;
    //std::cout << "************" << std::endl;
   for ( int s = start; s <= stop; s++ )
   {
      it1 = vmaxrmsd[0].find ( s );
      it2 = vmaxrmsd[1].find ( s );
      if ( it1 != vmaxrmsd[0].end() && it2 != vmaxrmsd[1].end() )
      {
        std::cout << s << " " << vmaxrmsd[0][s] << " " << vmaxrmsd[1][s] << " " << vmaxrmsd[0][s] - vmaxrmsd[1][s] << std::endl;
      }
   }  
    //std::cout << "************" << std::endl;
 }
  std::ofstream outfile;
  outfile.open ( "nCA.txt", std::ios::out | std::ios::app );
  outfile << nCA[0] << " " << nCA[1] << std::endl;
  return 0;
}

