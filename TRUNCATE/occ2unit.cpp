#include <iostream>
#include <cstring>
#include <fstream>
#include <stdlib.h>
#include <limits>
#include <cmath>
#include <vector>
#include <mmdb/mmdb_manager.h>

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
  if ( argc != 4 )
  {
    std::cout << "Usage: idmulti pdbin pdbout chainID." << std::endl;
    exit(1);
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
  myCMMDBManager.PDBCleanup ( PDBCLEAN_ATNAME | PDBCLEAN_INDEX | PDBCLEAN_ALTCODE_STRONG );
  myCMMDBManager.FinishStructEdit ( );

  if ( !pChn )
  {
    std::cout << "No Chain " << argv[3] << " in file." << std::endl;
    exit ( 2 );
  }

  pChn->GetResidueTable ( pRsdTbl, nRsds );
  
  for ( int i=0; i<nRsds; i++ )
  {
    if ( strncmp ( pRsdTbl[i]->GetResName ( ), "GLY", 3 ) == 0 || 
       ( strncmp ( pRsdTbl[i]->GetResName ( ), "MSE", 3 ) != 0 && !pRsdTbl[i]->isAminoacid() ) ) continue;
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
    //std::cout << " Processing Sequence Number :" << pRsdTbl[i]->GetSeqNum ( ) << std::endl;
       
    pRsdTbl[i]->GetAltLocations ( nAltLocs,pALoc,occupancy,alflag );
    
/*    for ( int ii=0; ii<nAltLocs; ++ii )
    {
      //std::cout << pALoc[ii][0] << std::endl;
      AltLoc aLoc;
      aLoc[0] = pALoc[ii][0];
      aLoc[1] = char(0);
      PCAtom pAtm = pRsdTbl[i]->GetAtom ( "CA", "*", aLoc );
      if ( pAtm )
      {
        if ( !CApos )
        {
	  CApos = new double[3];
	  CApos[0]= pAtm->x; 
	  CApos[1]= pAtm->y;
	  CApos[2]= pAtm->z;
        }
	else 
	{
	  double pos[3] = {pAtm->x, pAtm->y, pAtm->z };
	  d = Edist ( CApos, pos );
	  //std::cout << "d = " << d << std::endl;
	}
	tot_occ += pAtm->occupancy;
      }
      else
      {
        //std::cout << "Atom " << pALoc[ii][0] << "CA not found." << std::endl;
	//exit ( 1 );
      }
      if ( d > TOL )
        single_pos = false;
    }
     */     
    if (pALoc)  delete[] pALoc;
    FreeVectorMemory ( occupancy, 0 );
    delete [] CApos;
    CApos = NULL;
    
    pRsdTbl[i]->GetAtomTable ( pAtmTbl, nAtms );
    
    if ( nAltLocs == 1 ) /* Whole residue just 1 conformation */
    {
      for ( int ii=0; ii<nAtms; ii++ )
      {
        pAtmTbl[ii]->altLoc[0] = char(0);
	pAtmTbl[ii]->occupancy = 1.0;
      }
    }  
    else  
    {
      std::vector <std::pair <std::string,double> > vname_occ;
      double tocc = 0.0;
      for ( int ii=0; ii<nAtms; ii++ )
        if ( pAtmTbl[ii]->CheckID ( "CB", NULL, NULL ) )
          tocc += pAtmTbl[ii]->occupancy;

      for ( int ii=0; ii<nAtms; ii++ )
        if ( pAtmTbl[ii]->occupancy < 1.0 || pAtmTbl[ii]->occupancy > 1.0) /* exclude main-chain positions coincide */
          pAtmTbl[ii]->occupancy /= tocc;
      
      AltLoc aLoc;
      tocc = 0.0;
      for ( int ii=0; ii<nAtms; ii++ )
        if ( pAtmTbl[ii]->CheckID ( "CB", NULL, NULL ) )
	{
          tocc += floorf(pAtmTbl[ii]->occupancy*100+0.5)/100;
	  aLoc[0] = pAtmTbl[ii]->altLoc[0];
	}

      if ( tocc < 1.0 || tocc > 1.0 )
        for ( int ii=0; ii<nAtms; ii++ )
          if ( pAtmTbl[ii]->altLoc[0] == aLoc[0] ) /* Adjust only one alternate conformation` */
            pAtmTbl[ii]->occupancy += (1.0 - tocc);

      myCMMDBManager.FinishStructEdit ( );
     /* TEST : Output occupancies for all atoms 
      vname_occ.clear();
      for ( int ii=0; ii<nAtms; ii++ ) {
	bool found = false;
	if ( vname_occ.size() == 0 ) {
          vname_occ.push_back ( make_pair ( (std::string)pAtmTbl[ii]->GetAtomName(), floorf(pAtmTbl[ii]->occupancy*100+0.5)/100 ) );
	  continue;
	}  
        for ( int iii=0; iii<vname_occ.size(); ++iii) 
          if ( (std::string)pAtmTbl[ii]->GetAtomName() == vname_occ[iii].first ) {
            vname_occ[iii].second += floorf(pAtmTbl[ii]->occupancy*100+0.5)/100;;
	    found = true;
	  }
        if ( ! found )
          vname_occ.push_back ( make_pair ( (std::string)pAtmTbl[ii]->GetAtomName(), floorf(pAtmTbl[ii]->occupancy*100+0.5)/100 ) );
      }
      for ( int iii=0; iii<vname_occ.size(); ++iii)
        std::cout << vname_occ[iii].first << " " << vname_occ[iii].second << std::endl;
	*/
    }
  }
  
  myCMMDBManager.PDBCleanup ( PDBCLEAN_ATNAME | PDBCLEAN_INDEX | PDBCLEAN_ALTCODE_STRONG );
  myCMMDBManager.FinishStructEdit ( );

  myCMMDBManager.WritePDBASCII ( argv[2] );

  return 0;
}

