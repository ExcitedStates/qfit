#include <iostream>
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
  if ( argc != 5)
  {
    std::cout << "Usage: idmulti tol pdbin pdbout chainID" << std::endl;
    exit(1);
  }
  
  PCChain pChn;
  PPCAtom pAtmTbl;
  int nAtms;
  PPCResidue pRsdTbl;
  int nRsds;
  const double TOL = atof (argv[1]);
  std::string split_atm;
  std::vector<int> multi;

  CMMDBManager myCMMDBManager;

  std::string inpdb = argv[2];
  int RC = myCMMDBManager.ReadPDBASCII ( const_cast<char*>( inpdb.c_str() ) );

  pChn = myCMMDBManager.GetChain ( 1, argv[4] );
  myCMMDBManager.PDBCleanup ( PDBCLEAN_ATNAME | PDBCLEAN_INDEX | PDBCLEAN_ALTCODE_STRONG );
  myCMMDBManager.FinishStructEdit ( );

  if ( !pChn )
  {
    std::cout << "No Chain " << argv[4] << " in file." << std::endl;
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
    CApos = NULL;
    double* Opos;
    Opos = NULL;
    bool CA_single_pos = true;
    double CA_d = 0.0;
    bool O_single_pos = true;
    double O_d = 0.0;
    double tot_occ = 0.0;
    //std::cout << " Processing Sequence Number :" << pRsdTbl[i]->GetSeqNum ( ) << std::endl;
       
    pRsdTbl[i]->GetAltLocations ( nAltLocs,pALoc,occupancy,alflag );
    
    for ( int ii=0; ii<nAltLocs; ++ii )
    {
      //std::cout << pALoc[ii][0] << std::endl;
      AltLoc aLoc;
      aLoc[0] = pALoc[ii][0];
      aLoc[1] = char(0);
      
      ////////////////////////////////////////////////////////////////////////
      // CA (use this to tally occupancy)
      split_atm = "CA";
      PCAtom pAtm = pRsdTbl[i]->GetAtom ( split_atm.c_str(), "*", aLoc );
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
      	  CA_d = Edist ( CApos, pos );
      	  //std::cout << "d = " << d << std::endl;
      	}
      	tot_occ += pAtm->occupancy;
      }
      else
      {
        //std::cout << "Atom " << pALoc[ii][0] << "CA not found." << std::endl;
	      //exit ( 1 );
      }
      if ( CA_d > TOL )
        CA_single_pos = false;

      ////////////////////////////////////////////////////////////////////////
      // O
      split_atm = "O";
      pAtm = pRsdTbl[i]->GetAtom ( split_atm.c_str(), "*", aLoc );
      if ( pAtm )
      {
        if ( !Opos )
        {
          Opos = new double[3];
          Opos[0]= pAtm->x; 
          Opos[1]= pAtm->y;
          Opos[2]= pAtm->z;
        }
        else 
        {
          double pos[3] = {pAtm->x, pAtm->y, pAtm->z };
          O_d = Edist ( Opos, pos );
          //std::cout << "d = " << d << std::endl;
        }
      }
      else
      {
        //std::cout << "Atom " << pALoc[ii][0] << "CA not found." << std::endl;
        //exit ( 1 );
      }
      if ( O_d > TOL )
        O_single_pos = false;
    }
       
    if (pALoc)  delete[] pALoc;
    FreeVectorMemory ( occupancy, 0 );
    delete [] CApos;
    CApos = NULL;
    delete [] Opos;
    Opos = NULL;
    
    pRsdTbl[i]->GetAtomTable ( pAtmTbl, nAtms );
    
    if ( nAltLocs == 1 ) /* Whole residue just 1 conformation */
    {
      for ( int ii=0; ii<nAtms; ii++ )
      {
        pAtmTbl[ii]->altLoc[0] = char(0);
	      pAtmTbl[ii]->occupancy = 1.0;
      }
    }
    else if ( CA_single_pos && O_single_pos ) /* main-chain positions coincide */
    {
      for ( int ii=0; ii<nAtms; ii++ )
      {
        if ( isMC ( pAtmTbl[ii] ) )
      	{
      	  if ( pAtmTbl[ii]->altLoc[0] > 'A' )
            pRsdTbl[i]->DeleteAtom ( ii );
          else
          {
            pAtmTbl[ii]->altLoc[0] = char(0);
            pAtmTbl[ii]->occupancy = 1.0; //tot_occ;
          }
        }
      }
    }
    else 
      multi.push_back ( pRsdTbl[i]->GetSeqNum ( ) );
  }
  
  bool cons;
  for ( int i=1; i<multi.size(); ++i )
  {
    if ( multi[i-1] + 1 == multi[i] )
    {
      std::cout << multi[i-1] << " ";
      cons = true;
      if ( i == multi.size()-1 )
        std::cout << multi[i] << std::endl;
    }
    else if ( cons )
    {
      std::cout << multi[i-1] << std::endl;
      cons = false;
    }
  }
  //std::cout << std::endl;
  myCMMDBManager.PDBCleanup ( PDBCLEAN_ATNAME | PDBCLEAN_INDEX | PDBCLEAN_ALTCODE_STRONG );
  myCMMDBManager.FinishStructEdit ( );

  myCMMDBManager.WritePDBASCII ( argv[3] );

  return 0;
}
