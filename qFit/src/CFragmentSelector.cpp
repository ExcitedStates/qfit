/*
    qFit: Multiconformer modeling by constrained fitting of rotamer occupancies
    Henry van den Bedem, Ankur Dhanik, Jean-Claude Latombe, Ashley Deacon. Acta Cryst. D65:1107–1117 (2009)
    e-mail: vdbedem@slac.stanford.edu

        Copyright (C) 2009-2012 Stanford University

	Permission is hereby granted, free of charge, to any person obtaining a copy of
	this software and associated documentation files (the "Software"), to deal in
	the Software without restriction, including without limitation the rights to
	use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
	of the Software, and to permit persons to whom the Software is furnished to do
	so, subject to the following conditions: 

	This entire text, including the above copyright notice and this permission notice
	shall be included in all copies or substantial portions of the Software. 

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
	OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
	IN THE SOFTWARE.

    
*/


#include "CFragmentSelector.h"

bool WritePDB ( const PCChain pChain, const std::string fn );

// just for the sample -- print the data sets
std::ostream&
operator<<(std::ostream& os, const Vi& vi)
{
  os << "(";
  std::copy(vi.begin(), vi.end(), std::ostream_iterator<int>(os, ", "));
  os << ")";
  return os;
}
std::ostream&
operator<<(std::ostream& os, const Vvi& vvi)
{
  os << "(\n";
  for(Vvi::const_iterator it = vvi.begin();
      it != vvi.end();
      it++) {
      os << "  " << *it << "\n";
  }
  os << ")";
  return os;
}

CFragmentSelector::CFragmentSelector()
{
}

CFragmentSelector::~CFragmentSelector()
{
}

void CFragmentSelector::Init(
        std::string inpdb, std::string inmtz, const std::string ChnID, 
        const int lSeqNum, const int hSeqNum
        )
{
  clipper::Atom_list allatm_lst;
  clipper::ftype32 mean, sigsig;
  
  clipper::Xmap<float> xmap;
  clipper::HKL_info hkl_info;
  
//  clipper::CCP4MAPfile mapin;
//  mapin.open_read(inmap.c_str());
//  mapin.import_xmap(xmap);
//  mapin.close_read();

  srand ( 100 );
   
  int RC = myCMMDBManager.ReadPDBASCII ( const_cast<char*>( inpdb.c_str() ) );
  
  // get a list of all the atoms
  clipper::mmdb::PPCAtom psel;
  int hndl, nsel;
  hndl = myCMMDBManager.NewSelection();
  myCMMDBManager.SelectAtoms( hndl, 0, 0, SKEY_NEW );
  myCMMDBManager.SelectAtoms ( hndl,0,   
     (pstr)ChnID.c_str(),
     lSeqNum,"*",
     hSeqNum,"*",
     "*", "*", "*", "*", SKEY_XOR );

  myCMMDBManager.GetSelIndex( hndl, psel, nsel );
  allatm_lst = clipper::MMDBAtom_list ( psel, nsel );
  myCMMDBManager.DeleteSelection( hndl );
  
  std::cout << "Atom list size: " << allatm_lst.size() << std::endl;

  PScaleXMap mPScaleXMap;
  mPScaleXMap.DoScaleXMap ( xmap, hkl_info, inmtz, allatm_lst, mean, sigsig );
//  std::cout << "Resolution range: " << 1/sqrt(hkl_info.invresolsq_range ().max()) << "A - " << 1/sqrt(hkl_info.invresolsq_range ().min()) << "A" << std::endl;
  
//  clipper::CCP4MAPfile mapout;
//  mapout.open_write( "scaled.map");
//  mapout.export_xmap(xmap);
//  mapout.close_write();
  
  mCSolver.Init ( allatm_lst, xmap, hkl_info, mean, sqrt(sigsig) );
}


void CFragmentSelector::generate_paths_recursive(
        const std::vector<std::vector<char> > &a, 
        std::vector<char> &s, const int D 
        )
{
  if (D >= a.size()) {
//    for (int i=0; i<s.size(); ++i )
//      std::cout << s[i] << " ";
//    std::cout << std::endl;
    orderings.push_back(s);
  } else
    for ( int i=0; i<a[D].size(); ++i ) {
      s[D] = a[D][i];
      generate_paths_recursive ( a, s, D+1 );
    }
}

void CFragmentSelector::GetAltLocs (
        const int lowseq, const int highseq, const PCChain pChn, 
        std::vector<std::vector<char> >& altChars
        )
{
  for ( int seqnum=lowseq; seqnum<=highseq; seqnum++ ) {
    int nAltLocs, alflag;
    PAltLoc pALoc;
    rvector  occupancy;
    occupancy = NULL;
    pALoc = NULL;

    pChn->GetResidue(seqnum, "")->GetAltLocations(nAltLocs,pALoc,occupancy,alflag);

    for ( int k=0; k<nAltLocs; k++ )
      altChars[seqnum-lowseq].push_back( char('A' + k ) );

    if (pALoc)  delete[] pALoc;
    FreeVectorMemory ( occupancy, 0 );
  }
}


void CFragmentSelector::DoFindConfs(
        std::string inpdb, const std::string ChnID, 
        const int lSeqNum, const int hSeqNum
        )
{
  CMMDBManager mCMMDBManager;
//  PPCChain pChnTbl;
//  int nChns;
  PPCAtom pAtmTbl;
  int nAtms, lcount;
  char S[128];
  bool isSULFUR;

  int RC = mCMMDBManager.ReadPDBASCII ( const_cast<char*>( inpdb.c_str() ) );
  mCMMDBManager.FinishStructEdit ();
  
  if (RC) {
    //  An error was encountered. MMDB provides an error messenger
    //  function for easy error message printing.
    printf ( " ***** ERROR #%i READ:\n\n %s\n\n",RC,GetErrorDescription(RC) );
    //  Location of the error may be identified as precise as line
    //  number and the line itself. This information is now
    //  retrieved from MMDB input buffer:
    mCMMDBManager.GetInputBuffer ( S, lcount );
    if (lcount>=0) 
      printf ( "       LINE #%i:\n%s\n\n",lcount,S );
    else if (lcount==-1)
      printf ( "       CIF ITEM: %s\n\n",S );
    exit(1);
  }

  //mCMMDBManager.GetChainTable ( 1, pChnTbl, nChns );
  PCChain pMMDBChn = mCMMDBManager.GetChain(1, ChnID.c_str());
  if (!pMMDBChn->isAminoacidChain()) {
    std::cout << "No Chain " << ChnID.c_str() << std::endl; /*contains at least one AA*/
    exit(1);
  }
  std::vector<int> mse;
  PPCResidue pRsdTbl;
  int nRsds;
  pMMDBChn->GetResidueTable(pRsdTbl, nRsds);
  
  // Take care of MSE residues
  for ( int ii=lSeqNum; ii<=hSeqNum; ++ii ) 
    isMSE.push_back(strncmp(pMMDBChn->GetResidue(ii, "")->GetResName(), "MSE", 3) == 0);
  
  for ( int j=0; j<nRsds; j++ )
    if ( strncmp ( pRsdTbl[j]->GetResName(), "MSE", 3 ) == 0 ) {
      pRsdTbl[j]->SetResName ( "MET" );
      mse.push_back ( j );
    }

  std::vector<clipper::Atom_list> atm_lst;
  std::vector<double> occ;
  std::vector<std::vector<char> >::const_iterator o_iter;
  mCSolver.SetATMS_MC_OFFSET(0);
  double a = 0.205;
  //mCSolver.SetMILPThreshold ( a );
  const int Lsegment = 4;
  if ( ( hSeqNum - lSeqNum + 1 ) >= Lsegment ) {
    int Nsegments = ceil (double(hSeqNum - lSeqNum + 1)/double(Lsegment+1));
    std::cout << Nsegments << std::endl;
    int hSeq, lSeq;
    hSeq = lSeqNum-1;
    std::vector<std::vector<char> > order[Nsegments];
    for ( int dq = 0; dq<Nsegments; ++dq ) {
        lSeq = hSeq+1;
        hSeq = std::min (lSeq+Lsegment, hSeqNum);
        std::cout << "lSeq: " << lSeq << " hSeq: " << hSeq << std::endl;
        std::vector<std::vector<char> > alts ( hSeq - lSeq + 1 );
        std::vector<char> charseqs ( hSeq - lSeq + 1 );
        GetAltLocs ( lSeq, hSeq, pMMDBChn, alts );
        GetAtomIDs ( pMMDBChn, lSeq, hSeq );
        orderings.clear();
        generate_paths_recursive ( alts, charseqs, 0 );
        atm_lst.clear();
        for ( o_iter=orderings.begin(); o_iter!= orderings.end(); ++o_iter ) {
          atm_lst.push_back ( GetAtoms ( pMMDBChn, lSeq, hSeq, *o_iter ) );
        }
        isSULFUR = (strncmp(pMMDBChn->GetResidue(lSeqNum, "")->GetResName(), "MET", 3) == 0 ||
                    strncmp(pMMDBChn->GetResidue(lSeqNum, "")->GetResName(), "CYS", 3) == 0 ); 
     
        resnum = pMMDBChn->GetResidueNo ( lSeqNum, "" );
     
        mCSolver.SetATMS_SIZE ( atm_lst[0].size ( ) );
        mCSolver.CreateMaskMap ( atm_lst );
        std::cout << "paths to estimate: " << atm_lst.size() << std::endl;
 
        std::vector<int> ch ( atm_lst.size() );
        occ = mCSolver.SetupLPandSolve( atm_lst, ch, QP );

        for ( int j=0; j<occ.size(); ++j )
           if ( occ[j] > 0 )
             order[dq].push_back ( orderings[j] );

        if ( dq == Nsegments - 1 ) {
          orderings.clear();

          /*
           *  Use indices of order, i.e, each nseg 0, 1, ..., Nsegments -1;
           *  Each nseg has indices 0, 1, ..., order[nseg].size() of orders
           *  Need sequence of consective orders, 0, 2, 3, 1, etc.
           */
           Vvi vvi;
          for(int i = 0; i < Nsegments; ++i) {
              Vi vi;
              for(int j = 0; j < order[i].size(); ++j)
                  vi.push_back(j);
              vvi.push_back(vi);
          }

          Vvi output;
           Vi outputTemp;
           cart_product(output, outputTemp, vvi.begin(), vvi.end());

           for ( int i=0; i<output.size(); ++i ) {
               for ( int j=0; j<output[i].size(); ++j )
                   std::cout << output[i][j] << " ";
               for ( int j=0; j<output[i].size(); ++j )
                   for ( int k=0; k<order[j][output[i][j]].size();++k)
                     std::cout << order[j][output[i][j]][k] << " ";
               std::cout << std::endl;
           }

           std::cout << "charseq size: " << hSeqNum - lSeqNum + 1  << std::endl;

           std::vector<char> charseq ( hSeqNum - lSeqNum + 1 );
           for ( int i=0; i<output.size(); ++i )
           {
               int kidx = 0;
               for ( int j=0; j<output[i].size(); ++j )
               {
                   for ( int k=0; k<order[j][output[i][j]].size(); ++k)
                     charseq[kidx++] = order[j][output[i][j]][k];
               }
               orderings.push_back(charseq);
           }
        }
     }
     GetAtomIDs ( pMMDBChn, lSeqNum, hSeqNum );
     for (o_iter=orderings.begin(); o_iter!= orderings.end(); ++o_iter) {
    	for ( int k=0; k<o_iter->size(); ++k)
            std::cout << (*o_iter)[k];
    	std::cout << std::endl;
     }
   } else {
     std::vector<std::vector<char> > alts ( hSeqNum - lSeqNum + 1 );
     std::vector<char> charseqs ( hSeqNum - lSeqNum + 1 );
     GetAltLocs ( lSeqNum, hSeqNum, pMMDBChn, alts );
     GetAtomIDs ( pMMDBChn, lSeqNum, hSeqNum );
     generate_paths_recursive ( alts, charseqs, 0 );
   }

   atm_lst.clear();

    for (o_iter=orderings.begin(); o_iter!= orderings.end(); ++o_iter) {
      atm_lst.push_back ( GetAtoms ( pMMDBChn, lSeqNum, hSeqNum, *o_iter));
    }

    if ( 1.1 < mCSolver.RSLN && mCSolver.RSLN <= 1.4 )
      a = 0.251;
    if ( 1.4 < mCSolver.RSLN && mCSolver.RSLN <= 1.6 )
      a = 0.33;
    else if ( 1.6 < mCSolver.RSLN/* && mCSolver.RSLN <= 1.8 */)
      a = 0.34;
    mCSolver.SetMILPThreshold (a);

    std::vector<int> cho ( atm_lst.size() );
    mCSolver.SetATMS_SIZE ( atm_lst[0].size ( ) );
    mCSolver.CreateMaskMap ( atm_lst );
    occ.clear();
    occ = mCSolver.SetupLPandSolve( atm_lst, cho, MIQP );

    pstr resname = pMMDBChn->GetResidue ( lSeqNum, "" )->GetResName ( );
    std::cout << lSeqNum << " " << resname << std::endl;

    for ( int seqnum=lSeqNum; seqnum<=hSeqNum; seqnum++ )
      pMMDBChn->GetResidue ( seqnum, "" )->DeleteAllAtoms ( );

    const double LOWOCC = 0; //0.1429*mCSolver.RSLN-0.0857;
  
    int NaltLoc=0;
    for ( int j=0; j<occ.size(); ++j )
      if ( occ[j] > LOWOCC )
        NaltLoc++;
    const char altLoc_lst = char('A' + NaltLoc - 1 );
    
    int k=0;
    for ( int j=0; j<occ.size(); ++j )
      if ( occ[j] > LOWOCC )
        SetAtoms ( pMMDBChn, lSeqNum, hSeqNum, atm_lst[j], char(altLoc_lst - k++), occ[j] );
    
    for ( int j=0; j<mse.size(); ++j)
      pRsdTbl[mse[j]]->SetResName ( "MSE " );
    
    mCMMDBManager.FinishStructEdit();
    mCMMDBManager.PDBCleanup ( PDBCLEAN_ATNAME | PDBCLEAN_INDEX | PDBCLEAN_ALTCODE_STRONG );
  
   std::ofstream occOutFile ( "occ.txt" );
   occOutFile << lSeqNum << " " << resname << " ";
   for ( int j=0; j<occ.size(); ++j )
     if ( occ[j] > LOWOCC )
       occOutFile << occ[j] << " ";
   occOutFile.flush();

   mCMMDBManager.SetFlag (MMDBF_IgnoreDuplSeqNum );

   CMMDBManager mRefmngr;
   PCModel pRefmdl = newCModel ( );
   mRefmngr.AddModel ( pRefmdl );
   PCChain pRefChn = newCChain ( );
   pRefmdl->AddChain ( pRefChn );
   for ( int seqnum=lSeqNum; seqnum<=hSeqNum; seqnum++ )
     pRefChn->AddResidue ( pMMDBChn->GetResidue ( seqnum, "" ) );
   mRefmngr.FinishStructEdit();
   WritePDB(pRefChn, "frag.pdb");

   mCMMDBManager.FinishStructEdit();
   mCMMDBManager.WritePDBASCII( "fragFinal.pdb" );
}

// recursive algorithm to to produce cart. prod.
// At any given moment, "me" points to some Vi in the middle of the
// input data set.
//   for int i in *me:
//      add i to current result
//      recurse on next "me"
//
void CFragmentSelector::cart_product(
    Vvi& rvvi,  // final result
    Vi&  rvi,   // current result
    Vvi::const_iterator me, // current input
    Vvi::const_iterator end) // final input
{
    if(me == end) {
        // terminal condition of the recursion. We no longer have
        // any input vectors to manipulate. Add the current result (rvi)
        // to the total set of results (rvvvi).
        rvvi.push_back(rvi);
        return;
    }

    // need an easy name for my vector-of-ints
    const Vi& mevi = *me;
    for (Vi::const_iterator it = mevi.begin(); it != mevi.end(); it++) {
        // final rvi will look like "a, b, c, ME, d, e, f"
        // At the moment, rvi already has "a, b, c"
        rvi.push_back(*it);  // add ME
        cart_product(rvvi, rvi, me+1, end); //add "d, e, f"
        rvi.pop_back(); // clean ME off for next round
    }
}

void CFragmentSelector::GetAtomIDs ( const PCChain pChn, const int lowseq, const int highseq )
{
  const std::string mcAtms[4] = {"N", "CA", "C", "O"};
  atm_ids.clear ( );
  
  for ( int seqnum=lowseq; seqnum<=highseq; ++seqnum )
  {
  	std::string name = pChn->GetResidue ( seqnum, "" )->GetResName ( );
    
    int j=-1;
    int r;
    while ( j < 19 && ( r = strncmp ( name.c_str(), SCdhdrls[++j][0], 3 ) ) != 0 );
   
    if ( r != 0 )
    {
      std::cout << "ERROR: unknown residue." << std::endl;
      exit (1);
    } 
  
    int na = atoi ( SCdhdrls[j][1] );

    for ( int i=0; i<4; i++ )
      atm_ids.push_back ( mcAtms[i] );

    for ( int i=0; i<na; i++ )
    {
      std::string id = SCdhdrls[j][2+i];
      if ( id[0] != 'H' )
        atm_ids.push_back ( SCdhdrls[j][2+i] );
    }
  }
    
//  for ( int i=0; i<atm_ids.size(); ++i )
//    std::cout << atm_ids[i] << std::endl;
//  exit(1);
}

Atom_list CFragmentSelector::GetAtoms ( const PCChain pChn, const int lowseq, const int highseq, const std::vector<char> path )
{
  Atom_list atm_lst;
  AltLoc aLoc;
  int j=0;

  for ( int seqnum=lowseq; seqnum<=highseq; ++seqnum )
  {
  	PCResidue pRsd = pChn->GetResidue ( seqnum, "" );
  	
  	PPCAtom pAtmTbl;
  	int nAtms;
  	char S[128];
  	pRsd->GetAtomTable ( pAtmTbl, nAtms );
//  	for ( int ii=0; ii<nAtms; ++ii )
//  	  std::cout << pAtmTbl[ii]->GetAtomID ( S ) << std::endl; 

    do {
        {
            std::string id = ( atm_ids[j] != "SE" ) ? atm_ids[j] : "SD"; /*tmp fix for MSE*/
            aLoc[0] = path[seqnum-lowseq];
            aLoc[1] = char(0);
            const PCAtom pRsdAtm = pRsd->GetAtom ( id.c_str(), NULL, aLoc );
            Atom atm;
            std::string el = ( atm_ids[j] != "SE" ) ? pRsdAtm->element : "Se"; /*tmp fix for MSE*/
            atm.set_element(el);
            double ranT =  0.09*pRsdAtm->tempFactor*((double)rand()/RAND_MAX-0.5);
            atm.set_u_iso( clipper::Util::b2u ( pRsdAtm->tempFactor + ranT) );
            if ( atm_ids[j] == "N" || atm_ids[j] == "CA" || atm_ids[j] == "C" )
                atm.set_occupancy ( 0.0 );
            else
                atm.set_occupancy ( 1.0 );
            atm.set_coord_orth(Coord_orth(pRsdAtm->x,pRsdAtm->y,pRsdAtm->z));
            atm_lst.push_back ( atm );
        }
        j++;
    } while ( j < atm_ids.size() && atm_ids[j] != "N" );
  }
  return atm_lst;
}

void CFragmentSelector::SetAtoms ( const PCChain pChn, const int lowseq, const int highseq, const Atom_list& atm_lst, const char altLoc, const double occupancy )
{
  PPCAtom pAtmTbl;
  int nAtms;
  int j = 0;
  
//  std::cout << "From SetAtoms:" << std::endl;
  char S[128];
  
  for ( int seqnum=lowseq; seqnum<=highseq; ++seqnum )
  {
  	PCResidue pRsd = pChn->GetResidue ( seqnum, "" );
  
    do
    {
      //if ( isMC ( pAtmTbl[j] ) ) continue;
      double q = occupancy;
      PCAtom pAtm = new CAtom ( );
      pRsd->InsertAtom ( pAtm, atm_ids[j].c_str() );
       
      if ( isMSE[seqnum-lowseq] && atm_lst[j].element() == "S" )
        q *= 0.4706;
      
      const clipper::Coord_orth atm_crds ( atm_lst[j].coord_orth ( ) );
      pAtm->SetCoordinates ( atm_crds.x(), atm_crds.y(), atm_crds.z(), q, clipper::Util::u2b ( atm_lst[j].u_iso() ) );
      
     if ( isMSE[seqnum-lowseq] && atm_lst[j].element() == "S" )
       pAtm->SetElementName ( "SE" );
     else
       pAtm->SetElementName ( atm_lst[j].element().c_str() );

      //pAtm->SetCoordinates ( 0.0, 0.0, 0.0, 1.0, 10.0 );
      pAtm->SetAtomName ( atm_ids[j].c_str() );
      pAtm->altLoc[0] = altLoc;
      pAtm->altLoc[1] = char(0); 
      //res->InsertAtom ( pAtm, atm_ids[j].c_str() );
      //PCMMDBFile(res->GetCoordHierarchy())->FinishStructEdit ( );
      j++;
    } while ( j < atm_ids.size() && atm_ids[j] != "N" );
  }
//  res->GetAtomTable ( pAtmTbl, nAtms );
//  std::cout << std::endl;
//  for ( int i=0; i<nAtms; i++ )
//  {
//    std::cout << pAtmTbl[i]->GetAtomID( S) << " " << pAtmTbl[i]->serNum  << std::endl;
//  }
//  std::cout << std::endl;
  PCMMDBFile(pChn->GetCoordHierarchy())->FinishStructEdit ( );
//  WritePDB ( res->GetChain ( ), "Final_.pdb" );
 
//  std::cout << std::endl << std::endl;
}
