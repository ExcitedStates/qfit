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


#ifndef CSIDECHAINSAMPLER_H_
#include "CSideChainSampler.h"
#endif


#define B2U 78.9568

struct MyRotater: AtomFunctor {
public:
  MyRotater(Vector3 orig, Vector3 axis, float degrees) {
    rotMat = PMath::FindRotationMatrix(axis, DtoR(degrees));
    origin = orig;
  }

  void operator()(PAtom *atom, PBond *bondFrom) {
    Vector3 currPos = atom->getPos()-origin;
    //cout<<atom->getPos()<<endl;
    atom->changePosition(rotMat*currPos+origin);
//    cout<<atom->getPos()<<endl<<endl;
  }
  
private:
  Matrix3 rotMat;
  Vector3 origin;
};


struct MyTSetter : AtomFunctor 
{
public:
  MyTSetter(Real tempfactor) : TempFactor(tempfactor) {}
  
  void operator()(
          PAtom *atom, 
          PBond *bondFrom
          )
  {
    atom->setTempFactor(TempFactor);
  }
  
private:
  Real TempFactor;
};
  

template<class T, class U> 
struct ComparePair1stDsc 
{
  bool operator()(
          const std::pair<T,U>& a, 
          const std::pair<T,U>& b
          ) const 
  {
    return a.first > b.first; 
  }
};


void leftouter_product(
        const std::vector<double>& U, 
        const std::vector<double>& V, 
        std::vector<double>& UxV
        )
{
  UxV[0] = V[1]*U[2] - U[1]*V[2];  
  UxV[1] = V[2]*U[0] - U[2]*V[0];
  UxV[2] = V[0]*U[1] - U[0]*V[1];
}


void normalize3D ( std::vector<double>& U, double& _l )
{
  _l = sqrt ( U[0]*U[0]+U[1]*U[1]+U[2]*U[2] );
  if ( _l > 0 )
    for ( int i=0; i<U.size(); ++i )
    U[i] /= _l;
}


bool WritePDB ( const PCChain pChain, const std::string fn )
{
  CFile f;
  bool success = false;
  f.assign ( const_cast<char*>(fn.c_str()), true, false, GZM_CHECK );
  if ( ( success = f.rewrite() ) )
  {
    pChain->PDBASCIIAtomDump ( f );
    f.WriteLine ( pstr("END") );
    f.shut();
  }
  else
    cout << "Error: Can not open file " << fn << " for writing." << endl;
  return success;
}


CSideChainSampler::CSideChainSampler()
{
}


void CSideChainSampler::Init(
        std::string inpdb, 
        std::string inmtz, 
        const std::string ChnID, 
        const int seqnum, 
        const bool fp
        )
{
  clipper::Atom_list allatm_lst;
  clipper::ftype32 mean, sigsig;
  
  clipper::Xmap<float> xmap;
  clipper::HKL_info hkl_info;

  int RC = myCMMDBManager.ReadPDBASCII(
          const_cast<char*>(inpdb.c_str())
          );
  
  flip_peptide = fp;

  srand ( 100 );
   
  // get a list of all the atoms except N, CA and C
  clipper::mmdb::PPCAtom psel;
  int hndl, nsel;
  hndl = myCMMDBManager.NewSelection();
  myCMMDBManager.SelectAtoms( hndl, 0, 0, SKEY_NEW );
  myCMMDBManager.SelectAtoms(hndl, 0, (pstr)ChnID.c_str(), 
          seqnum, "*", seqnum,"*", "*", 
          "!N,CA,C", "*", "*", SKEY_XOR );
  myCMMDBManager.GetSelIndex(hndl, psel, nsel);
  allatm_lst = clipper::MMDBAtom_list ( psel, nsel );
  myCMMDBManager.DeleteSelection( hndl );
  
  std::cout << "Atom list size: " << allatm_lst.size() << std::endl;

  // Scale the observed map
  PScaleXMap mPScaleXMap;
  mPScaleXMap.DoScaleXMap(xmap, hkl_info, inmtz, allatm_lst, mean, sigsig);
  
  // Get some tensors to determine direction for ADPATOM later on
  Metric_tensor_direct = hkl_info.cell().metric_real();
  Metric_tensor_direct2 = 
      hkl_info.cell().matrix_orth().transpose() *
      hkl_info.cell().matrix_orth();
  Ueigen = clipper::Matrix<>(3, 3, 0.0);

  // Initialize the Solver
  mCSolver.Init(allatm_lst, xmap, hkl_info, mean, sqrt(sigsig));
}


CSideChainSampler::~CSideChainSampler()
{
}


void CSideChainSampler::CalcBndvctrRotPlaneIntsct(
        std::vector<double>& intsct, 
        std::vector<double>& bndvctr,
        const std::vector<double>& bndtail, 
        const std::vector<double>& trmnlatm
        ) const
{
  std::vector<double> convctr(3); 
  double t, init0;
  init0 = 0;

  for ( int i=0; i<trmnlatm.size(); ++i )
      convctr[i] = trmnlatm[i] - bndtail[i];
  
  t = std::inner_product (&convctr[0], &convctr[3], &bndvctr[0], init0 ) ;
  t /= std::inner_product (&bndvctr[0], &bndvctr[3], &bndvctr[0], init0 ) ;

  for ( int i=0; i<trmnlatm.size(); ++i )
    intsct[i] = bndtail[i] + t*bndvctr[i];  
}


double CSideChainSampler::CalcMinimizingAngle(
        std::vector<double>& bndvctr, 
        const std::vector<double>& bndtail, 
        const std::vector<double>& trmnlatms, 
        const std::vector<double>& anchrAtms
        ) const
{
  std::vector<double> orgn(3), anchrvctr(3), plnvctr(3), rotatm(3);
   std::vector<double> rhat (3), shat(3);
  
  double b, c, r, init0, l0;
  init0 = 0; b = 0; c = 0; r = 0; l0 = 0;
  for ( std::size_t i=0; i<1; i++ ) {
    for ( int j=3*i; j<3*(i+1); ++j )
      plnvctr[j-3*i] = trmnlatms[j];
    CalcBndvctrRotPlaneIntsct ( orgn, bndvctr, bndtail, plnvctr );
    for ( int j=0; j<3; ++j )
      plnvctr[j] -= orgn[j];
    leftouter_product ( plnvctr, bndvctr, shat );
    rhat = plnvctr;
    normalize3D ( shat, r );
    normalize3D ( rhat, r );
    for ( int j=3*i; j<3*(i+1); ++j )
      anchrvctr[j-3*i] = anchrAtms[j];
    for ( int j=0; j<3; ++j )
      anchrvctr[j] -= orgn[j];
    b += r* std::inner_product ( &rhat[0], &rhat[3], &anchrvctr[0], init0 );
    c += r* std::inner_product ( &shat[0], &shat[3], &anchrvctr[0], init0 );
  }
  return atan2 ( c, b );
}


void CSideChainSampler::SetADPAxis_Direct(const PCAtom pAtm)
{
  clipper::Mat33<> Ustar, U;
  Ustar( 0,0 ) = pAtm->u11;
  Ustar( 1,1 ) = pAtm->u22;
  Ustar( 2,2 ) = pAtm->u33;
  Ustar( 0,1 ) = Ustar ( 1,0 ) = pAtm->u12;
  Ustar( 0,2 ) = Ustar ( 2,0 ) = pAtm->u13;
  Ustar( 1,2 ) = Ustar ( 2,1 ) = pAtm->u23;
  
  std::cout << "Ustar: " << Ustar.format ( ) << std::endl;
  U = Metric_tensor_direct2 * Ustar * Metric_tensor_direct2;
  std::cout << "U: " << U.format ( ) << std::endl;
  for ( int i=0; i<3; ++i )
    for ( int j=0; j<3; ++j )
      Ueigen (i,j) = Ustar(i, j); 
  
  std::vector<clipper::ftype> evals = Ueigen.eigen();
  
  std::cout << "Ueigen: ";
  for ( int i=0; i<3; ++i ) {
    for ( int j=0; j<3; ++j )
      std::cout << Ueigen ( i, j ) << " ";
    std::cout << std::endl;
  }

  for ( int i=0; i<evals.size(); ++i )
    std::cout << evals[i] << " ";
  std::cout << std::endl;
    
  clipper::Vec3<> ESum;
  for (int i=0; i<3; i++ ) {
    ESum[i] = 0;
    for (int j=0; j<3; j++ )
      ESum[i] += evals[j]*Ueigen ( i, j );
  }
  EVectorSum = ESum.unit();
}


void CSideChainSampler::DoFindConfs(
        const std::string ChnID, 
        const int lSeqNum, 
        const int hSeqNum
        )
{
    PPCAtom pAtmTbl;
    int nAtms;
    char S[128];
    bool isSULFUR;
    
    PCChain pMMDBChn = myCMMDBManager.GetChain(1, ChnID.c_str());
  
    if (!pMMDBChn->isAminoacidChain()) {
      /*contains at least one AA*/
      std::cout << "Not an AA chain." << std::endl; 
      exit(1);
    }

    // Check for MSE residues, and change to MET. Keep track of residue index.
    std::vector<int> mse;
    PPCResidue pRsdTbl;
    int nRsds;
    pMMDBChn->GetResidueTable(pRsdTbl, nRsds);
    for ( int j=0; j<nRsds; j++ )
        if (strncmp(pRsdTbl[j]->GetResName(), "MSE", 3) == 0) {
            pRsdTbl[j]->SetResName ( "MET" );
            mse.push_back(j);
        }
     
    const bool isGLY = strncmp(
            pMMDBChn->GetResidue(lSeqNum, "")->GetResName(), "GLY", 3 
            ) == 0;

    std::string ADP_ATOM;
    if ( isGLY ) {
       ADP_ATOM = "O";
       mCSolver.SetATMS_MC_OFFSET(3);
    } else {
       ADP_ATOM = "CB";
       mCSolver.SetATMS_MC_OFFSET(3);
    }
    
    SetADPAxis_Direct(
            pMMDBChn->GetAtom(lSeqNum, "", ADP_ATOM.c_str(), NULL, NULL )
            );
    CATempFactor = 
        pMMDBChn->GetAtom(lSeqNum, "", ADP_ATOM.c_str(), NULL, NULL)->tempFactor;
    
    // Write chain to file and read in with looptk-based reader.
    std::ostringstream oss_seq;
    oss_seq << lSeqNum;
    std::string fn = ".chn_" + oss_seq.str() + ".pdb";
    WritePDB(pMMDBChn, fn.c_str());
    PProtein *pChn = PDBIO::readFromFile(fn);

    std::vector<clipper::Atom_list> atm_lst;
    pMMDBChn->GetResidue(lSeqNum, "")->GetAtomTable(pAtmTbl, nAtms);
    
    isMSE = strncmp(
            pMMDBChn->GetResidue(lSeqNum, "")->GetResName(), "MET", 3
            ) == 0;
    
    // Unused boolean
    isSULFUR = (strncmp(
                    pMMDBChn->GetResidue(lSeqNum, "")->GetResName(), "MET", 3
                    ) == 0 ||
                strncmp(
                    pMMDBChn->GetResidue(lSeqNum, "")->GetResName(), "MSE", 3
                    ) == 0 ||
                strncmp(
                    pMMDBChn->GetResidue(lSeqNum, "")->GetResName(), "CYS", 3
                    ) == 0
                ); 

    if ( flip_peptide )
      std::cout << "flip_peptide = true" << std::endl;
    else
      std::cout << "flip_peptide = false" << std::endl;
     
    resnum = pMMDBChn->GetResidueNo(lSeqNum, "" );
    std::cout << "resnum " << resnum << std::endl;
      
    GetAtomIDs(pChn->getResidue(resnum));
    
    std::vector<double> occ;
    vector<PAtom*> vResAtms = *(pChn->getResidue ( resnum )->getAtoms());
    
    double a = 0.2;
    if ( 1.4 < mCSolver.RSLN && mCSolver.RSLN <= 1.6 )
      a = 0.25;
    else if ( 1.6 < mCSolver.RSLN && mCSolver.RSLN <= 1.9 )
      a = 0.33;
    else if ( 1.9 < mCSolver.RSLN )
      a = 0.34;

    MC_AMPL = 0.20;
    MC_SIG = 0.125;
    gSAMPLE_SIZE = 13;
  
    bool occfromfile = false;
    std::ifstream occfile("minocc.txt");
    double b = 0;
    if (occfile.is_open()) {
       occfile >> b >> MC_AMPL >> MC_SIG >> gSAMPLE_SIZE;
    }

    std::cout << "occupancy MC_AMPL MCSIG: " 
        << a << " " << MC_AMPL << " " << MC_SIG << std::endl;

    mCSolver.SetMILPThreshold(a);
    const bool isLong_Ring = (
            isLongSC(pChn->getResidue(resnum)) || 
            isRing(pChn->getResidue(resnum))
            );

    if ( isLong_Ring ) {
      /*long sc*/
      PProtein *pChnSmall = new PProtein(pChn, resnum - 3, resnum + 3);
      unsigned chiMax = 
          PResources::numChiIndices(pChn->getResidue(resnum)->getName());
      std::cout << "chiMax: " << chiMax << std::endl;
     
      SampleCB_MC_AtChi(pChnSmall->getResidue(3), atm_lst);

      for (int k=1; k<chiMax; k++ ) {
        /*Select correct number of atoms to include for this chi angle*/
        occ.clear();

        vector<string> rotList = 
            PResources::GetChiIndex(pChnSmall->getResidue(3)->getName(), k+1);
        std::string atm_name = rotList[2] == "SD" ? "SE" : rotList[2];

        if (pChnSmall->getResidue(3)->getName() == PID::PHE ||
            pChnSmall->getResidue(3)->getName() == PID::TYR) {
            atm_name = "CZ";
        }

        std::vector<std::string>::iterator 
            atm_ids_iter = std::find(atm_ids.begin(), atm_ids.end(), atm_name);
        
        if (atm_ids_iter == atm_ids.end()) {
          std::cout << "ERROR: Atom not found" << std::endl;
          exit(1);
        }
         
        size_t NChiAtms = atm_ids_iter - atm_ids.begin();
        
        std::cout << "Looked for: " << atm_name << std::endl;
        std::cout << "NChiAtms: " << NChiAtms << std::endl;
        std::cout << "Conformers to estimate: " << atm_lst.size() << std::endl;
        
        mCSolver.SetATMS_SIZE(NChiAtms + 1);

        mCSolver.CreateMaskMap(atm_lst);
     
        std::vector<int> ch(atm_lst.size());
        occ = mCSolver.SetupLPandSolve(atm_lst, ch);
      
        std::vector<double> occtmp;
        std::vector<clipper::Atom_list> atm_lst_;
        for ( int j=0; j<occ.size(); ++j )
            if ( occ[j] > 0.05 ) {
                atm_lst_.push_back ( atm_lst[j] );
                occtmp.push_back(occ[j]);
            }
        occ.clear();
        atm_lst.clear();
        for (int j=0; j<atm_lst_.size(); ++j) {
          atm_lst.push_back(atm_lst_[j]);
          occ.push_back(occtmp[j]);
        }
        Sample_NextChi(pChnSmall->getResidue(3), k+1, occ, atm_lst);
      }
      pChnSmall->Obliterate();

    // For short chains
    } else {

      PProtein *pChnSmall = new PProtein ( pChn, resnum - 3, resnum + 3 );
      SampleCB_MC(pChnSmall->getResidue(3), atm_lst);
      pChnSmall->Obliterate();
    }

    if (mCSolver.RSLN <= 1.1)
      a = 0.201;
    if ( 1.1 < mCSolver.RSLN && mCSolver.RSLN <= 1.4 )
      a = 0.251;
    if ( 1.4 < mCSolver.RSLN && mCSolver.RSLN <= 1.6 )
      a = 0.33;
    else if ( 1.6 < mCSolver.RSLN)
      a = 0.34;
    if ( !occfromfile /*&& !isGLY */ )
      mCSolver.SetMILPThreshold(a);

    occ.clear();
    mCSolver.SetATMS_SIZE(atm_lst[0].size());
    mCSolver.CreateMaskMap(atm_lst);
    std::cout << "SetATMS_SIZE: " << atm_lst[0].size() << std::endl;
    std::cout << "Conformers to estimate: " << atm_lst.size() << std::endl;
    std::vector<int> ch(atm_lst.size());
    occ = mCSolver.SetupLPandSolve(atm_lst, ch, QP);
      
    mCSolver.GoodnessOfFit(atm_lst, occ);
      
    std::vector<clipper::Atom_list> atm_lst_;
 
    for ( int j=0; j<occ.size(); ++j )
        if ( occ[j] > 0 ) 
          atm_lst_.push_back(atm_lst[j]);

    mCSolver.CreateMaskMap ( atm_lst_ );
    std::cout << "Conformers to estimate: " << atm_lst_.size() << std::endl;
    std::vector<int> cho ( atm_lst_.size() );
    
    occ.clear();
    occ = mCSolver.SetupLPandSolve(atm_lst_, cho, MIQP);

    pstr resname = pMMDBChn->GetResidue(lSeqNum, "")->GetResName();
    std::cout << lSeqNum << " " << resname << std::endl;
    PCResidue pRsd = pMMDBChn->GetResidue(lSeqNum, "");
    pRsd->DeleteAllAtoms();
    
    const double LOWOCC = 0; //0.1429*mCSolver.RSLN-0.0857;
  
    int NaltLoc=0;
    for ( int j=0; j<occ.size(); ++j )
      if ( occ[j] > LOWOCC )
        NaltLoc++;
    
    if ( NaltLoc == 0 ) {
      std::cout << "Residue " << lSeqNum << " " << 
          resname << " should probably be truncated." << std::endl;
      NaltLoc = 1;
      occ[0] = a;
    }
    
    const char altLoc_lst = char('A' + NaltLoc - 1 );
    
    int k=0;
    for ( int j=0; j<occ.size(); ++j )
      if ( occ[j] > LOWOCC )
        SetAtoms ( pRsd, atm_lst_[j], char(altLoc_lst - k++), occ[j] );
    
    for ( int j=0; j<mse.size(); ++j)
      pRsdTbl[mse[j]]->SetResName ( "MSE " );
    
    myCMMDBManager.FinishStructEdit();
    myCMMDBManager.PDBCleanup(
            PDBCLEAN_ATNAME | PDBCLEAN_INDEX | PDBCLEAN_ALTCODE_STRONG
            );
  
   std::ofstream occOutFile("occ.txt");
   occOutFile << lSeqNum << " " << resname << " ";
   for ( int j=0; j<occ.size(); ++j )
     if ( occ[j] > LOWOCC )
       occOutFile << occ[j] << " ";
   occOutFile.flush();

   myCMMDBManager.SetFlag(MMDBF_IgnoreDuplSeqNum);

   CMMDBManager mRefmngr;
   PCModel pRefmdl = newCModel();
   mRefmngr.AddModel(pRefmdl);
   PCChain pRefChn = newCChain();
   pRefmdl->AddChain(pRefChn);
   pRefChn->AddResidue(pRsd);
   mRefmngr.FinishStructEdit();
   WritePDB(pRefChn, ("Res" + oss_seq.str() + ".pdb").c_str());

   pRsd->DeleteAtom ( NULL, "H", NULL );
   myCMMDBManager.FinishStructEdit();
   myCMMDBManager.WritePDBASCII(
           const_cast<char*>(("Final" + oss_seq.str() + ".pdb").c_str())
           );
}


// Unused function
void CSideChainSampler::CalcConfs(
        PProtein* pChn, 
        const int start_res_idx, 
        const int stop_res_idx
        )
{
  std::vector<clipper::Atom_list> mAtom_list;

  for ( int nres=start_res_idx; nres<=stop_res_idx; nres++ ) {
    PProteinResidue* pPRes_cur = pChn->getResidue ( nres );
    GenTrialPositions ( pPRes_cur, mAtom_list );
  }
}


void CSideChainSampler::GetAtomIDs(PResidue *res)
{
  const std::string mcAtms[4] = {"N", "CA", "C", "O"};
  vector<PAtom*> vResAtms = *(res->getAtoms());
  std::string name = res->getName();
  atm_ids.clear();
  
  int j=-1;
  int r;
  while (j < 19 && (r = strncmp(name.c_str(), SCdhdrls[++j][0], 3)) != 0 );
      if ( r != 0 ) {
        std::cout << "ERROR: unknown residue." << std::endl;
        exit (1);
      }
  
  int na = atoi(SCdhdrls[j][1]);
  
  /*First MC atoms ... */
  for (int i=0; i<4; i++ )
    atm_ids.push_back(mcAtms[i]);

  /* then SC */
  for (int i=0; i<na; i++)
    atm_ids.push_back(SCdhdrls[j][2+i]);
}


/* To fit atoms in prvs or next residue, change:
 * 1. GetAtoms below
 * 2. The number of atoms in DoFindConf (NChiAtoms++)
 * 3. SetAtoms (finla)
 * 4. Setatoms (intermediate)
 */

Atom_list CSideChainSampler::GetAtoms(
        PResidue *res, 
        const bool SConly
        )
{
  Atom_list atm_lst;
  vector<PAtom*> vResAtms = *(res->getAtoms());
  atm_lst.reserve( vResAtms.size());

  for ( int j=0; j<atm_ids.size(); ++j) {
      /* tmp fix for MSE */
      std::string id = (atm_ids[j] != "SE" ) ? atm_ids[j] : "SD";
      PAtom* pResAtm = res->getAtom(id);
      Atom atm;
      /* tmp fix for MSE */
      std::string el = ( atm_ids[j] != "SE" ) ? pResAtm->getName() : "Se";
      atm.set_element(el);     
      atm.set_u_iso(clipper::Util::b2u(pResAtm->getTempFactor()));
      atm.set_occupancy (1.0);
      const Vector3 post = pResAtm->getPos();
      atm.set_coord_orth(Coord_orth(post.x,post.y,post.z));
      atm_lst.push_back(atm);
  }
  return atm_lst;
}


void CSideChainSampler::SetAtoms(
        PCResidue res, 
        const Atom_list& atm_lst, 
        const char altLoc, 
        const double occupancy
        )
{
    PPCAtom pAtmTbl;
    int nAtms;

    int k = 0;
    //  std::cout << "From SetAtoms:" << std::endl;
    char S[128];
    for ( int j=0; j<atm_lst.size(); ++j) {
         /* if ( j == 4 ) {continue;} VDB: change to get next N*/
        double q = occupancy;
        //if ( isMC ( pAtmTbl[j] ) ) continue;
        PCAtom pAtm = new CAtom ( );
        std::cout << atm_ids[k].c_str() << std::endl;
        res->InsertAtom ( pAtm, atm_ids[k].c_str() );
            
        const clipper::Coord_orth atm_crds ( atm_lst[j].coord_orth ( ) );
            
        if ( isMSE && atm_lst[j].element() == "S" )
              q *= 0.4706;
              
        pAtm->SetCoordinates(atm_crds.x(), atm_crds.y(), atm_crds.z(), 
                             q, clipper::Util::u2b(atm_lst[j].u_iso()));
        //pAtm->SetCoordinates ( 0.0, 0.0, 0.0, 1.0, 10.0 );
        if (isMSE && atm_lst[j].element() == "S")
            pAtm->SetElementName ( "Se" );
        else
            pAtm->SetElementName ( atm_lst[j].element().c_str() );
        //std::cout << atm_ids[j].c_str() << endl;
        pAtm->SetAtomName(atm_ids[k].c_str());
        pAtm->altLoc[0] = altLoc;
        pAtm->altLoc[1] = char(0); 
        //res->InsertAtom ( pAtm, atm_ids[j].c_str() );
        //PCMMDBFile(res->GetCoordHierarchy())->FinishStructEdit ( );
        k++;
    }
    res->GetAtomTable ( pAtmTbl, nAtms );
    PCMMDBFile(res->GetCoordHierarchy())->FinishStructEdit ( );
}


void CSideChainSampler::SetAtoms(
        PResidue *res, 
        const Atom_list& atm_lst, 
        const int c_idx
        )
{
    const int KNMTCS_SIZE = 3;
    //if ( atm_lst.size () != atm_ids.size ( ) + 1 ) /*VDB: change to get next N*/
    if ( atm_lst.size () != atm_ids.size ( ) )
        std::cout << "ERROR: Atom lists differ in size!" << std::endl; 
    PAtom *patm;

    int k=0;
    for ( int i=0; i<atm_lst.size (); ++i ) {
        std::string id = ( atm_ids[k] != "SE" ) ? atm_ids[k] : "SD";
        patm = res->getAtom ( id );
        clipper::Coord_orth crd_orth = atm_lst[i].coord_orth();
        Vector3 newPos ( crd_orth.x(),crd_orth.y(), crd_orth.z() );
        patm->changePosition ( newPos );
        patm->setTempFactor ( clipper::Util::u2b ( atm_lst[i].u_iso() ) );
        patm->setOccupancy ( 1.0 );
        k++;
    }
    int r = rand();
    std::ostringstream oss;
    oss << r;
    std::ostringstream oss2;
    oss2 << c_idx;
    PDBIO::writeToFile(res->getChain(), 
            "afterSetAtoms_" + oss.str()  +"c-idx" + oss2.str() + ".pdb");
} 
          
Atom_list CSideChainSampler::GetAtoms(
        PProtein *pProt, const int excl_res
        )
{
    Atom_list atm_lst;
    for ( int i=0; i<pProt->size(); ++i ) {
        if ( i == excl_res ) continue;
        PProteinResidue *res = pProt->getResidue ( i );
        vector<PAtom*> vResAtms = *(res->getAtoms());

        for ( int j=0; j<vResAtms.size(); ++j) {
            Atom atm;
            atm.set_element(vResAtms[j]->getName());
            atm.set_u_iso(clipper::Util::b2u ( vResAtms[j]->getTempFactor() ) );
            atm.set_occupancy(1.0);
            const Vector3 post = vResAtms[j]->getPos();
            atm.set_coord_orth(Coord_orth(post.x,post.y,post.z));
            atm_lst.push_back ( atm );
        }
    }
    return atm_lst;
}


// Unused function
/* Sample CB position of residue i by adjusting psi_{i-1} and phi_{i} to leave
 * the position CA_{i+1} nearly unchanged.
 * Sample psi_{i-1} and minimize phi_i.
 * This method takes 3-mers?
 */
void CSideChainSampler::SamplePsiPhi ( std::vector<PProtein*> pChnTbl, const int res_idx )
{
  const int insize = pChnTbl.size();
  const int rot_dir = 1;
  Vector3 hd, tl, efftor[2];
  std::vector<double> bnd ( 3 ), bndtl ( 3 ), endefftor ( 6 ), anchrs ( 6);
  
  //PDBIO::writeToFile( pChnTbl.back(), "before.pdb" );
  
  efftor[0] = pChnTbl.back()->getResidue ( 2 )->getAtom ( PID::C_ALPHA )->getPos();
  efftor[1] = pChnTbl.back()->getResidue ( 2 )->getAtom ( PID::C )->getPos();    
  anchrs[0] = efftor[0].x; anchrs[1] = efftor[0].y; anchrs[2] = efftor[0].z;
  anchrs[3] = efftor[1].x; anchrs[4] = efftor[1].y; anchrs[5] = efftor[1].z;
  
//  for ( int i=0; i<insize; ++i)
//  {
//    pChnTbl.push_back ( pChnTbl[i]->Clone() );
    pChnTbl.back()->RotateBackbone ( 1 /*psi*/, BondDirection(rot_dir), -10 );
//  }

  //PDBIO::writeToFile( pChnTbl.back(), "mid.pdb" );
  hd = pChnTbl.back()->getResidue ( 1 )->getAtom ( PID::C_ALPHA )->getPos();
  tl = pChnTbl.back()->getResidue ( 1 )->getAtom ( PID::N )->getPos();
  bnd[0] = hd.x-tl.x; bnd[1] = hd.y-tl.y; bnd[2] = hd.z-tl.z;
  bndtl[0] = tl.x; bndtl[1] = tl.y; bndtl[2] = tl.z;

  efftor[0] = pChnTbl.back()->getResidue ( 2 )->getAtom ( PID::C_ALPHA )->getPos();
  efftor[1] = pChnTbl.back()->getResidue ( 2 )->getAtom ( PID::C )->getPos();    
  endefftor[0] = efftor[0].x; endefftor[1] = efftor[0].y; endefftor[2] = efftor[0].z;
  endefftor[3] = efftor[1].x; endefftor[4] = efftor[1].y; endefftor[5] = efftor[1].z;
  
  std::cout << CalcMinimizingAngle ( bnd, bndtl, endefftor, anchrs ) * 180./M_PI << "\n";
  
  pChnTbl.back()->RotateBackbone ( 2 /*psi*/, BondDirection(rot_dir), 
   -CalcMinimizingAngle ( bnd, bndtl, endefftor, anchrs ) * 180./M_PI);
    
  //PDBIO::writeToFile( pChnTbl.back(), "after.pdb" );
  /* Rotate phi angle to match anchor positions.*/
}


bool CSideChainSampler::BondAngle(
        PProteinResidue* pPres, 
        std::vector<clipper::Atom_list>& Atm_lst
        )
{
    /*Change bondangles simultaneouly for 2 dofs */
    PAtom *patm1, *patm2, *patm3;
    Vector3 bnd1, bnd2, axis;
    //const float dAngle = 10.;
    const float dAngle = 7.5;
    
    patm1 = pPres->getAtom(PID::C_ALPHA );
    patm2 = pPres->getAtom(PID::C_BETA );
    patm3 = pPres->getAtom("CG" );
    
    bnd1 = patm2->getPos() - patm1->getPos();
    bnd2 = patm3->getPos() - patm2->getPos();
    PBond *pBond = pPres->getBond ( PID::C_ALPHA, PID::C_BETA );

    axis  = cross(bnd1, bnd2);
 
    MyRotater rinit1(patm2->getPos(), axis, -dAngle );
    pBond->traverseAtoms(forward, &rinit1, pPres->getChain());
    Atm_lst.push_back(GetAtoms(pPres, true));

    pBond->traverseAtoms(forward, &rinit1, pPres->getChain());
    Atm_lst.push_back(GetAtoms(pPres, true));

    MyRotater rinit2 (patm2->getPos(), axis, 4*dAngle );
    pBond->traverseAtoms(forward, &rinit2, pPres->getChain());
    Atm_lst.push_back ( GetAtoms ( pPres, true ) );

    pBond->traverseAtoms(forward, &rinit1, pPres->getChain());
    Atm_lst.push_back ( GetAtoms ( pPres, true ) );

    pBond->traverseAtoms(forward, &rinit1, pPres->getChain());
    return true;
}


void CSideChainSampler::GenTrialPositions(
        PProteinResidue* pPres, 
        std::vector<clipper::Atom_list>& Atm_lst
        )
{
    if ((pPres->getName()==PID::ALA) || (pPres->getName()==PID::GLY)) {
        Atm_lst.push_back(GetAtoms(pPres, true));
        return;
    }

    vector<vector<Real> > rotamerAngles = PResources::GetRotamer(pPres->getName());
    //std::cout << "NRotamers: " <<  rotamerAngles.size() << std::endl;
    const Real Tf = pPres->getAtom(PID::C_ALPHA)->getTempFactor();
    for (int j=0; j<rotamerAngles.size(); ++j) {
        vector<vector<Real> > rotamerNeighborhood;
        if ( pPres->getName()==PID::PRO )
            rotamerNeighborhood.push_back ( rotamerAngles[j] );
        else
            SampleRotamerNeighborhood(
                    pPres, rotamerAngles[j], rotamerNeighborhood
                    );

        for ( int k=0; k<rotamerNeighborhood.size(); ++k ) {
            RotamerApply(pPres, rotamerNeighborhood[k]);
            for ( int l=0; l<1; l++ ) {
                SetTFactor(pPres, Tf + 0.2 * l * Tf);
                Atm_lst.push_back(GetAtoms(pPres, true ));
                if (isRing(pPres))
                    BondAngle(pPres, Atm_lst);
            }
        }
    }
}

void CSideChainSampler::GenTrialPositions_AtChi(
        PProteinResidue* pPres, 
        const int chi_idx, 
        std::vector<clipper::Atom_list>& Atm_lst
        )
{
    if ((pPres->getName()==PID::ALA) || (pPres->getName()==PID::GLY) || chi_idx == 0 ) {
        Atm_lst.push_back ( GetAtoms ( pPres, true ) );
    return;
    }

    std::cout << "gSAMPLE_SIZE: " << gSAMPLE_SIZE << std::endl;
    const int SAMPLE_SIZE = gSAMPLE_SIZE;
    int bnd = (int) gSAMPLE_SIZE / 2;
    const double delta = (double) 25 / (double) (bnd);
    Real rotamer_delta[SAMPLE_SIZE];

    for ( int i=-bnd; i<bnd+1; i++ ) {
        rotamer_delta[i+bnd] = i*delta;
    }

    std::cout << "SAMPLE_SIZE: " << SAMPLE_SIZE
              << " bnd: " << bnd << " delta: " << delta << std::endl;

    vector<vector<Real> > rotamerAngles = PResources::GetRotamer(pPres->getName());
    std::vector<int> angles_added;
       
    for ( int j=0; j<rotamerAngles.size(); ++j ) {
        if ( chi_idx > 1 && mCSolver.RSLN > 1.1 && 
             pPres->getName()!=PID::TRP ) {
            Real Chi_last = pPres->GetChi(chi_idx - 1);
            // The number here is 1 / 360
            Chi_last -= rint(Chi_last * 0.0027777778) * 360;
            Real dchi = Chi_last - rotamerAngles[j][chi_idx - 2];

            std::cout << "last chi: " << Chi_last 
                      << " rotamer angle: " << rotamerAngles[j][chi_idx-2];

            if ( dchi > 180 ) dchi -= 360;
            if ( dchi <= -180 ) dchi += 360;

            std::cout << " distance: " << dchi << std::endl;

            if ( std::abs (dchi) > (bnd + .5)*delta ) continue;

            std::cout << "Rotamer angle corresponding to current configuration: " 
                      << rotamerAngles[j][chi_idx-2] << std::endl;
            std::cout << "Setting next angle to : " 
                      << rotamerAngles[j][chi_idx-1] << std::endl;
        }
        
        if (std::find(angles_added.begin(), 
                      angles_added.end(), 
                      (int)(rotamerAngles[j][chi_idx-1])) 
                != angles_added.end()) 
            continue; 
        
        angles_added.push_back(rotamerAngles[j][chi_idx - 1]);
        std::ostringstream ossj;
        ossj << j + rand();
        const Real Tf = pPres->getAtom(PID::C_ALPHA)->getTempFactor();
        for ( int k=0; k<SAMPLE_SIZE; k++ ) {
            SetChi(pPres, chi_idx, rotamer_delta[k] + rotamerAngles[j][chi_idx-1]);
            for ( int l=0; l < 1; l++ ) {
                SetTFactor(pPres, Tf + 0.2 * l * Tf, chi_idx);
                Atm_lst.push_back(GetAtoms(pPres, true));
                if (isMSE && chi_idx >= 2) {
                    Atm_lst.push_back(Atm_lst.back());
                    Atm_lst.back()[10].set_element("S");
                }
                if (isRing(pPres) && chi_idx == 1 ) {
                    BondAngle(pPres, Atm_lst);
                }
            }
        }
    }
}


// Unused function
bool CSideChainSampler::Calc_GradientMove(
        PResidue *res, 
        std::vector<clipper::ftype>& grad, 
        const double occ
        )
{
  /*Packae up atoms */
  Atom_list atm_lst;
  PAtom* resAtm = ( res->getAtom(PID::C_BETA) );
          
  Atom atm;
  atm.set_element(resAtm->getName());
  atm.set_u_iso( clipper::Util::b2u ( resAtm->getTempFactor() ) );
  atm.set_occupancy ( 1.0 );
  const Vector3 post = resAtm->getPos();
  atm.set_coord_orth(Coord_orth(post.x,post.y,post.z));
  atm_lst.push_back ( atm );
  
  std::cout << atm.element ( ) << " " << atm.coord_orth().format() << " " << atm.u_iso() << " " << atm.occupancy() << std::endl;
  
  mCSolver.CalcAtomGradient ( atm_lst, grad, occ );
  
  double norm =0.0;
  for ( int i=0; i<grad.size(); ++i )
    norm += (grad[i]*grad[i]);
  norm = sqrt(norm);
  for ( int i=0; i<grad.size(); ++i )
    grad[i] /= norm;

  return true;
}


//Unused function
float CSideChainSampler::Calc_CC ( PResidue *res )
{
    /*Package up atoms */
    Atom_list atm_lst;
    PAtom* resAtm = res->getAtom(PID::O);
    Atom atm;
    atm.set_element(resAtm->getName());
    atm.set_u_iso( clipper::Util::b2u ( resAtm->getTempFactor() ) );
    atm.set_occupancy ( 1.0 );
    const Vector3 post = resAtm->getPos();
    atm.set_coord_orth(Coord_orth(post.x,post.y,post.z));
    atm_lst.push_back ( atm );
    return  mCSolver.CalcCC ( atm_lst );
}


bool CSideChainSampler::SampleCB_MC(
        PProteinResidue* pPres, 
        std::vector<clipper::Atom_list>& Atm_lst
        )
{
   const int KNMTCS_SIZE = 3;
   std::vector<std::pair<Real, std::pair< Vector3, bool > > > mvDirs;

  // Create directions to sample CB-atom
  for (float ampl = 0.10; ampl < 0.31; ampl += 0.10) {
      MC_AMPL = ampl;
      mvDirs.push_back(
              std::make_pair(0.0, 
                  std::make_pair(
                      (MC_AMPL + MC_SIG * ((double)rand() / RAND_MAX)) * 
                      Vector3(Ueigen(0,0), Ueigen(1,0), Ueigen(2,0)), false 
                      )
                  )
              );
      mvDirs.push_back(std::make_pair(0.0,
                  std::make_pair(
                      (-MC_AMPL - MC_SIG *((double)rand()/RAND_MAX)) * 
                      Vector3(Ueigen(0,0), Ueigen(1,0), Ueigen(2,0)), false)
                  )
              );
      mvDirs.push_back (std::make_pair ( 0.0, std::make_pair (
                      (MC_AMPL + MC_SIG * ((double)rand()/RAND_MAX) ) * 
                      Vector3(Ueigen(0,1), Ueigen(1,1), Ueigen(2,1)), false)) 
              );
      mvDirs.push_back ( std::make_pair ( 0.0, std::make_pair ( 
                      ( -MC_AMPL-MC_SIG*((double)rand()/RAND_MAX) ) * 
                      Vector3 ( Ueigen (0,1), Ueigen (1,1), Ueigen (2,1) ), false ) ) 
              );
      mvDirs.push_back(std::make_pair(0.0, std::make_pair(
                      (MC_AMPL+MC_SIG*((double)rand()/RAND_MAX) ) * 
                      Vector3 ( Ueigen (0,2), Ueigen (1,2), Ueigen (2,2) ), false ) ) 
              );
      mvDirs.push_back ( std::make_pair ( 0.0, std::make_pair(
                      ( -MC_AMPL-MC_SIG*((double)rand()/RAND_MAX) ) * 
                      Vector3 ( Ueigen (0,2), Ueigen (1,2), Ueigen (2,2) ), false ) ) 
              );
      mvDirs.push_back ( std::make_pair ( 0.0, std::make_pair ( 
                      ( MC_AMPL+MC_SIG*((double)rand()/RAND_MAX) ) * 
                      Vector3 ( EVectorSum[0], EVectorSum[1], EVectorSum[2] ), false ) ) 
              );
      mvDirs.push_back ( std::make_pair ( 0.0, std::make_pair(
                      ( -MC_AMPL-MC_SIG*((double)rand()/RAND_MAX) ) * 
                      Vector3 ( EVectorSum[0], EVectorSum[1], EVectorSum[2] ), false ) ) 
              );
  }
  
  PProtein *pChn = pPres->getProtein()->Clone();
  /*Check for leak*/

  const bool isGLY = pPres->getName() == PID::GLY;
  const string atomid = ( isGLY ) ? PID::O : PID::C_BETA;
  const int nbb = mvDirs.size();
  const int nflips = 4;
  int k = 0;
  for ( int j=-1; j<nbb; ++j) {
      PProtein *pChnCln = pChn->Clone();

      if (j >= 0) {
         IKSolutions myIKSolns = MoveAtom(
                 pChnCln, 
                 pChnCln->getResidue(KNMTCS_SIZE)->getAtom(atomid), 
                 mvDirs[j].second.first
                 );
      }

      std::ostringstream oss;
      oss << k;
      k++;
      PDBIO::writeToFile(pChnCln, "after_" + oss.str() + ".pdb");

      GenTrialPositions(pChnCln->getResidue(KNMTCS_SIZE), Atm_lst);

      if (flip_peptide) {
          for ( int l=0; l<nflips; ++l ) {
              PProtein *pChnFlipCln = pChnCln->Clone ( );
              if ( l == 0 )
                  FlipTransformPeptide(pChnFlipCln, FD0, 3 );
              else if ( l == 1 )
                  FlipTransformPeptide(pChnFlipCln, FD1, 3 );
              else if ( l == 2 )
                  FlipTransformPeptide(pChnFlipCln, FD2, 3 );
              else if ( l == 3 )
                  FlipTransformPeptide(pChnFlipCln, FD3, 3 );

              PProtein *pChnSmall = new PProtein(
                      pChnFlipCln, KNMTCS_SIZE, KNMTCS_SIZE
                      );

              PDBIO::writeMainchainToFile(pChnSmall, ".t.pdb");
              pChnFlipCln->Obliterate();
              pChnFlipCln = PDBIO::readFromFile(".t.pdb");
              PDBIO::writeToFile(pChnFlipCln, ".tsc.pdb");

              Vector3 hd, tl, efftor;
              std::vector<double> 
                  bnd(3), bndtl(3), endefftor(3), anchrs (3);

              hd = pChnFlipCln->getResidue ( 0 )->getAtom ( PID::C )->getPos();
              tl = pChnFlipCln->getResidue ( 0 )->getAtom ( PID::C_ALPHA )->getPos();
              bnd[0] = hd.x - tl.x; 
              bnd[1] = hd.y - tl.y; 
              bnd[2] = hd.z - tl.z;
              bndtl[0] = tl.x; 
              bndtl[1] = tl.y; 
              bndtl[2] = tl.z;

              efftor = pChnCln->getResidue(KNMTCS_SIZE)->getAtom(PID::N)->getPos();
              anchrs[0] = efftor.x; 
              anchrs[1] = efftor.y; 
              anchrs[2] = efftor.z;

              efftor = pChnFlipCln->getResidue ( 0 )->getAtom ( PID::N )->getPos();
              endefftor[0] = efftor.x; 
              endefftor[1] = efftor.y; 
              endefftor[2] = efftor.z;

              const Real dAngle = CalcMinimizingAngle(
                      bnd, bndtl, endefftor, anchrs
                      ) * 180. / M_PI;
              pChnFlipCln->RotateBackbone(1 /*psi*/, BondDirection(0), dAngle);
              std::cout << dAngle << "\n";

              std::ostringstream oss;
              oss << k;
              k++;
              PDBIO::writeToFile(pChnFlipCln, "afterFlip_" + oss.str() + ".pdb");

              GenTrialPositions( pChnFlipCln->getResidue( 0 ), Atm_lst );

              pChnFlipCln->Obliterate(); /*deallocates pChnSmall too */
          }
      }
      pChnCln->Obliterate();
  }

  pChn->Obliterate();
  std::cout << "SampleCB_MC: Conformers to estimate: " << Atm_lst.size() << std::endl;
  return true;
}


// Unused function
bool CSideChainSampler::Sample_OnMC(
        PProteinResidue* pPres, 
        const std::vector<double>& oc,  
        std::vector<clipper::Atom_list>& Atm_lst
        )
{
  /* 1. Extract coords from Atm_lst
   * 2. Sample new chi
   */
  const int KNMTCS_SIZE = 3;
  std::vector<clipper::Atom_list> atm_lst_;

  PProtein *pChn = pPres->getProtein()->Clone();

  double occ_thrshld = 0.09;
  if ( find_if (oc.begin(), oc.end(),bind2nd(greater<double>(),occ_thrshld) ) == oc.end() )
  {
    std::cout << "WARNING: minimum occupancy not found." << std::endl;
    occ_thrshld = 0.0;
  }

  for ( int j=0; j<oc.size(); ++j )
    if ( oc[j] > occ_thrshld )
      atm_lst_.push_back ( Atm_lst[j] );

  Atm_lst.clear();

  /* Output below looks funky if MC atoms came from flipped peptide.
   *Shouldn't matter as long as 'N' of next atom is not part of fit
   */

  for ( int j=0; j<atm_lst_.size(); ++j )
  {
    SetAtoms ( pChn->getResidue( KNMTCS_SIZE ), atm_lst_[j] );
    GenTrialPositions( pChn->getResidue( KNMTCS_SIZE ), Atm_lst );
    std::ostringstream oss;
    oss << j;
    PDBIO::writeToFile ( pChn, "after_sampled"+oss.str()+".pdb" );
  }

  pChn->Obliterate();
  return true;
}


bool CSideChainSampler::Sample_NextChi(
        PProteinResidue* pPres, 
        const int next_chi_idx, 
        const std::vector<double>& oc,  
        std::vector<clipper::Atom_list>& Atm_lst
        )
{
    /* 1. Extract coords from Atm_lst
     * 2. Sample new chi
     */
    const int KNMTCS_SIZE = 3;
    std::vector<clipper::Atom_list> atm_lst_;
     
    PProtein *pChn = pPres->getProtein()->Clone();

    double occ_thrshld = 0.0;
    if ( find_if (oc.begin(), oc.end(),bind2nd(greater<double>(),occ_thrshld) ) == oc.end() )
    {
        std::cout << "WARNING: minimum occupancy not found." << std::endl;
        occ_thrshld = 0.0;
    }
          
    for ( int j=0; j<oc.size(); ++j ) {
        if ( oc[j] > occ_thrshld ) {
            atm_lst_.push_back ( Atm_lst[j] );
        }
    }
    Atm_lst.clear();
          
    //  std::cout << "atm_lst_[0]: " << atm_lst_[0].size() << std::endl;
    //  std::cout << "atm_ids: " << atm_ids.size() << std::endl;
          
    for ( int j=0; j<atm_lst_.size(); ++j ) {
        SetAtoms ( pChn->getResidue( KNMTCS_SIZE ), atm_lst_[j], next_chi_idx );
        GenTrialPositions_AtChi( pChn->getResidue( KNMTCS_SIZE ), next_chi_idx, Atm_lst );
    }
          
    pChn->Obliterate();
    return true;
}


bool CSideChainSampler::SampleCB_MC_AtChi(
        PProteinResidue* pPres, 
        std::vector<clipper::Atom_list>& Atm_lst
        )
{
    const int KNMTCS_SIZE = 3;

    std::vector<std::pair<Real, std::pair< Vector3, bool > > > mvDirs;

    for ( float ampl = 0.10; ampl < 0.31; ampl += 0.1 ) {
      MC_AMPL = ampl;
      mvDirs.push_back(
              std::make_pair(
                  0.0, std::make_pair(
                      (MC_AMPL + MC_SIG*((double)rand()/RAND_MAX)) * 
                          Vector3(Ueigen (0,0), Ueigen (1,0), Ueigen (2,0))
                      , false)
                  )
              );
      mvDirs.push_back ( std::make_pair ( 0.0, std::make_pair ( ( -MC_AMPL-MC_SIG*((double)rand()/RAND_MAX) )*Vector3 ( Ueigen (0,0), Ueigen (1,0), Ueigen (2,0) ), false ) ) );
      mvDirs.push_back ( std::make_pair ( 0.0, std::make_pair ( ( MC_AMPL+MC_SIG*((double)rand()/RAND_MAX) )*Vector3 ( Ueigen (0,1), Ueigen (1,1), Ueigen (2,1) ), false ) ) );
      mvDirs.push_back ( std::make_pair ( 0.0, std::make_pair ( ( -MC_AMPL-MC_SIG*((double)rand()/RAND_MAX) )*Vector3 ( Ueigen (0,1), Ueigen (1,1), Ueigen (2,1) ), false ) ) );
      mvDirs.push_back ( std::make_pair ( 0.0, std::make_pair ( ( MC_AMPL+MC_SIG*((double)rand()/RAND_MAX) )*Vector3 ( Ueigen (0,2), Ueigen (1,2), Ueigen (2,2) ), false ) ) );
      mvDirs.push_back ( std::make_pair ( 0.0, std::make_pair ( ( -MC_AMPL-MC_SIG*((double)rand()/RAND_MAX) )*Vector3 ( Ueigen (0,2), Ueigen (1,2), Ueigen (2,2) ), false ) ) );
      mvDirs.push_back ( std::make_pair ( 0.0, std::make_pair ( ( MC_AMPL+MC_SIG*((double)rand()/RAND_MAX) )*Vector3 ( EVectorSum[0], EVectorSum[1], EVectorSum[2] ), false ) ) );
      mvDirs.push_back ( std::make_pair ( 0.0, std::make_pair ( ( -MC_AMPL-MC_SIG*((double)rand()/RAND_MAX) )*Vector3 ( EVectorSum[0], EVectorSum[1], EVectorSum[2] ), false ) ) );
    }

    PProtein *pChn = pPres->getProtein()->Clone();
    int chi_idx = 1;

    const int nbb = mvDirs.size();
    const int nflips = 4;
    int k = 0;
    for ( int j=-1; j<nbb; ++j) {
      PProtein *pChnCln = pChn->Clone();

      if ( j>= 0 ) {
        IKSolutions myIKSolns = MoveAtom(
                pChnCln, 
                pChnCln->getResidue(KNMTCS_SIZE)->getAtom(PID::C_BETA), 
                mvDirs[j].second.first
                );
      }

      std::ostringstream oss;
      oss << k;
      k++;
      PDBIO::writeToFile(pChnCln, "after_" + oss.str() + ".pdb" );

      GenTrialPositions_AtChi(pChnCln->getResidue(KNMTCS_SIZE), chi_idx, Atm_lst);

      if ( flip_peptide && ( j == -1 || j == 6 || j == 7 ) )
          for ( int l=0; l<nflips; ++l )
          {
              PProtein *pChnFlipCln = pChnCln->Clone ( );
              if ( l == 0 )
                  FlipTransformPeptide ( pChnFlipCln, FD0, 3 );
              else if ( l == 1 )
                  FlipTransformPeptide ( pChnFlipCln, FD1, 3 );
              else if ( l == 2 )
                  FlipTransformPeptide ( pChnFlipCln, FD2, 3 );
              else if ( l == 3 )
                  FlipTransformPeptide ( pChnFlipCln, FD3, 3 );

              PProtein *pChnSmall = new PProtein ( pChnFlipCln, KNMTCS_SIZE,  KNMTCS_SIZE );

              PDBIO::writeMainchainToFile ( pChnSmall, ".t.pdb" );
              pChnFlipCln->Obliterate();
              // Rereading in the main chain builds automatically the side chains.
              pChnFlipCln = PDBIO::readFromFile( ".t.pdb" );
              PDBIO::writeToFile ( pChnFlipCln, ".tsc.pdb" );

              Vector3 hd, tl, efftor;
              std::vector<double> bnd ( 3 ), bndtl ( 3 ), endefftor ( 3 ), anchrs ( 3 );

              hd = pChnFlipCln->getResidue ( 0 )->getAtom ( PID::C )->getPos();
              tl = pChnFlipCln->getResidue ( 0 )->getAtom ( PID::C_ALPHA )->getPos();
              bnd[0] = hd.x - tl.x; 
              bnd[1] = hd.y-tl.y; 
              bnd[2] = hd.z-tl.z;
              bndtl[0] = tl.x; 
              bndtl[1] = tl.y; 
              bndtl[2] = tl.z;

              efftor = pChnCln->getResidue ( KNMTCS_SIZE )->getAtom ( PID::N )->getPos();
              anchrs[0] = efftor.x; 
              anchrs[1] = efftor.y; 
              anchrs[2] = efftor.z;

              efftor = pChnFlipCln->getResidue ( 0 )->getAtom ( PID::N )->getPos();
              endefftor[0] = efftor.x; 
              endefftor[1] = efftor.y; 
              endefftor[2] = efftor.z;

              const Real dAngle = CalcMinimizingAngle ( bnd, bndtl, endefftor, anchrs ) * 180./M_PI;
              pChnFlipCln->RotateBackbone ( 1 /*psi*/, BondDirection(0), dAngle );
              std::cout << dAngle << "\n";

              std::ostringstream oss;
              oss << k;
              k++;
              PDBIO::writeToFile ( pChnFlipCln, "afterFlip_"+oss.str()+".pdb" );

              GenTrialPositions_AtChi( pChnFlipCln->getResidue( 0 ), chi_idx, Atm_lst );

              pChnFlipCln->Obliterate(); /*deallocates pChnSmall too */
          }
      pChnCln->Obliterate();
    }
    pChn->Obliterate();
    return true;
}


bool CSideChainSampler::FlipPeptide(
        PProtein* pChn, const FlipDirection isForw, 
        const float dAngle
        )
{
   const int KNMTCS_SIZE = 3;
   PAtom *startV, *endV;
   Vector3 axis;

   const int start = (isForw == FlipForward) ? 0 : -1;

   PResidue* pPres = pChn->getResidue ( KNMTCS_SIZE + start );
   PResidue* pPresp1 = pChn->getResidue ( KNMTCS_SIZE + start + 1 );

   startV = pPres->getAtom ( PID::C_ALPHA );
   endV = pPresp1->getAtom ( PID::C_ALPHA );
   axis = endV->getPos()-startV->getPos();

   PBond *pBond = pPresp1->getBond ( PID::N, PID::C_ALPHA );

   MyRotater rinit(startV->getPos(),axis, dAngle);

   pBond->traverseAtoms(backward, &rinit, pChn, pPres->getAtom ( PID::C_ALPHA) );

   Matrix3 M = PMath::FindRotationMatrix(endV->getPos()-startV->getPos(),DtoR(dAngle));
   Vector3 P = startV->getPos()+M*(pPresp1->getAtom ( PID::N )->getPos()-startV->getPos());
   pPresp1->getAtom ( PID::N )->changePosition(P);

   return true;
}


/* Sample CB position by rotating along phi (forward) and psi (backward).
 * This method deforms the N-CA-C-CB tetrahedron.
 */
// Unused function
bool CSideChainSampler::SampleCB(
        PProteinResidue* pPres, 
        std::vector<clipper::Atom_list>& Atm_lst
        )
{
  if ( pPres->getName()==PID::GLY ) {
    Atm_lst.push_back ( GetAtoms ( pPres, true ) );
    return true;
  }
  
  const float dAngle = 5.0;
  PAtom *startV, *endV;
  Vector3 axis[2];
  vector<PProtein*> temp_Lps;
  
  startV = pPres->getAtom ( PID::C_ALPHA );
  endV = pPres->getAtom ( PID::N );
  axis[0] = endV->getPos()-startV->getPos();
  axis[1] = pPres->getAtom ( PID::C ) - startV;
  
  PBond *pBond = pPres->getBond ( PID::C_ALPHA, PID::C_BETA );
    
  MyRotater rinit1(startV->getPos(),axis[0], -5.0);
  MyRotater rinit2(startV->getPos(),axis[1], -5.0);
  pBond->traverseAtoms(forward, &rinit1, pPres->getChain());
  pBond->traverseAtoms(forward, &rinit2, pPres->getChain());
  
  for ( int d1=0; d1<3; d1++ ) {
    for ( int d2=0; d2<3; d2++ ) {
      GenTrialPositions( pPres, Atm_lst );
      MyRotater r2(endV->getPos(),axis[1], dAngle);
      pBond->traverseAtoms(forward, &r2, pPres->getChain());
    }  
    MyRotater rr(endV->getPos(),axis[1], -3*dAngle);
    pBond->traverseAtoms(forward, &rr, pPres->getChain());
      
    MyRotater r1(endV->getPos(),axis[0], dAngle);
    pBond->traverseAtoms(forward, &r1, pPres->getChain());
  }
  return true;
}
  

void CSideChainSampler::RotamerApply(
        PResidue *res, const vector<Real> &rotamer
        )
{
  // see how many chi angles in this give residue.
  unsigned chiMax = PResources::numChiIndices(res->getName());
  // store the name of atoms that defines the chi angle.
  for (unsigned i = 1; i<=chiMax; i++)
    SetChi(res, i, rotamer[i-1]);
}


void CSideChainSampler::SetChi(
        PResidue *res, const int chi_idx, const Real chi
        )
{
  vector<string> rotList;
  Real curDihedral, toRotate;

  rotList = PResources::GetChiIndex(res->getName(), chi_idx );
  
  curDihedral = PMath::AngleBetweenPlanes(
          res->getAtomPosition(rotList[0]), 
          res->getAtomPosition(rotList[1]), 
          res->getAtomPosition(rotList[2]), 
          res->getAtomPosition(rotList[3])
          );
  toRotate = (chi - curDihedral);

  if (res->getName()==PID::PRO) {
    Matrix3 M = PMath::FindRotationMatrix(
            res->getAtomPosition(rotList[2]) -
            res->getAtomPosition(rotList[1]),
            -toRotate * M_PI / 180.0
            );
    Vector3 P = res->getAtomPosition(rotList[2]) + M * (
            res->getAtomPosition(rotList[3]) - 
            res->getAtomPosition(rotList[2])
            );
    res->getAtom(rotList[3])->changePosition(P);
  } else
    res->getDOF("sidechain", rotList[1], rotList[2])->Rotate(forward, (-toRotate));
}


void CSideChainSampler::SetTFactor(
        PResidue *res, const Real TempFactor
        )
{
  const Real dT = 0.09*TempFactor;
  //see how many chi angles in this give residue.
  unsigned chiMax = PResources::numChiIndices(res->getName());
  //store the name of atoms that defines the chi angle.
  vector<string> rotList;
  //amount to rotate.
  for(unsigned i = 1; i<=chiMax; i++) {
      MyTSetter mTSetter(TempFactor + i*dT + 
              0.09 * TempFactor * ((double) rand() / RAND_MAX - 0.5)
              );
      //get the atoms used to define dihedral angles
      rotList = PResources::GetChiIndex(res->getName(), i);
      //calculate current dihedral angle (chi angle) and calculate the delta
      //needed to achive the goal.
      res->getDOF("sidechain", rotList[1], rotList[2])->
          traverseAtoms(forward, &mTSetter, res->getChain());
  }
  res->getAtom(PID::C_BETA)->setTempFactor(
          TempFactor + dT*((double)rand()/RAND_MAX-0.5)
          );
  res->getAtom (PID::O)->setTempFactor(
          TempFactor + dT*((double)rand()/RAND_MAX-0.5)
          );
}


void CSideChainSampler::SetTFactor(
        PResidue *res, const Real TempFactor, const int chi_idx 
        )
{
  const Real dT = 0.09*TempFactor;
  //if (res->getName()==PID::PRO) return;
  //see how many chi angles in this give residue.
  unsigned chiMax = PResources::numChiIndices(res->getName());
  //store the name of atoms that defines the chi angle.
  vector<string> rotList;
  //amount to rotate.
  for(unsigned i = chi_idx; i<=chiMax; i++) {
    MyTSetter mTSetter(TempFactor + i*dT + 0.09*TempFactor*((double)rand()/RAND_MAX-0.5));
    //get the atoms used to define dihedral angles
    rotList = PResources::GetChiIndex(res->getName(), i);
    //calculate current dihedral angle (chi angle) and calculate the delta needed to achive the goal.
    res->getDOF("sidechain",rotList[1],rotList[2])->traverseAtoms(forward, &mTSetter, res->getChain());
  }
  if (chi_idx == 1)
  {
    res->getAtom (PID::C_BETA)->setTempFactor ( TempFactor + 0.09*TempFactor*((double)rand()/RAND_MAX-0.5) );
    res->getAtom (PID::O)->setTempFactor ( TempFactor + 0.09*TempFactor*((double)rand()/RAND_MAX-0.5) );
  }
}


void CSideChainSampler::SampleRotamerNeighborhood(
        const PResidue *res, 
        const vector<Real> &rotamer, 
        vector<vector<Real> >& rotamer_nbrhd
        )
{
  
  const double delta = 10.;
  const int SAMPLE_SIZE = 5;
  vector<Real> rotamer_delta[SAMPLE_SIZE];
  
  for ( int i=-2; i<3; i++ )
    rotamer_delta[i+2] = std::vector<Real> ( rotamer.size(), i*delta );

  for ( int i=0; i<SAMPLE_SIZE; i++ ) {
    for ( int j=0; j<rotamer.size(); ++j ) {
      rotamer_delta[i][j] += rotamer[j];
    }
  }

  //see how many chi angles in this give residue.
  unsigned chiMax = PResources::numChiIndices(res->getName());
  if ( chiMax > 4 )
    std::cout << "ERROR: chiMax > 4. Fix needed." << std::endl;
  
  rotamer_nbrhd.resize ( (int)pow ( (float)SAMPLE_SIZE, (float)chiMax ) );
  int p=0;
  for ( int i=0; i<SAMPLE_SIZE; i++ ) {
    if ( chiMax == 1 ) {
      rotamer_nbrhd[p].push_back ( rotamer_delta[i][0] );
      p++; 
      continue;
    }
    for ( int j=0; j<SAMPLE_SIZE; j++ ) {
      if ( chiMax == 2 ) {
        rotamer_nbrhd[p].push_back ( rotamer_delta[i][0] );
        rotamer_nbrhd[p].push_back ( rotamer_delta[j][1] );
        p++; 
        continue;
      }
      for ( int k=0; k<SAMPLE_SIZE; k++ )
      {
        if ( chiMax == 3 )
        {
          rotamer_nbrhd[p].push_back ( rotamer_delta[i][0] );
          rotamer_nbrhd[p].push_back ( rotamer_delta[j][1] );
          rotamer_nbrhd[p].push_back ( rotamer_delta[k][2] );
          p++; 
          continue;
        }
        for ( int l=0; l<SAMPLE_SIZE; l++ )
        {
          rotamer_nbrhd[p].push_back ( rotamer_delta[i][0] );
          rotamer_nbrhd[p].push_back ( rotamer_delta[j][1] );
          rotamer_nbrhd[p].push_back ( rotamer_delta[k][2] );
          rotamer_nbrhd[p].push_back ( rotamer_delta[l][3] );
          p++;
        }
      }
    }
  }
}


/***** PEPTIDE FLIPS *****/

struct Transformer: AtomFunctor {
public:

  Transformer(Matrix4 toApply)
  {
    m_toApply = toApply;
    cerr<<toApply<<endl;
  }

  void operator()(
          PAtom *atom, 
          PBond *bondFrom
          ) 
  {
    atom->ApplyTransform(m_toApply);
  }

private:

  Matrix4 m_toApply;

};


void CSideChainSampler::MakeHTMatrix(
        const Matrix3& rFrame, 
        const Vector3& org, 
        Matrix4& htMatrix
        )
{
    for ( int i=0;i<3;++i )
        for ( int j=0;j<3;++j )
            htMatrix.data[i][j] = rFrame.data[i][j];

    Vector3 tmp = rFrame*org;

    tmp *= -1;

    for ( int i=0; i<3; ++i)
       htMatrix.data[3][i] = tmp[i];
    htMatrix.data[3][3] = 1;

    for ( int i=0; i<3; ++i )
        htMatrix.data[i][3] = 0;
}


/*Takes atom positions of CA, CA+1 , and O */
void CSideChainSampler::CoordFrame(
        const Vector3& a,  
        const Vector3& b, 
        const Vector3& c, 
        Matrix3& T, 
        Vector3& org
        )
{
    Vector3 xa, xb, xc;
    org = a;
    xa = b - org;
    xb = c - org;

    Real gs = dot ( xb, xa ) / dot ( xa, xa);

    for ( int i=0; i<3; ++i )
      xb[i] -= gs * xa[i];

    normalize ( xa );

    normalize ( xb );

    xc = cross ( xa, xb );

    normalize(xc);

    T = Matrix3 (xa, xb, xc);

    T.inplaceTranspose();

}


bool CSideChainSampler::PepCoordFrame(
        PProtein* pChn, 
        Matrix3& T, 
        Vector3& org, 
        const int resno
        )
{
    const int KNMTCS_SIZE = resno;
    PAtom *startV, *endV, *carbonylO;

    PResidue* pPres = pChn->getResidue ( KNMTCS_SIZE );
    PResidue* pPresp1 = pChn->getResidue ( KNMTCS_SIZE + 1 );

    startV = pPres->getAtom ( PID::C_ALPHA );
    endV = pPresp1->getAtom ( PID::C_ALPHA );
    carbonylO = pPres->getAtom ( PID::O );

    CoordFrame(startV->getPos(), endV->getPos(), carbonylO->getPos(), T, org );

    return true;
}


bool CSideChainSampler::FlipTransformPeptide(
        PProtein* pChn, 
        const double (*FD)[4], 
        const int resno
        )
{
    const int KNMTCS_SIZE = resno;
    Vector3 refx ( ReferencePeptide[0] );
    Vector3 refy ( ReferencePeptide[1] );
    Vector3 refz ( ReferencePeptide[2] );

    /*Compute Transform from Origin to Reference */

    Vector3 org_o2r (0,0,0);
    Matrix3 R_o2r;
    Matrix4 T_o2r, T_o2rInv;

    CoordFrame ( refx, refy, refz, R_o2r, org_o2r );
    MakeHTMatrix ( R_o2r, org_o2r, T_o2r );
    T_o2r.getInverse(T_o2rInv);

    /*Compute Transform from Origin to Reference */

    Vector3 org_t2o (0,0,0);
    Matrix3 R_t2o;
    Matrix4 T_t2o, T_t2oInv;

    PepCoordFrame ( pChn, R_t2o, org_t2o, resno );
    MakeHTMatrix ( R_t2o, org_t2o, T_t2o);
    T_t2o.getInverse(T_t2oInv);

    Matrix4 F1 (FD);
    F1.inplaceTranspose();

    Matrix4 T = T_t2oInv * T_o2r * F1 * T_o2rInv * T_t2o;

    vector<PAtom *> vAtoms;
    vAtoms.push_back ( pChn->getResidue(KNMTCS_SIZE)->getAtom ( PID::N));
    vAtoms.push_back ( pChn->getResidue(KNMTCS_SIZE)->getAtom ( PID::C_ALPHA));
    vAtoms.push_back ( pChn->getResidue(KNMTCS_SIZE)->getAtom ( PID::C));
    vAtoms.push_back ( pChn->getResidue(KNMTCS_SIZE)->getAtom ( PID::O));
    vAtoms.push_back ( pChn->getResidue(KNMTCS_SIZE+1)->getAtom ( PID::N));

    for ( int i=0; i<vAtoms.size(); ++i )
    {
        vAtoms[i]->ApplyTransform(T);
    }

//    PBond* pB = pChn->getResidue(KNMTCS_SIZE)->getBond ( PID::C_ALPHA, PID::C_BETA );
//    Transformer Tformer ( T );
//    pB->traverseChain ( forward, &Tformer, NULL, pChn );

    return true;
}
/***** end PEPTIDE FLIPS *****/
