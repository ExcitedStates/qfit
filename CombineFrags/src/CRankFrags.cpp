#include "CRankFrags.h"

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
    std::cout << "Error: Can not open file " << fn << " for writing." << std::endl;
  return success;
}

CRankFrags::CRankFrags() : mainCorrHandle (-1), sideCorrHandle(-1), registerUD(false)
{
}

bool CRankFrags::CalcMap (  const std::string& mtzFile )
{
  //std::cout << "Calculating map for: " << mtzFile.c_str() << "\n";
  const char* F_phi_strings[4] = { "*/*/[2FOFCWT,PH2FOFCWT]", "*/*/[FSHFT,PHISHFT]", "*/*/[FP,PHIM]", "*/*/[FC,PHIC]" };
  const int F_phi_trials =4;
  
     
  clipper::MTZcrystal xtal;
  clipper::MTZdataset dset;
  clipper::CCP4MTZfile file;

  file.open_read ( mtzFile );
  hkl_info = clipper::HKL_info ( );
  //hkl_info.init ( file.spacegroup(), file.cell(), file.resolution ( ) );
  
  clipper::HKL_data<clipper::data32::F_phi> f_phi ( hkl_info );
  
//  clipper::HKL_data<clipper::data32::Phi_fom> phi_fom ( hkl_info );
  file.import_hkl_info( hkl_info, true );
  
//  file.import_chkl_data( f_phi, "*/*/[FP,PHIM]" );

  mapres = hkl_info.resolution().limit();
  RADIUS = mapres < 3.0 ? 0.7 + ( mapres - 0.6 )/3.0 : 0.5 * mapres;
  BSMEAR = 5.6595*pow ( mapres, 2.315 );

  bool F_phi_found = false;
  int j = 0; 
  while (j < F_phi_trials && !F_phi_found )
  {
    try
    {
      //std::cout << "Trying: " << F_phi_strings[j] << "\n";
      file.import_hkl_data( f_phi, F_phi_strings[j] );
      //std::cout << F_phi_strings[j] << "found!\n";
      F_phi_found = true;;
    }
    catch ( clipper::Message_fatal msg ) 
    { 
      //std::cout << msg.text() << "\n";
      j++;
    }
  }
  if ( !F_phi_found  )
  {
    std::cout << "No Amplitude/Phase info found in file: " << mtzFile.c_str() << "\n";
    return false;
  }
  
//  file.import_hkl_data( f_sigf, dset, xtal, "*/*/[FP,SIGFP]" );
//  file.import_chkl_data( phi_fom, "*/*/[PHIM,FOMM]" );

  file.close_read();
  

  clipper::Grid_sampling grid( file.spacegroup(), file.cell(), file.resolution ( ) );
  xmap.init ( file.spacegroup(), file.cell(), grid );
  xmap.fft_from( f_phi, clipper::Xmap_base::Normal );
//  clipper::CCP4MAPfile mapout;
//  mapout.open_write( "out1FOMM.map" );
//  mapout.export_xmap( xmap );
//  mapout.close_write();

  return true;
}

bool CRankFrags::Register ( PCMMDBManager pCMMDBManager )
{
  mainCorrHandle=pCMMDBManager->RegisterUDReal(UDR_RESIDUE,"mainCoor");
  sideCorrHandle=pCMMDBManager->RegisterUDReal(UDR_RESIDUE,"sideCoor");
  registerUD = true;
  return registerUD;
}

int CRankFrags::Rankum ( const PCMMDBManager pCMMDBManager )
{
  PPCResidue pResTable;
  int nRsds;
  int best = 0;
  double best_cc, total_cc, cc;
  best_cc = -1;
  for ( int i=1; i<=pCMMDBManager->GetNumberOfModels(); i++ )
  {
    total_cc = 0;
    pCMMDBManager->GetModel(i)->GetResidueTable ( 0, pResTable, nRsds );
    if ( pResTable )
    {
      for ( int j=0; j<nRsds; j++ )
      {
        pResTable[j]->GetUDData ( mainCorrHandle, cc );
        total_cc += cc;
      }
      if ( best_cc < total_cc )
      {
        best = i;
        best_cc = total_cc;
      }
    }
  }
  return best;
}

bool CRankFrags::InsertBestTerminal ( PCMMDBManager pMMDBManager, PCMMDBManager pFragManager, const bool isForward )
{
  PPCResidue pResTable;
  int nRsds;
  int best = 0;
  double best_cc, total_cc, cc;
  best_cc = -1;
  nRsds = pFragManager->GetModel(1)->GetChain(0)->GetNumberOfResidues();
  const int start = isForward ? 0 : nRsds-3;
  const int end = isForward ? 3 : nRsds;
  
  for ( int i=1; i<=pFragManager->GetNumberOfModels(); i++ )
  {
    pFragManager->GetModel(i)->GetResidueTable ( 0, pResTable, nRsds );
    total_cc = 0;
    if ( pResTable )
    {
      for ( int j=start; j<end; j++ )
      {
        pResTable[j]->GetUDData ( mainCorrHandle, cc );
        total_cc += cc;
      }
      std::cout << "Total cc = " << total_cc << std::endl;
      if ( best_cc < total_cc )
      {
        best = i;
        best_cc = total_cc;
      }
    }
  }
  CMMDBManager bestManager;
  CModel mMdl;
  bestManager.AddModel ( &mMdl );
  CChain mChn;
  mMdl.AddChain ( &mChn );
  for ( int i=start; i<end; i++ )
    mChn.AddResidue ( pFragManager->GetModel ( best )->GetChain ( 0 )->GetResidue ( i ) );
  mChn.SetChainID ( pFragManager->GetModel ( best )->GetChain ( 0 )->GetChainID ( ) );
  bestManager.FinishStructEdit ( );
  std::cout << "Adding residues from fragment " << best << std::endl;
  InsertFragment ( pMMDBManager->GetModel ( 1 )->GetChain ( mChn.GetChainID ( ) ), &mChn );
  pMMDBManager->FinishStructEdit ( );
  if ( isForward )
  {
  	mChn.DeleteResidue ( end - 1 );
    WritePDB ( &mChn, "bestN.pdb" );
  }
  else
  {
  	mChn.DeleteResidue ( 0 );
    WritePDB ( &mChn, "bestC.pdb" );
  }  
  return true;
}

bool CRankFrags::InsertFragment ( PCChain pChn, const PCChain pFrag, const double ccTHRESHOLD )
{
  std::ofstream ccfile ( "ccbest.log"); 
  if (!ccfile) // were there any errors on opening?  
    std::cout << "ERROR: Unable to open file ccbest.log." << std::endl;
  ccfile.setf(std::ios::fixed,std::ios::floatfield); 
  ccfile.precision ( 2 );
  const int gapStartSeqNum = pFrag->GetResidue ( 0 )->GetSeqNum ( );
  const int gapEndSeqNum = gapStartSeqNum + pFrag->GetNumberOfResidues( );
  int resno = pChn->GetResidueNo ( gapStartSeqNum, "" );
  if ( resno == -1 )
    resno = pChn->GetResidueNo ( gapEndSeqNum, "" ) - 1;
  double cc;
  for ( int l=gapStartSeqNum; l<gapEndSeqNum; l++ )
  {
  pChn->DeleteResidue ( l, "");  
  pFrag->GetResidue ( l-gapStartSeqNum )->GetUDData ( mainCorrHandle, cc );
  if ( cc > ccTHRESHOLD ) 
  {
    PCResidue pRes = newCResidue ( );
    ccfile << pChn->GetChainID() << l << " ";
    ccfile << pFrag->GetResidue ( l-gapStartSeqNum )->GetResName ( ) << " ";
    ccfile << cc << std::endl;
    pRes->Copy ( pFrag->GetResidue ( l-gapStartSeqNum ) );
    pRes->SetResID ( pFrag->GetResidue ( l-gapStartSeqNum )->GetResName(), l, "" );
    pRes->SetChain ( pChn );
    pChn->InsResidue ( pRes, resno );
  }
  resno++;
  }
//  for ( int l=gapStartSeqNum; l<gapEndSeqNum; l++ )
//  {
//    PCAtom patm[4];
//    for ( int m=0; m<4; m++ )
//      patm[m] = pChn->GetResidue ( l-m-3, "" )->GetAtom ("CA");
//    rvector vctrs[3];
//    for ( int m=0; m<3; m++ )
//      GetVectorMemory ( vctrs[m], 3, 0 );
//    vctrs[0][0] = patm[0]->x - patm[1]->x; 
//    vctrs[0][1] = patm[0]->y - patm[1]->y;
//    vctrs[0][2] = patm[0]->z - patm[1]->z;
//    
//    vctrs[1][0] = patm[2]->x - patm[1]->x; 
//    vctrs[1][1] = patm[2]->y - patm[1]->y;
//    vctrs[1][2] = patm[2]->z - patm[1]->z;
//
//    vctrs[2][0] = patm[3]->x - patm[2]->x; 
//    vctrs[2][1] = patm[3]->y - patm[2]->y;
//    vctrs[2][2] = patm[3]->z - patm[2]->z;
//
//    std::cout << clipper::Util::rad2d ( GetAngle ( vctrs[0], vctrs[1] ) ) << " " << 
//      clipper::Util::rad2d ( GetTorsion ( vctrs[0], vctrs[1], vctrs[2] ) )<< std::endl;
//    
//    for ( int m=0; m<3; m++ )
//      FreeVectorMemory ( vctrs[m], 0 );
//  }
//  std::cout << std::endl;
  return true;
}

bool CRankFrags::CreateSkeleton ( const clipper::Xmap<float>& xmap )
{
  xskl.init ( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling( ) );
  xskl = 1;
  clipper::Skeleton_basic mSkel;
  
  mSkel (xskl, xmap); 
  
  clipper::CCP4MAPfile xmapfile;
  xmapfile.open_write( "skel.map" );
  xmapfile.export_xmap( xskl );
  xmapfile.close_write();
  return true;
}

void CRankFrags::Normalize ( clipper::Xmap<clipper::ftype32>& xm,  const clipper::Atom_list& atmList )
{
	/*clipper::Map_stats mstats( xm );
	
	float m = mstats.mean ( );
	float s = mstats.std_dev ( );

	clipper::Xmap_base::Map_reference_index ix;
	
	for ( ix = xm.first(); !ix.last(); ix.next() )
	{
		xm[ ix ] -= m;
		xm[ ix ] /= s;
	}*/

  float m, s;
  std::vector<clipper::Xmap<float>::Map_reference_coord> ref_coords = CreateCarve ( xm, atmList );

   m = s = 0.0;
   for ( int i=0; i<ref_coords.size (); i++ )
   {
      clipper::Xmap<float>::Map_reference_coord iw = ref_coords[i];
      m += xm[iw];
   }

   m /= ref_coords.size ();
   
   std::cout << "m = " << m << std::endl; 

   for ( int i=0; i<ref_coords.size (); i++ )
   {
      clipper::Xmap<float>::Map_reference_coord iw = ref_coords[i];
      s += ( xm[iw] - m )*( xm[iw] - m );
   }
   s /= ref_coords.size ();
   std::cout << "s = " << s << std::endl;
   for ( int i=0; i<ref_coords.size (); i++ )
   {
      clipper::Xmap<float>::Map_reference_coord iw = ref_coords[i];
      xm[iw] -= m;
      xm[iw] /= s;
   }
}

std::vector<clipper::Xmap<clipper::ftype32>::Map_reference_coord> CRankFrags::CreateCarve ( const clipper::Xmap<clipper::ftype32>& xm, 
		const clipper::Atom_list& atmList ) const
{
   std::vector<clipper::Xmap<clipper::ftype32>::Map_reference_coord> ref_coords;
	std::vector<int> index;
   clipper::Xmap<unsigned char> bexmap( xm.spacegroup(), xm.cell(), xm.grid_sampling() );
   bexmap = 0;
   ref_coords.reserve ( 5000 );
   index.reserve ( 5000 );
   
   const clipper::Grid_range gr( xm.cell(), xm.grid_sampling ( ), 2.5 );

   clipper::Xmap<clipper::ftype32>::Map_reference_coord i0, iu, iv, iw;

   for ( int i=0; i<atmList.size(); i++ )
   {
      const clipper::Coord_grid g0 = xm.coord_map( atmList[i].coord_orth() ).coord_grid() + gr.min();
      const clipper::Coord_grid g1 = xm.coord_map( atmList[i].coord_orth() ).coord_grid() + gr.max();
                
      i0 = clipper::Xmap_base::Map_reference_coord ( xm, g0 );
      for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
         for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
            for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
			   if ( bexmap[iw] == 0 )
			   {
				   ref_coords.push_back ( iw );
                   bexmap[iw] = 1;
               }
   }

   return ref_coords;
}

bool CRankFrags::calcEDMCC( const PCMMDBManager pCMMDBManager, const int nMdl )
{
  if ( !registerUD )
    Register ( pCMMDBManager );
  
  PCModel pModel = pCMMDBManager->GetModel ( nMdl );
//  clipper::CCP4MAPfile mapout;
//    mapout.open_write( "out1FOMM.map" );
//    mapout.export_xmap( xmap );
//    mapout.close_write();

  // now test ed calc
  //
  // get a list of all the atoms
  clipper::mmdb::PPCAtom psel;
  int hndl, nsel;
  hndl = pCMMDBManager->NewSelection();
          pCMMDBManager->SelectAtoms ( hndl,nMdl,   
          "*", 
          ANY_RES,"*",
          ANY_RES,"*", 
          "*", "*", "*", "*",  SKEY_NEW  );
  pCMMDBManager->GetSelIndex( hndl, psel, nsel );
  clipper::MMDBAtom_list atoms( psel, nsel );
  pCMMDBManager->DeleteSelection( hndl );
  
  clipper::Xmap<float> exmap( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() );
  clipper::Xmap<unsigned char> bexmap( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() );
//  clipper::EDcalc_iso<float> ediso;
//  ediso( exmap, atoms );
//  mapout.open_write( "out_1.map" );
//  mapout.export_xmap( exmap );
//  mapout.close_write();
 
  clipper::HKL_data<clipper::data32::F_phi> f_phi_calc ( hkl_info );
  clipper::SFcalc_iso_fft<float> sfc;
  sfc ( f_phi_calc, atoms );
/*  clipper::SFcalc_obs_bulk<float> sfcb;
  sfcb (f_phi_calc, f_sigf, atoms );
  std::cout <<  "Bulk Correction Volume: " << sfcb.bulk_frac() << " Bulk Correction Factor: " << sfcb.bulk_scale() << std::endl;
  clipper::SFscale_aniso<float> sc;
  sc ( f_phi_calc, f_sigf );
  std::cout << "\nAnisotropic scaling:\n"
              << sc.u_aniso_orth().format() << "\n";
*/  
  exmap.fft_from( f_phi_calc, clipper::Xmap_base::Normal );
  
  
//  mapout.open_write( "out_2.map" );
//  mapout.export_xmap( exmap );
//  mapout.close_write();
  
  float mo, mc, sso, ssc, ssoc, st;
  std::vector<clipper::Xmap<float>::Map_reference_coord> ref_coords = CreateCarve ( xmap, atoms );

   mo = mc = 0.0;
   for ( int i=0; i<ref_coords.size (); i++ )
   {
      clipper::Xmap<float>::Map_reference_coord iw = ref_coords[i];
      mo += xmap[iw];
      mc += exmap[iw];
   }

   mo /= ref_coords.size ();
   mc /= ref_coords.size ();

   sso = ssc = ssoc = 0.0;

   for ( int i=0; i<ref_coords.size (); i++ )
   {
      clipper::Xmap<float>::Map_reference_coord iw = ref_coords[i];
      sso += (xmap[iw] - mo)*(xmap[iw] - mo);
      ssc += (exmap[iw] - mc)*(exmap[iw] - mc);
      ssoc += (exmap[iw] - mo)*(exmap[iw] - mc);
	}
	
   st = (ssoc*ssoc);
   st /= (sso*ssc);
  
//  std::cout << st << std::endl;
  
 // Normalize ( xmap, atoms );
 // Normalize ( exmap, atoms );

//  std::cout << xmap.cell().format () << " " << xmap.grid_sampling().format() << std::endl;

  // work out how big a box we need to calc density over for each atom
  clipper::Grid_sampling grid( hkl_info.spacegroup(), hkl_info.cell(), hkl_info.resolution ( ) );
  clipper::Grid_range gd( hkl_info.cell(), xmap.grid_sampling(), RADIUS );
  clipper::Xmap<float>::Map_reference_coord i0, iu, iv, iw;
  std::cout << hkl_info.cell().format () << " " << grid.format() << std::endl;

  clipper::mmdb::PPCChain pChnTable;
  clipper::mmdb::PPCResidue pRsdTable;
  clipper::mmdb::PPCAtom pAtmTable;
  int nChns, nRsds, nAtms;

  pModel->GetChainTable ( pChnTable,nChns );
  char S[128];
  float x,xx,xy,y,yy,n;
    for ( int chn_idx=0; chn_idx<nChns; chn_idx++ )
    {
      pChnTable[chn_idx]->GetResidueTable ( pRsdTable, nRsds );
      for ( int rsd_idx=0; rsd_idx<nRsds; rsd_idx++ )
      {
        float cc[2] = {0.0,0.0};
        //std::cout << pRsdTable[rsd_idx]->GetResidueID ( S ) << " ";
        int selHnd;
        selHnd = pCMMDBManager->NewSelection();    
        pCMMDBManager->SelectAtoms ( selHnd,0,   
          pChnTable[chn_idx]->GetChainID(), 
          pRsdTable[rsd_idx]->GetSeqNum(),"*",
          pRsdTable[rsd_idx]->GetSeqNum(),"*", 
          "*", "CA", "C", "*",  SKEY_NEW  );
        pCMMDBManager->SelectAtoms ( selHnd,0,   
          pChnTable[chn_idx]->GetChainID(),
          pRsdTable[rsd_idx]->GetSeqNum(),"*",
          pRsdTable[rsd_idx]->GetSeqNum(),"*",
          "*", "*", "*", "*", SKEY_OR );
        for ( int mnsd=0; mnsd<2; mnsd++ )
        { 
          bexmap = 0;
          x=xx=y=yy=xy=n=0;
          if ( mnsd==1 )
          {
            pCMMDBManager->SelectAtoms ( selHnd,0,   
            pChnTable[chn_idx]->GetChainID(),
            pRsdTable[rsd_idx]->GetSeqNum(),"*",
            pRsdTable[rsd_idx]->GetSeqNum(),"*",
            "*", "*", "*", "*", SKEY_XOR );
          }
          pCMMDBManager->GetSelIndex( selHnd, pAtmTable, nAtms );
//          char S[128];
//          for ( int kk=0; kk<nAtms; kk++ )
//            std::cout << pAtmTable[kk]->GetAtomID ( S ) << std::endl;
          clipper::MMDBAtom_list atoms( pAtmTable,nAtms );

          for ( int i = 0; i < atoms.size(); i++ )
          {
            if ( !atoms[i].is_null() ) 
            {
              clipper::Coord_frac uvw = atoms[i].coord_orth().coord_frac( hkl_info.cell() );
              clipper::Coord_grid g0 = uvw.coord_grid( grid ) + (gd.min)();
              clipper::Coord_grid g1 = uvw.coord_grid( grid ) + (gd.max)();
              i0 = clipper::Xmap<float>::Map_reference_coord( xmap, g0 );
              // sum all map contributions from this atoms
              for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
                for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
                  for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
                    if ( bexmap[iw]==0 )
                    {
                      bexmap[iw] = 1;
                      x += xmap[iw];
                      xx += xmap[iw]*xmap[iw];
                      yy += exmap[iw]*exmap[iw];
                      y += exmap[iw];
                      xy += xmap[iw]*exmap[iw];
                      n++;
                    }
            }
          }
          if ( atoms.size() > 0 )
            cc[mnsd] = (xy/n-x*y/(n*n))/(sqrt(xx/n-(x/n)*(x/n))*sqrt(yy/n-(y/n)*(y/n)));
	  cc[mnsd]*=cc[mnsd];
      }
      pRsdTable[rsd_idx]->PutUDData(mainCorrHandle,cc[0]);
      pRsdTable[rsd_idx]->PutUDData(sideCorrHandle,cc[1]);      
//      std::cout << pRsdTable[rsd_idx]->GetChainID() << pRsdTable[rsd_idx]->GetSeqNum() << " "<< cc[0] << " " << cc[1] << "\n";
      pCMMDBManager->DeleteSelection(selHnd);
    }
  }
  std::cout << std::endl;
  return true;
}

CRankFrags::~CRankFrags()
{
}
