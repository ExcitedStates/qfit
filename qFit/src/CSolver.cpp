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


#ifndef CSOLVER_H_
#include "CSolver.h"
#endif

#include <CGAL/MP_Float.h>
#include <CGAL/Timer.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
//#include <gsl/gsl_linalg.h>

#include "MyMIQP.hpp"

typedef CGAL::MP_Float ET;

typedef CGAL::Quadratic_program_from_iterators
    <double**,                                                // for A
     double*,                                                 // for b
     CGAL::Const_oneset_iterator<CGAL::Comparison_result>, // for r
     bool*,                                                // for fl
     double*,                                                 // for l
     bool*,                                                // for fu
     double*,                                                 // for u
     double**,                                                // for D
     double*>                                                 // for c
     Program;


typedef CGAL::Quadratic_program_solution<ET> Solution;


template<class T, class U> struct ComparePair1stDsc 
{
  bool operator()( const std::pair<T, const std::pair<U,U> >& a, 
                   const std::pair<T, const std::pair<U,U> >& b ) const 
  {
    return a.first > b.first; 
  }
};


CSolver::CSolver()
{
}


CSolver::~CSolver()
{
}


bool CSolver::Init (
        const clipper::Atom_list& Atm_lst, 
        const clipper::Xmap<float>& xm, 
        const clipper::HKL_info& hkli, 
        const double mu, 
        const double sd
        )
{
  hkl_info = hkli;
  exmap = xm;
  RSLN = hkl_info.resolution ( ).limit();
  
  // Unused variable
  BSMEAR = 5.6595 * pow(RSLN, 2.315);
  RADIUS = RSLN < 3.0 ? 0.7 + ( RSLN - 0.6 )/3.0 : 0.5 * RSLN;
  MU = mu;
  SD = sd;
  
  std::cout << "Mean, SD " << MU << " " << SD << std::endl;
  
  for (clipper::Xmap<float>::Map_reference_index ix = exmap.first(); 
          !ix.last(); ix.next() )
    if ( ( exmap[ix] - MU )/SD < -2 ) 
      exmap[ix] -= 10;
  
  mskmap.init(exmap.spacegroup(), exmap.cell(), exmap.grid_sampling());
  allAtm_list = clipper::Atom_list(Atm_lst);
  
  return true;
}


bool CSolver::CalcAtomGradient(
        clipper::Atom_list& atms, 
        std::vector<clipper::ftype>& grad, 
        const double occ 
        )
{
  std::vector<clipper::AtomShapeFn::TYPE> params;
  params.push_back( clipper::AtomShapeFn::X );
  params.push_back( clipper::AtomShapeFn::Y );
  params.push_back( clipper::AtomShapeFn::Z );  
  
  clipper::ftype rho_o, rho_c;
  std::vector<clipper::ftype> rho_c_grad ( params.size() );
  clipper::AtomShapeFn myAtomShapeFn(
          atms[0].coord_orth(), atms[0].element(), atms[0].u_iso(),
          atms[0].occupancy());
  myAtomShapeFn.agarwal_params() = params;
  
  for ( int i=0; i<rho_c_grad.size(); ++i )
    grad[i] = 0.0;
  
  const clipper::Grid_range gr(exmap.cell(), exmap.grid_sampling(), RADIUS);

  clipper::Xmap<clipper::ftype32>::Map_reference_coord i0, iu, iv, iw;

  for ( int i=0; i<atms.size(); i++ ) {
    const clipper::Coord_grid g0 = 
        exmap.coord_map( atms[i].coord_orth() ).coord_grid() + gr.min();
    const clipper::Coord_grid g1 = 
        exmap.coord_map( atms[i].coord_orth() ).coord_grid() + gr.max();
                
    i0 = clipper::Xmap_base::Map_reference_coord ( exmap, g0 );
    for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
      for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
        for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
          myAtomShapeFn.rho_grad ( iw.coord_orth ( ), rho_c, rho_c_grad );
          for ( int j=0; j<rho_c_grad.size(); ++j )
            grad[j] += 2*(exmap[iw] - occ*rho_c)*occ*rho_c_grad[j];
        }
  }
  
  return true;
}


double CSolver::CalcCCforRanking ( const clipper::Atom_list& atoms )
{
  double x, xx, y, yy, xy, cc, n;
  std::vector<float> mapO_msk, mapC_msk;
  clipper::Xmap_base::Map_reference_index ix;

  if ( ATMS_SIZE <= ATMS_MC_OFFSET || atoms.size() < ATMS_SIZE )
  {
  	std::cout << "ERROR: ATMS_SIZE not set correctly." << std::endl;
  	exit ( 10 );
  }

  clipper::Xmap<float> msk ( exmap.spacegroup(), exmap.cell(), exmap.grid_sampling());
  msk = 0.0;
  
  clipper::EDcalc_mask<float> edmask ( RADIUS );
//  clipper::Xmap<float> xmap_( exmap.spacegroup(), exmap.cell(), exmap.grid_sampling());
  
  edmask(msk, atoms );
  
  clipper::Xmap<float> mapC = CalculateMap( atoms );
  
  for ( ix = exmap.first(); !ix.last(); ix.next() )
    if ( msk[ix] > 0. )
    {
      mapO_msk.push_back ( exmap[ix] );
      mapC_msk.push_back ( mapC[ix] );
    }

  /* Calculate CC */
  x = xx = y = yy = xy = 0.0;
  
  for ( int l=0; l<mapO_msk.size(); ++l )
  {
    x += mapO_msk[l];
    xx += mapO_msk[l] * mapO_msk[l];
    yy += mapC_msk[l] * mapC_msk[l];
    y += mapC_msk[l];
    xy += mapO_msk[l] * mapC_msk[l];
  }    
  
  n = mapO_msk.size();
   
  cc = (xy/n - x*y/(n*n)) /
      (sqrt(xx/n-(x/n)*(x/n))*sqrt(yy/n-(y/n)*(y/n)));
    
  cc *= cc;
  
  std::cout << "CC = " << cc << std::endl;
  return cc;
}


double CSolver::CalcR2 ( const std::vector<double>& mapO, const clipper::Atom_list& atoms, const double coeff )
{
  double ss_err, ss_tot, rr, m;
  std::vector<float> mapO_msk, mapC_msk;
  clipper::Xmap_base::Map_reference_index ix;
  
  if ( ATMS_SIZE <= ATMS_MC_OFFSET || atoms.size() < ATMS_SIZE )
  {
  	std::cout << "ERROR: ATMS_SIZE not set correctly." << std::endl;
  	exit ( 10 );
  }

  clipper::Xmap<float> msk ( exmap.spacegroup(), exmap.cell(), exmap.grid_sampling());
  msk = 0.0;
  
  clipper::EDcalc_mask<float> edmask ( RADIUS );
  clipper::Xmap<float> xmap_( exmap.spacegroup(), exmap.cell(), exmap.grid_sampling());
  
  edmask (msk, atoms );
  
  clipper::Xmap<float> mapC = CalculateMap( atoms );
  
  for ( ix = exmap.first(); !ix.last(); ix.next() )
    if ( msk[ix] > 0. )
    {
      mapO_msk.push_back ( exmap[ix] );
      mapC_msk.push_back ( mapC[ix] );
    }
  
  /*Calculate ss_tot*/
  
  m = std::accumulate ( mapO_msk.begin(), mapO_msk.end(), 0.0)/mapO_msk.size();
  ss_tot = 0.0;
    
  for ( int l=0; l<mapO_msk.size(); ++l )
    ss_tot += (mapO_msk[l] - m)*(mapO_msk[l] - m);

  /*Subtract model */
  for ( int l=0; l<mapO_msk.size(); ++l )
    mapO_msk[l] -= coeff*mapC_msk[l];
      
  /*Calculate ss_err (residual sum of squares) */
  
  ss_err = 0.0;
    
  for ( int l=0; l<mapO_msk.size(); ++l )
    ss_err += (mapO_msk[l]*mapO_msk[l]);      

  rr  = 1-ss_err/ss_tot;


  return rr;
}  

float CSolver::CalcCC ( clipper::Atom_list& atms )
{
  clipper::Xmap<float> xmap(exmap.spacegroup(),exmap.cell(), exmap.grid_sampling());
  clipper::Xmap<unsigned char> bexmap(exmap.spacegroup(),exmap.cell(), exmap.grid_sampling());
  xmap = 0.0;
  bexmap = 0;
  
  float cc, x, xx, y, yy, xy, n;
  const clipper::Grid_range gr( exmap.cell(), exmap.grid_sampling ( ), RADIUS );

  clipper::Xmap<clipper::ftype32>::Map_reference_coord i0, iu, iv, iw;
/*
  for ( int i=0; i<atms.size(); i++ )
  {
  	clipper::AtomShapeFn sf ( atms[i].coord_orth(), atms[i].element ( ), atms[i].u_iso(), atms[i].occupancy() );
  	
    const clipper::Coord_grid g0 = exmap.coord_map( atms[i].coord_orth() ).coord_grid() + gr.min();
    const clipper::Coord_grid g1 = exmap.coord_map( atms[i].coord_orth() ).coord_grid() + gr.max();
                
    i0 = clipper::Xmap_base::Map_reference_coord ( exmap, g0 );
    for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
      for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
        for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
          xmap[iw] += sf.rho( iw.coord_orth() );
  }
  */
   clipper::Atom_list al ( allAtm_list );
   al.insert ( al.end(), atms.begin(), atms.end() );
    
   //clipper::HKL_info hkl_info ( exmap.spacegroup(), exmap.cell(), clipper::Resolution ( RSLN ), true );
   clipper::HKL_data<clipper::data32::F_phi> f_phi_calc ( hkl_info );
   clipper::SFcalc_iso_fft<float> (f_phi_calc, al );
   xmap.fft_from( f_phi_calc, clipper::Xmap_base::Normal );
  
  
/*  clipper::CCP4MAPfile mapin,mapout;
  mapout.open_write("calc_cc.map");
  mapout.export_xmap(xmap);
  mapout.close_write();
  exit(1); */
  
  cc = x = xx =y = yy = xy = n = 0.0; 
  
  for ( int i=0; i<atms.size(); ++i )
    if ( !atms[i].is_null() ) 
    {
      const clipper::Coord_grid g0 = exmap.coord_map( atms[i].coord_orth() ).coord_grid() + gr.min();
      const clipper::Coord_grid g1 = exmap.coord_map( atms[i].coord_orth() ).coord_grid() + gr.max();
      
      i0 = clipper::Xmap<float>::Map_reference_coord( exmap, g0 );      // sum all map contributions from this atoms
      for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
        for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
          for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
           if ( bexmap[iw] == 0 )
           {
             x += exmap[iw];
             xx += exmap[iw]*exmap[iw];
             yy += xmap[iw]*xmap[iw];
             y += xmap[iw];
             xy += exmap[iw]*xmap[iw];
             n++;   
           	 bexmap[iw] = 1;
           }
    }
  
  if ( atms.size() > 0 )
    cc = (xy/n-x*y/(n*n))/(sqrt(xx/n-(x/n)*(x/n))*sqrt(yy/n-(y/n)*(y/n)));
    
  cc *= cc;
  
  return cc;
}


void CSolver::CreateMaskMap(
        const std::vector<clipper::Atom_list>& vAtms
        )
{
  /* Create envelope for calculated configurations. 
   * Exclude those that do not scatter to at least xcl_threshold 
   * at all atomic centers in observed density.
   */

  if ( ATMS_SIZE <= ATMS_MC_OFFSET || vAtms[0].size() < ATMS_SIZE ) {
      std::cout << "ERROR: ATMS_SIZE not set correctly." << std::endl;
      exit ( 10 );
  }

  clipper::EDcalc_mask<float> edmask(RADIUS);
  clipper::Xmap<float> xmap_(
          exmap.spacegroup(), exmap.cell(), exmap.grid_sampling()
          );
  mskmap = 0.0;
  
  for ( int i=0; i<vAtms.size(); ++i ) {
      clipper::Atom_list atm_lst;
      for (int j=ATMS_MC_OFFSET; j<ATMS_SIZE; ++j)
          if (vAtms[i][j].occupancy() > 0)
              atm_lst.push_back(vAtms[i][j]);
      edmask(xmap_, atm_lst);
      mskmap += xmap_;
  }
}


clipper::Xmap<float> CSolver::CalculateMap(const clipper::Atom_list& atoms)
{
   clipper::Xmap<clipper::ftype32> 
       xmap(exmap.spacegroup(),exmap.cell(), exmap.grid_sampling());

   clipper::Atom_list al ( allAtm_list );
   clipper::Atom_list::const_iterator begin = atoms.begin() + ATMS_MC_OFFSET;
   clipper::Atom_list::const_iterator end = atoms.begin()+ ATMS_SIZE;
   al.insert ( al.end(), begin, end );
    
   clipper::HKL_data<clipper::data32::F_phi> f_phi_calc ( hkl_info );
   clipper::SFcalc_iso_fft<clipper::ftype32> (f_phi_calc, al );

   xmap.fft_from(f_phi_calc, clipper::Xmap_base::Normal);

   return xmap;
}


void CSolver::GoodnessOfFit(
        const std::vector<clipper::Atom_list>& vAtms, 
        std::vector<double>& sv 
        )
{
  double rr, ss_tot, ss_err, siga, m;
  int K =0;
  clipper::Xmap_base::Map_reference_index ix;
  
  for ( int i=0; i<sv.size(); ++i )
    if ( sv[i] > 0 )
      K++;
  
  double coeff[K], sigB[K];
  std::vector<double> dAIC (K), AIC (K), R2 ( K );
  std::vector<double> B[K], maps[K], map_o, AIC_map;
  std::vector<std::pair<double, std::pair <int,int> > > coeff_idx;
  std::vector<std::pair<double, std::pair <int,int> > > cc_idx;
  
  for ( ix = exmap.first(); !ix.last(); ix.next() )
    if ( mskmap[ix] > edmin )
       map_o.push_back(exmap[ix]);

  AIC_map = map_o;
  
  int map_idx = 0;
  for ( int i=0; i<sv.size(); ++i )
    if ( sv[i] > 0 )
    {
      clipper::Xmap<float> map_l = CalculateMap( vAtms[i] );
      for ( ix = map_l.first(); !ix.last(); ix.next() )
        if ( mskmap[ix] > edmin )
          maps[map_idx].push_back(map_l[ix]);
      coeff[map_idx] = sv[i];
      std::pair<int,int> idx_pair = std::make_pair (map_idx, i );
      coeff_idx.push_back ( std::make_pair ( coeff[map_idx], idx_pair ) );
      cc_idx.push_back ( std::make_pair ( CalcCCforRanking ( vAtms[i] ), idx_pair ) );
      map_idx++;
    }

  std::stable_sort( coeff_idx.begin(), coeff_idx.end(), ComparePair1stDsc<double,int>() );
  std::stable_sort( cc_idx.begin(), cc_idx.end(), ComparePair1stDsc<double,int>() );
  
  std::cout << "CHECK--Sorted Coefficient list:\n";
  for ( int i=0; i<K; ++i )
    std::cout << coeff_idx[i].first << " " << (coeff_idx[i].second).first << " " << (coeff_idx[i].second).second << "\n";
  std::cout << "END CHECK--Sorted list\n"; 
  
  std::cout << "CHECK--Sorted CC list:\n";
  for ( int i=0; i<K; ++i )
    std::cout << cc_idx[i].first << " " << (cc_idx[i].second).first << " " << (cc_idx[i].second).second << "\n";
  std::cout << "END CHECK--Sorted list\n"; 

  for ( int i=0; i<K; ++i ) {
    double RSS = 0.0;
    for ( int l=0; l<AIC_map.size(); ++l ) {
      AIC_map[l] -= sv[(cc_idx[i].second).second] *
                    maps[(cc_idx[i].second).first][l];
      RSS += AIC_map[l] * AIC_map[l];
    }
    AIC[i] = 2*(i+2) + AIC_map.size()* ( log (2*M_PI*RSS/AIC_map.size() ) + 1 );
    R2[i] = CalcR2 ( map_o, vAtms[(cc_idx[i].second).second], sv[(cc_idx[i].second).second] );
    double bic = AIC_map.size()* ( log (RSS/AIC_map.size() ) + (i+2)*log(AIC_map.size()) );
    std::cout << "AIC " << i << " = " << AIC[i] << " ";
    std::cout << "R^2 " << i << " = " << R2[i] << " ";
    std::cout << "BIC " << i << " = " << bic << std::endl;
  }
  
  std::vector<double>::const_iterator max_iter = AIC.begin();
  std::vector<double>::const_iterator min_iter = std::min_element ( AIC.begin(), AIC.end() );
  
  double maxAIC = * max_iter;
  double minAIC = * min_iter;
  const double totAIC = maxAIC - minAIC;

  double stop_crit = RSLN < 2.0 ? 0.0379*exp(1.2696*RSLN)-1 : -0.50;
  //stop_crit = -0.75;
  stop_crit = -0.95;
  std::cout << "stop_crit: " << stop_crit << std::endl;
  
  int stop;
  
  if ( max_iter == min_iter )
    stop = 0;
  else
  { 
  	stop = K-1;

    dAIC[0] = 0;  
    for ( int i=1; i<K; ++i )
    {
      dAIC[i] = ( AIC[i] - maxAIC )/ totAIC;
      std::cout << "dAIC " << i << " = " << dAIC[i] << std::endl;
    } 
  
    for ( int i=1; i<K; ++i )
      if ( dAIC[i] < stop_crit )
      {
        if ( std::abs ( stop_crit - dAIC[i-1] ) < std::abs ( dAIC[i] - stop_crit ) )
          stop = i-1;
        else
          stop = i;
        break;
      }
  }
  
  for ( int i=stop+1; i<K; ++i )
    sv[(cc_idx[i].second).second] = 0.0; 
  
  for ( int i=0; i<=stop; ++i )
  { 
    B[i] = maps[(cc_idx[i].second).first];
    for ( int j=0; j<=stop; ++j )
    {
      if ( i==j ) continue;
      for ( int l=0; l<maps[(cc_idx[i].second).first].size(); ++l )
  	    B[i][l] -=  sv[(cc_idx[i].second).second]*maps[(coeff_idx[i].second).first][l];
    }
  }
  
  /*Calculate ss_tot*/
  
  m = std::accumulate ( map_o.begin(), map_o.end(), 0.0)/map_o.size();
  ss_tot = 0.0;
    
  for ( int l=0; l<map_o.size(); ++l )
    ss_tot += (map_o[l] - m)*(map_o[l] - m);

  /*Subtract model */
  for ( int i=0; i<=stop; ++i )
    for ( int l=0; l<map_o.size(); ++l )
      map_o[l] -= sv[(cc_idx[i].second).second]*maps[(coeff_idx[i].second).first][l];
      
  /*Calculate ss_err (residual sum of squares) */
  
  ss_err = 0.0;
    
  for ( int l=0; l<map_o.size(); ++l )
    ss_err += (map_o[l]*map_o[l]);      

  std::cout << "R^2: " << 1-ss_err/ss_tot << std::endl;
  
  /*Calculate sd of B[i]*/
  
  for ( int i=0; i<=stop; ++i )
  {
    m = std::accumulate ( B[i].begin(), B[i].end(), 0.0)/B[i].size();
    double var = 0.0;
    
    for ( int l=0; l<B[i].size(); ++l )
      var += (B[i][l] - m)*(B[i][l] - m);
 
    var /= B[i].size();
    sigB[i] = sqrt(var);
  }

  /*Calculate sd of a*/
  
//  m = std::accumulate ( map_o.begin(), map_o.end(), 0.0)/map_o.size();
//  siga = 0;
//  for ( int l=0; l<map_o.size(); ++l )
//    siga += (map_o[l]-m)*(map_o[l]-m);
//  
  siga = sqrt(ss_err/map_o.size());
  
  std::cout << "sigs: ";
  for ( int i=0; i<=stop; ++i )
  {
    sigB[i] = siga/(sigB[i]*sqrt(map_o.size()-stop-1)); /*changed from map_o.size()*/
    std::cout << sigB[i] << " ";
  }
  std::cout << std::endl;
} 	


std::vector<double> CSolver::SetupLPandSolve(
        const std::vector<clipper::Atom_list>& vAtms, 
        const std::vector<int>& chosen, 
        const QP_t qp_t )
{
  std::vector<double> lower;
  const int coeff_size = chosen.size() + 1;
  std::vector<double> alpha_coeff[coeff_size];
  clipper::Xmap_base::Map_reference_index ix;
  
  int len = 0;
  for (ix = mskmap.first(); !ix.last(); ix.next())
    if (mskmap[ix] > edmin) len++;

  for (int i = 0; i < chosen.size(); ++i)
    alpha_coeff[i].reserve(len);

  for (int i = 0; i < chosen.size(); ++i) {
    clipper::Xmap<float> map_l = CalculateMap( vAtms[i] );

    for ( ix = map_l.first(); !ix.last(); ix.next()) {
      if ( mskmap[ix] > edmin )
         alpha_coeff[i].push_back(map_l[ix]);
    }
  }

  for ( ix = exmap.first(); !ix.last(); ix.next() )
    if ( mskmap[ix] > edmin )
       lower.push_back(-exmap[ix]);

  double sol[coeff_size];
  /*offset*/
  alpha_coeff[chosen.size()] = 
      std::vector<double> ( alpha_coeff[0].size(), 1.0 );

  double **rows_of_D = new double *[coeff_size];
  for ( int i=0; i<coeff_size; ++i )
    rows_of_D[i] = new double [coeff_size];

  double *ATbb = new double [coeff_size];

  for ( int i=0; i<coeff_size; ++i )
    for ( int j=0; j<coeff_size; ++j )
      rows_of_D[i][j] = 
          inner_product(alpha_coeff[i].begin(), 
                        alpha_coeff[i].end(), 
                        alpha_coeff[j].begin(), 0.0 );

  for ( int i=0; i<coeff_size; ++i )
    ATbb[i] = inner_product(
            alpha_coeff[i].begin(), alpha_coeff[i].end(), 
            lower.begin(), 0.0 );

  if (qp_t == QP_MIQP )
  {
    my_QP ( ATbb, rows_of_D, coeff_size, sol );
    my_MIQP ( ATbb, rows_of_D, coeff_size, sol );
  }
  else if ( qp_t == QP )
    my_QP ( ATbb, rows_of_D, coeff_size, sol );
  else if ( qp_t == MIQP )
  {
    for ( int i=0; i<coeff_size; ++i ) sol[i] =  0.01;
    my_MIQP ( ATbb, rows_of_D, coeff_size, sol );
  }

  /*   For QuadLP */
  std::vector<double> solvec;
  for ( int i=0;i<chosen.size();i++ )
    solvec.push_back(sol[i]);

  for ( int i=0; i<coeff_size; ++i ) {
    delete [] rows_of_D[i];
  }

  delete [] ATbb;
  delete [] rows_of_D;

  return solvec;
}


void CSolver::my_QP(
        double* ATbb_, double** rows_of_D_, 
        const int n, double s[] )
{
  double **cols_of_A = new double *[n];
  for ( int i=0; i<n; ++i )
    cols_of_A[i] = new double [1];
    
  bool *blower =  new bool [n];
  bool *bupper = new bool [n];
  double *lower = new double [n];
  double *upper = new double [n];
    
  for ( int i=0; i<n; i++ ) {
    blower[i] = true;
    bupper[i] = true;
    lower[i] = 0.0;
    upper[i] = 1.0;
  }
  lower[n-1] = -1.0;
    
  for ( int i=0; i<n; ++i )
    cols_of_A[i][0] = 1.0;
   
  cols_of_A[n-1][0] = 0.0;
    
  double b[]       = {1.0};

  CGAL::Const_oneset_iterator<CGAL::Comparison_result> rt( CGAL::SMALLER); 
    
  Program qp(n, 1, cols_of_A, b, rt, blower, lower, bupper, upper, rows_of_D_, ATbb_);

    // solve the program, using ET as the exact type
  CGAL::Timer tmr;
  tmr.reset(); tmr.start();
  Solution solver = CGAL::solve_quadratic_program(qp,ET());
  tmr.stop();
  std::cout << "Time (s) spent in QP = " << tmr.time() << std::endl;

  // we know that, don't we?
  if (solver.is_optimal()) { 
    std::cout << "Basic variables: ";
    for (Solution::Index_iterator it = solver.basic_variable_indices_begin();
	 it != solver.basic_variable_indices_end(); ++it)
      std::cout << *it << "  ";
    std::cout << std::endl;
  	
    std::cout << "Optimal feasible solution x: ";
    int k =0;
    for (Solution::Variable_value_iterator it = solver.variable_values_begin(); 
         it != solver.variable_values_end(); ++it) {
      s[k] = CGAL::to_double(*it);
      k++;
      std::cout << CGAL::to_double(*it) << " "; // variable values
    }
    std::cout << "f(x): " << CGAL::to_double(solver.objective_value()) << std::endl;
  } else {
    std::cout << "SOLUTION NOT OPTIMAL. STOP.";
    exit(1);
  }

  for ( int i=0; i<n; ++i ) {
    delete [] cols_of_A[i];
  }

  delete [] blower;
  delete [] bupper;
  delete [] lower;
  delete [] upper;
  delete [] cols_of_A;

  return;
}


void CSolver::my_MIQP(
        double* ATbb_, double** rows_of_D_, const int n, double s[]
        )
{
  int Num = 1;
  bool Indicator[n];
  Indicator[n - 1] = true;
  for(int i = 0; i < n - 1; ++i) {
  	if(s[i] > 1e-5) {
  		Indicator[i] = true;
  		++Num;
  	} else 
  		Indicator[i] = false;
  }
  double** upD = new double* [Num];
  for(int i = 0; i < Num; ++i) {
  	upD[i] = new double [Num];
  }
  double* upB = new double [Num];
  double* ups = new double [2 * Num - 1];
  
  int row = 0;
  for(int i = 0; i < n; ++i) {
  	if(Indicator[i]) {
  		int col = 0;
  		for(int j = 0; j < n; ++j) {
  			if(Indicator[j]) {
  				upD[row][col] = rows_of_D_[i][j];
  				std::cout << upD[row][col] << " ";
				++col;
  			}
  		}
  		std::cout << "\n";
  		++row;
  	}
  }

  for(int i = 0; i < Num; ++i) {
	for(int j = i + 1; j < Num; ++j) {
		upD[i][j] = upD[j][i];
	}
  }
    
  row = 0;
  for(int i = 0; i < n; ++i) {
  	if(Indicator[i]) {
  		upB[row] = ATbb_[i];
  		std::cout << upB[row] << " ";
  		++row;
  	}
  }
  
  std::cout << "Num of Unknows: " << Num << std::endl;
  std::cout << "Threshold: " << MILPthrshld << std::endl;
  MyMIQP(upD, upB, Num, ups, this->MILPthrshld);

  row = 0;
  for(int i = 0; i < n; ++i) {
  	if(Indicator[i]) {
  		s[i] = ups[row];
  		++row;
  	}
  	else
  		s[i] = 0.;
  }
  
  for(int i = 0; i < Num; ++i) {
  	delete [] upD[i];
  }
  delete upD;
  delete upB;
  delete ups;

}


// Unused function
void CSolver::CGAL_QuadLP ( const std::vector<double>& r_obs, std::vector<double> r_calc[], const int n, double s[] )
{
  double **rows_of_D = new double *[n];
  for ( int i=0; i<n; ++i )
    rows_of_D[i] = new double [n];
  double **cols_of_A = new double *[n];
  for ( int i=0; i<n; ++i )
    cols_of_A[i] = new double [1];
  double *ATbb = new double [n];
    
  bool *blower =  new bool [n];
  bool *bupper = new bool [n];
  double *lower = new double [n];
  double *upper = new double [n];
    
  for ( int i=0; i<n; i++ )
  {
    blower[i] = true;
    bupper[i] = true;
    lower[i] = 0.0;
    upper[i] = 1.0;
  }
  lower[n-1] = -1.0;
    
  //  rows_of_D[0][0] = 1.0; rows_of_D[0][1] = rows_of_D[1][0] = rows_of_D[1][1] = 0.0;
  //  cols_of_A[0][0] = 0.5; cols_of_A[1][0] = 1.0;
  //  ATbb[0] = 0.0; ATbb[1] = 3.0;
    

//    ofstream mat_file("matrix.dat");
//    for ( int i=0; i<r_calc[0].size(); ++i )
//    {
//      for ( int j=0; j<n; ++j )
//        mat_file << r_calc[j][i] << " ";
//      mat_file << "\n";
//    }
//    mat_file.flush();
    /* Use that in QR decomposition of A, R is Cholesky factor of ATA.
     * Compute QR of A, and then compute RTR = ATA.
     * Ideally, we'd pass R to QP_solver 
    gsl_matrix * A2 = gsl_matrix_calloc ( r_obs.size(), n );
    gsl_vector * t = gsl_vector_calloc ( n );

    for ( int i=0; i<n; ++i )
      for ( int j=0; j<r_calc[i].size(); ++j )
         gsl_matrix_set ( A2, j, i, r_calc[i][j] );
      
    gsl_linalg_QR_decomp ( A2, t );
    
    for ( int i=0; i<n; ++i )
      for ( int j=i; j<n; ++j )
        ATA[i][j] = gsl_matrix_get ( A2, i, j );
   
    gsl_matrix_free (A2);
    gsl_vector_free (t);
    
    for ( int r=0; r<n; ++r )
      for ( int c=0; c<n; ++c )
      {
         rows_of_D[r][c] = 0.0;
         for ( int k=0; k<n; ++k ) 
     rows_of_D[r][c] += ATA[k][r]*ATA[k][c];
      }
    End QR*/

  for ( int i=0; i<n; ++i )
    for ( int j=0; j<n; ++j )
      rows_of_D[i][j] = inner_product ( r_calc[i].begin(), r_calc[i].end(), r_calc[j].begin(), 0.0 );

  for ( int i=0; i<n; ++i )
    ATbb[i] = inner_product ( r_calc[i].begin(), r_calc[i].end(), r_obs.begin(), 0.0 );
      
  for ( int i=0; i<n; ++i )
    cols_of_A[i][0] = 1.0;
   
  cols_of_A[n-1][0] = 0.0;
    
  double b[]       = {1.0};

  CGAL::Const_oneset_iterator<CGAL::Comparison_result> rt( CGAL::SMALLER); 
    
  Program qp (n, 1, cols_of_A, b, rt, blower, lower, bupper, upper, rows_of_D, ATbb);

    // solve the program, using ET as the exact type
  CGAL::Timer tmr;
  tmr.reset(); tmr.start();
  Solution solver = CGAL::solve_quadratic_program(qp,ET());
  tmr.stop();
  std::cout << "Time (s) spent in QP = " << tmr.time() << std::endl;
    
  if (solver.is_optimal()) { // we know that, don't we?
  	std::cout << "Basic variables: ";
    for (Solution::Index_iterator it = solver.basic_variable_indices_begin();
	 it != solver.basic_variable_indices_end(); ++it)
      std::cout << *it << "  ";
    std::cout << std::endl;
  	
    std::cout << "Optimal feasible solution x: ";
    int k =0;
    for (Solution::Variable_value_iterator it = solver.variable_values_begin(); 
         it != solver.variable_values_end(); ++it) {
      s[k] = CGAL::to_double(*it);
      k++;
      std::cout << CGAL::to_double(*it) << " "; // variable values
    }
    std::cout << "f(x): " << CGAL::to_double(solver.objective_value()) << std::endl;
  }
  else
  {
    std::cout << "SOLUTION NOT OPTIMAL. STOP.";
    exit(1);
  }
    
  for ( int i=0; i<n; ++i )
  {
	delete [] rows_of_D;
    delete [] cols_of_A[i];
  }

  delete [] ATbb;
  delete [] blower;
  delete [] bupper;
  delete [] lower;
  delete [] upper;
  delete [] rows_of_D;
  delete [] cols_of_A;  
}


// Unused function
std::vector<clipper::Xmap<clipper::ftype32>::Map_reference_coord> CSolver::CreateCarve ( const clipper::Xmap<clipper::ftype32>& xm, 
		const clipper::Atom_list& atmList ) const
{
   std::vector<clipper::Xmap<clipper::ftype32>::Map_reference_coord> ref_coords;
   clipper::Xmap<unsigned char> bexmap( xm.spacegroup(), xm.cell(), xm.grid_sampling() );
   bexmap = 0;
   ref_coords.reserve ( 5000 );
   
   const clipper::Grid_range gr( xm.cell(), xm.grid_sampling ( ), 1.3 );

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


/*
 * Check to see if density is sufficiently strong around residue.
 */
// Unused function
bool CSolver::Include ( const clipper::Atom_list& atoms ) const
{
  bool include = false;
  double mo;
  std::vector<clipper::Xmap<float>::Map_reference_coord> ref_coords = CreateCarve ( exmap, atoms );

  mo = 0.0;
  for ( int i=0; i<ref_coords.size (); i++ )
  {
     clipper::Xmap<clipper::ftype32>::Map_reference_coord iw = ref_coords[i];
     mo += exmap[iw];
  }
  
  mo /= ref_coords.size();
  
  std::cout << "mo = " << mo << " MU, SD = " << MU << " " << SD << " " << "Z = " << (mo-MU)/SD << std::endl;
  
  if ( (mo-MU)/SD > -0.15 )
    include = true;
  return include;
}
