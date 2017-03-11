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


#include "PScaleXMap.h"

template<class T> bool CompareMapIndex ( const T& a, const T& b)
{
  return a.index() < b.index();
};

template<class T> bool EqualMapIndex ( const T& a, const T& b)
{
  return a.index() == b.index();
};

PScaleXMap::PScaleXMap ( )
{
}

PScaleXMap::~PScaleXMap()
{
}


std::vector<clipper::Xmap<clipper::ftype32>::Map_reference_coord> 
    PScaleXMap::CreateCarve(
            const clipper::Xmap<clipper::ftype32>& xm, 
            const clipper::Atom_list& atmList
            ) const
{
   std::vector<clipper::Xmap<clipper::ftype32>::Map_reference_coord> ref_coords;
   std::vector<int> index;
   clipper::Xmap<unsigned char> 
       bexmap(xm.spacegroup(), xm.cell(), xm.grid_sampling());
   bexmap = 0;
   ref_coords.reserve(5000);
   
   const clipper::Grid_range 
       gr(xm.cell(), xm.grid_sampling(), RADIUS);

   clipper::Xmap<clipper::ftype32>::Map_reference_coord i0, iu, iv, iw;

   for ( int i=0; i<atmList.size(); i++ ) {
      const clipper::Coord_grid 
          g0 = xm.coord_map(atmList[i].coord_orth()).coord_grid() + gr.min();
      const clipper::Coord_grid 
          g1 = xm.coord_map(atmList[i].coord_orth()).coord_grid() + gr.max();
                
      i0 = clipper::Xmap_base::Map_reference_coord ( xm, g0 );
      for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
         for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
            for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
                if ( bexmap[iw] == 0 ) {
                    ref_coords.push_back ( iw );
                    bexmap[iw] = 1;
                }
   }

   return ref_coords;
}


bool PScaleXMap::DoScaleXMap(
        clipper::Xmap<clipper::ftype32>& xm, 
        clipper::HKL_info& hkl_info, 
        const std::string inmtz, 
        const clipper::Atom_list& atmList, 
        clipper::ftype32& m, 
        clipper::ftype32& var
        )
{
   clipper::ftype32 S, k;
   std::vector<double> mo, mc;
   clipper::Xmap<clipper::ftype32> fxmap;
   bool success = false;
   
  clipper::CCP4MAPfile mapout;
  
  /* Calculate Observed map: */
  clipper::CCP4MTZfile mtzin;
  mtzin.open_read ( inmtz.c_str ( ) );
  mtzin.import_hkl_info ( hkl_info, false ); 
  clipper::HKL_data<clipper::data32::F_phi> f_phi ( hkl_info ); 
  mtzin.import_hkl_data( f_phi, "*/*/[FWT,PHWT]" );
  mtzin.close_read();
  
  xm.init(hkl_info.spacegroup(), hkl_info.cell(), 
          clipper::Grid_sampling(hkl_info.spacegroup(), 
              hkl_info.cell(), hkl_info.resolution()
              )
          );
  xm.fft_from( f_phi, clipper::Xmap_base::Normal );
  
  // Write map to file
  mapout.open_write( "observed.map");
  mapout.export_xmap(xm);
  mapout.close_write(); 
   
  /* Calculate Atom map */
  clipper::HKL_data<clipper::data32::F_phi> f_phi_calc ( hkl_info );
  clipper::SFcalc_iso_fft<float>(f_phi_calc, atmList);
  fxmap.init(hkl_info.spacegroup(), hkl_info.cell(), 
          clipper::Grid_sampling(hkl_info.spacegroup(),  
              hkl_info.cell(), hkl_info.resolution()
              )
          );
  fxmap.fft_from( f_phi_calc, clipper::Xmap_base::Normal );

  // Write map to file
  mapout.open_write( "calc.map");
  mapout.export_xmap(fxmap);
  mapout.close_write();

  const double rsln = hkl_info.resolution().limit();

  // Unused variable
  BSMEAR = 5.6595*pow ( rsln, 2.315 );
  RADIUS = rsln < 3.0 ? 0.7 + ( rsln - 0.6 )/3.0 : 0.5 * rsln;
 
  // Get all map coordinates of interest
  std::vector<clipper::Xmap<clipper::ftype32>::Map_reference_coord> 
      ref_coords = CreateCarve(xm, atmList);
    
  // Get all map values of interest
  for ( int i=0; i<ref_coords.size (); i++ ) {
     clipper::Xmap<clipper::ftype32>::Map_reference_coord iw = ref_coords[i];
     mo.push_back(xm[iw]);
     mc.push_back(fxmap[iw]);
  }
  
  // Get average and variance of masked observed map values
  Stats(mo, m, var);
  
  // Get optimal scaling factor and mean-difference. 
  if (L2Fit( mo, mc, S, k)) {
    clipper::Xmap_base::Map_reference_index ix;
     // Scale the observed map to the calculated map
    for ( ix = xm.first(); !ix.last(); ix.next() )
      xm[ix] = S * xm[ix] + k;
    // Set the mean and variance
    m =  S * m + k;
    var *= S * S;
    success = true;
  }
  
  return success;
}


bool PScaleXMap::Stats(
        const std::vector<double>& xmap_o, 
        clipper::ftype32& m, 
        clipper::ftype32& var
        )
{
  m = std::accumulate(xmap_o.begin(), xmap_o.end(), 0.0) / xmap_o.size();
  var = 0.0;
  
  for ( int i=0; i<xmap_o.size(); ++i )
    var += (xmap_o[i] - m) * (xmap_o[i] - m);
  
  var /= xmap_o.size();
  return true;
}


bool PScaleXMap::L2Fit(
        const std::vector<double>& xmap_o, 
        const std::vector<double>& xmap_c, 
        clipper::ftype32& S_, 
        clipper::ftype32& k_ 
        ) const
{
  double m1 = std::accumulate(xmap_o.begin(), xmap_o.end(), 0.0)/xmap_o.size();
  double m2 = std::accumulate(xmap_c.begin(), xmap_c.end(), 0.0)/xmap_c.size();
  
  // Variance and covariance
  double s1 = 0.0, s2 = 0.0;
  
  for ( int i=0; i<xmap_o.size(); ++i )
  {
    s1 += (xmap_o[i] - m1) * (xmap_o[i] - m1);
    s2 += (xmap_o[i] - m1) * (xmap_c[i] - m2);
  }

  S_ = s2 / s1;
  k_ = m2 - S_ * m1;
    
  std::cout << "L2 Scaling: S = " << S_ << " k = " << k_ << "\n";
  return true;
}


// Unused function
bool PScaleXMap::L1Fit(
        const std::vector<double>& xmap_o, 
        const std::vector<double>& xmap_c, 
        clipper::ftype32& S_, 
        clipper::ftype32& k_
        ) const
{
  if ( xmap_o.size() != xmap_c.size() ) return false;
  size_t it =0, alpha=0, j;
  const size_t maxit = 100;
  int p;
  int ip1 = -1;
  double xmap_o_alpha, xmap_c_alpha, sum, o, c;
  std::vector<double> T ( xmap_o.size() );
  std::vector<size_t> D ( xmap_o.size() );
  bool done = false;
  
  while ( it < maxit && !done )
  {
    sum = 0;
    xmap_o_alpha = xmap_o[alpha];
    xmap_c_alpha = xmap_c[alpha];
    for ( size_t i=0; i<xmap_o.size(); i++ )
    {
      o = xmap_o[i] - xmap_o_alpha;
      c = xmap_c[i] - xmap_c_alpha;
      T[i] = ( std::abs ( o ) > 1E-6 ) ? c/o : 1E19;
      D[i] = i;
      sum += std::abs ( o );
    }
    sum /= 2.0;

    sort_idxtbl ( T.begin(), T.end(), D.begin() );    

    o = 0;
    j=0;
    while ( o < sum )
    {
      p = D[j];
      j++;
      o += std::abs ( xmap_o[p]-xmap_o_alpha );
      if ( p == ip1 ) 
      {
        done = true;
        break;
      }
    }
    ip1 = alpha;
    alpha = p;
    it++;
  }
  
  S_ = T[p];
  k_ = xmap_c[alpha] - T[p]*xmap_o[alpha];
  
  std::cout << "L1 Scaling: S = " << S_ << " k = " << k_ << " iters = " << it << "\n";
  return done;
}
