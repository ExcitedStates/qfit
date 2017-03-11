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


/*! \file atomsf.h
  Header file for atomic scattering factors
  \ingroup g_fftmapbspl
*/
//C Copyright (C) 2000-2003 Kevin Cowtan and University of York
//L This library is free software; you can redistribute it and/or modify
//L it under the terms the CCP4 licence agreement as `Part 0' software,
//L which is defined as version 2.1 of the GNU Lesser General Public
//L License (LGPL). See the file 'license.txt' in the top CCP4 directory,
//L or the file 'COPYING' in this directory or its parent directory. 
//L
//L This library is distributed in the hope that it will be useful, but
//L WITHOUT ANY WARRANTY; without even the implied warranty of
//L MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L Lesser General Public License for more details. 
//L
//L You should have received a copy of the CCP4 license and/or GNU Lesser
//L General Public License along with this library; if not, write to the
//L CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.


#ifndef PMYATOMSF_H
#define PMYATOMSF_H


#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

using namespace clipper;
//using namespace clipper::data32;

  //! Atomic shape function object
  /*! The atomic scattering factor object is instantiated for each
    atom in turn, giving the atom parameters: position, element,
    occupancy and the isotropic or anisotropic U-value. (See
    clipper::Util for conversion from B-factors.). The methods of the
    class may then be called to return the scattering in reciprocal
    space or density in real space using either isotropic or
    anistropic models as required.

    If the atom only has an isotropic U, the faster isotropic methods
    will be used where available.

    This implementation uses the coefficients from Waasmaier & Kirfel
    (1995), Acta Cryst. A51, 416-431. The source data can be found at:
    ftp://wrzx02.rz.uni-wuerzburg.de/pub/local/Crystallography/sfac.dat
  */

class MyAtomShapeFn
{
 public:
  enum TYPE { X, Y, Z, Uiso, Occ, U11, U22, U33, U12, U13, U23 };
  //! null constructor
  MyAtomShapeFn() {}
  //! constructor: from atom object
  MyAtomShapeFn( const Atom& atom );
  //! constructor: from coord, element, isotropic U, occupancy
  MyAtomShapeFn( const Coord_orth& xyz, const String& element, const ftype u_iso = 0.0, const ftype occ = 1.0 );
  MyAtomShapeFn( ftype var, const Coord_orth& xyz, const String& element, const ftype u_iso = 0.0, const ftype occ = 1.0 );

  //! constructor: from coord, element, anisotropic U, occupancy
  MyAtomShapeFn( const Coord_orth& xyz, const String& element, const U_aniso_orth& u_aniso, const ftype occ = 1.0 );
  //! initialiser:  from atom object
  void init( const Atom& atom );
  
  //! initialiser: from coord, element, isotropic U, occupancy
  void init( const Coord_orth& xyz, const String& element, const ftype u_iso = 0.0, const ftype occ = 1.0 );
  void init( ftype var, const Coord_orth& xyz, const String& element, const ftype u_iso = 0.0, const ftype occ = 1.0 );
  
  //! initialiser: from coord, element, anisotropic U, occupancy
  void init( const Coord_orth& xyz, const String& element, const U_aniso_orth& u_aniso, const ftype occ = 1.0 );
  
  //! return scattering factor as a function of reflection posn
  ftype f( const Coord_reci_orth& rfl ) const;
  //! return electron density as a function of coordinate
  ftype rho( const Coord_orth& xyz ) const;
  
  //! return Agarwal density gradients as a function of coordinate
  bool rho_grad( const Coord_orth& xyz, std::vector<ftype>& grad ) const;
  
  
  bool rho_grad(const Coord_orth& xyz, ftype& grad) const;
  
  
  //! return Agarwal density curvatures as a function of coordinate
  bool rho_curv( const Coord_orth& xyz, Matrix<ftype>& curv ) const;
  
  //! return Agarwal density curvatures as a function of coordinate
  static bool rho_curv( const AtomShapeFn& sf1, const AtomShapeFn& sf2, const Coord_orth& uvw, Matrix<ftype>& curv );
  
  //! return (isotropic) scattering factor as a function of resolution
  ftype f( const ftype& invresolsq ) const;
  //! return (isotropic) electron density as a function of radius
  ftype rho( const ftype& rsq ) const;
  
  //! define parameters for Agarwal gradient/curvature calcs
  std::vector<TYPE>& agarwal_params() { return params; }
 private:
  //! look up atom coeffs
  void init( const String& element, const ftype& u_iso );
  void init( ftype var, const String& element, const ftype& u_iso );
  // members
  Coord_orth coord_;
  U_aniso_orth u_aniso_;
  ftype u_iso_, occ_;
  ftype a[6],  b[6];                //!< atom coeffs
  ftype aw[6], bw[6];               //!< intermediate results (iso)
  std::vector<Mat33sym<> > uaninv;  //!< intermediate results (aniso)
  bool is_iso;
  std::vector<TYPE> params;
};

#endif
