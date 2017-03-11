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
#define CSOLVER_H_

#include <cstring>
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-mmdb.h>

enum QP_t {
    QP_MIQP, 
    QP, 
    MIQP
};


class CSolver
{
public:

  double RSLN;

  CSolver();

  virtual ~CSolver();

  bool Init(
          const clipper::Atom_list& Atm_lst, 
          const clipper::Xmap<float>& xm, 
          const clipper::HKL_info& hkli, 
          const double mu, 
          const double sd
          );

  void CreateMaskMap(
          const std::vector<clipper::Atom_list>& vAtms
          );

  std::vector<double> SetupLPandSolve(
          const std::vector<clipper::Atom_list>& vAtms, 
          const std::vector<int>& chosen, 
          const QP_t = QP_MIQP
          );

  bool CalcAtomGradient(
          clipper::Atom_list& atms, 
          std::vector<clipper::ftype>& grad, 
          const double occ
          );

  bool Include(
          const clipper::Atom_list& atoms
          ) const;

  float CalcCC(clipper::Atom_list& atms);

  void GoodnessOfFit(
          const std::vector<clipper::Atom_list>& vAtms, 
          std::vector<double>& sv
          );
  
  void SetATMS_SIZE(const size_t N)
  {
    ATMS_SIZE = N;
  }
  
  void SetATMS_MC_OFFSET(const size_t N = 4)
  {
    ATMS_MC_OFFSET = N;
  }

  
  void SetMILPThreshold ( const double t )
  {
    MILPthrshld = t;
  }

private:
  
  static const double edmin = 0.0;
  size_t ATMS_MC_OFFSET;
  size_t ATMS_SIZE;
  clipper::HKL_info hkl_info;
  clipper::Atom_list allAtm_list;
  double BSMEAR, RADIUS;
  clipper::Xmap<float> mskmap, exmap;
  double MU, SD;
  double MILPthrshld;
  
  clipper::Xmap<float> CalculateMap(const clipper::Atom_list& atoms);

  void my_QP(
          double* ATbb_,  
          double** rows_of_D_, 
          const int n, 
          double s[]
          );

  void my_MIQP(
          double* ATbb_,  
          double** rows_of_D_, 
          const int n, 
          double s[]
          );

  void  CGAL_QuadLP(
          const std::vector<double>& r_obs, 
          std::vector<double> r_calc[], 
          const int n, 
          double s[]);

  std::vector<clipper::Xmap<clipper::ftype32>::Map_reference_coord> 
      CreateCarve(
              const clipper::Xmap<clipper::ftype32>& xm, 
              const clipper::Atom_list& atmList
              ) const;

  double CalcR2(
          const std::vector<double>& mapO, 
          const clipper::Atom_list& atoms, 
          const double coeff
          );

  double CalcCCforRanking(const clipper::Atom_list& atoms);
};

#endif /*CSOLVER_H_*/
