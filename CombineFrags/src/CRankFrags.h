#ifndef CRANKFRAGS_H_
#define CRANKFRAGS_H_

#include <vector>
#include <list>
#include <fstream>
#include <iomanip>
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-mmdb.h>
//#include <math_.h>
#include <math.h>

class CRankFrags
{
private:
  
  clipper::Xmap<int> xskl;
  clipper::Xmap<float> xmap;
  clipper::HKL_info hkl_info;
  clipper::HKL_data<clipper::data32::F_sigF> f_sigf;
  clipper::ftype mapres, RADIUS, BSMEAR;
  int mainCorrHandle, sideCorrHandle;
  bool registerUD;
  void Normalize ( clipper::Xmap<clipper::ftype32>& xm, const clipper::Atom_list& atmList );
  std::vector<clipper::Xmap<clipper::ftype32>::Map_reference_coord> CreateCarve ( const clipper::Xmap<clipper::ftype32>& xm, 
		const clipper::Atom_list& atmList ) const;
  
public:
  CRankFrags();
  virtual ~CRankFrags();

  bool CalcMap (  const std::string& mtzFile );
  bool calcEDMCC ( const PCMMDBManager pCMMDBManager, const int nMdl );
  int Rankum ( const PCMMDBManager );
  bool InsertFragment ( PCChain, const PCChain, const double ccTHRESHOLD = -1 );
  bool Register ( PCMMDBManager pCMMDBManager );
  bool InsertBestTerminal ( PCMMDBManager pMMDBManager, PCMMDBManager pFragManager, const bool isForward );
  bool CreateSkeleton ( const clipper::Xmap<float>& xmap );
};

#endif /*CRANKFRAGS_H_*/
