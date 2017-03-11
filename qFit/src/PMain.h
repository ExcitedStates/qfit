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


#ifndef PMAIN_H
#define PMAIN_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "PBasic.h"
#include "PChain.h"
#include "PExtension.h"
#include "PTools.h"
#include "PMyatomsf.h"
#include "PIKAlgorithms.h"
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-mmdb.h>
#include <mmdb_manager.h>
#include "ClpSimplex.hpp"
#include "ClpInterior.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include "CoinBuild.hpp"
#include "ClpPresolve.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "CoinSort.hpp"
#include "PSampMethods.h"
#include "POptimize.h"
#include "PScaleXMap.h"
#ifndef __RAMPROB_H_
#include "ramprob.h"
#include <mkl.h>
#include <gsl/gsl_linalg.h>
#endif
#include "QuadProg++.hh"



#define Min(a,b) a<b?a:b

struct sort_rank{
	double prob;
	int chosen;
};
struct protdist{
	int prot1;
	int prot2;
	double dist;
};

struct LPSolution{
	vector<double> soln;
	double obj_val;
};

class PXRayCrystApp{
 
	public:
		PXRayCrystApp(string paramfile, string debugfile);
		~PXRayCrystApp();
		void Init();
		void TestMap();
		void ScaleEDM();
		void ResolveEDM();
		
		vector<PProtein*> ResolveAdaptive(Xmap<float> map_t, bool forw);
		void ResolveBB(Xmap<float> map_t, vector<PProtein*> &loops, int dof, double prob_clust, bool forw );
		vector<double> ResolveSC(Xmap<float> map_t, vector<PProtein*> &loops_f, int numres, double prob_clust, bool forw);
		void ResolveStitchedProteins(Xmap<float> map_t,vector<PProtein*> &loops);
		void ResolveSCsWithFixedBB(Xmap<float> map_t, vector<PProtein*> &loops);
		void ResolveCluster(Xmap<float> map_t, vector<PProtein*> &loops);
		void ResolveWithPerturbed(Xmap<float> map_t, vector<PProtein*> &loops);
	
		void StitchProteins(vector<PProtein*> lpsF, vector<PProtein*> lpsB, vector<PProtein*> &loops, int stitchpt);
		vector<sort_rank> ClusterProteins(vector<PProtein*> lps, vector<sort_rank> srank, int idof, int K);
		vector<sort_rank> ClusterProteins(vector<PProtein*> lps, vector<sort_rank> srank, int idof, double prob_thresh, double dist_thresh,bool forw, bool WithSC);
		double RmsdOpt(PProtein *lp1, PProtein *lp2, int dof, bool forw, bool WithSC);
		LPSolution SetupLPandSolve(vector<PProtein*>& lps, const vector<int>& chosen, int dof, bool forw, bool WithSC);
		void ComputeScores(Xmap<float> map_t,vector<PProtein*>& orig_lps, vector<PProtein*>& calc_lps);
		
	private:
		string g_Paramfile;
		int loopS, loopE, g_Num_actual_loops;
		vector<PProtein*> g_Veclp;
		vector<PProtein*> g_RandFwLps;
		vector<PProtein*> g_RandBckLps;
 		vector<string> g_Lpfile;
		string g_Input_mapfile;
		string g_Input_pdbfile;
		string g_Final_mapfile;
		string g_Dir;
		string g_Id;
		string g_Dbgfile;
		ofstream g_Dbg;
		double g_Mapres;
		PProtein* g_Fullprot;
		clipper::Spgr_descr g_Spgrd;
		clipper::Cell_descr g_Celldescriptor;
		clipper::Spacegroup g_Spgr;
		clipper::Cell g_Cell;
		bool wSC;
		PProteinResidue *pres_anchN, *pres_anchC;
		static const int KERNEL_SIZE = 9;
		double kernelU[KERNEL_SIZE], kernelV[KERNEL_SIZE], kernelW[KERNEL_SIZE];
		static const int VAR = 0.0;
		Coord_grid box_min, box_max;
		vector<double> probs_f;
		Atom_list allAtom_list;
		
		Xmap<float> g_Finalmap;
		Xmap<float> mskmap, Zmap;
		int g_Num_anchN, g_Num_anchC;
		vector<double> g_TFac;
		double fmin, fmax;
		vector<clipper::Coord_grid> grid_F;
		int REAL, FULL_MAP;
		int g_Numresolve;
		int BATCH_SIZE;
		double ATA[MATRIX_DIM][MATRIX_DIM], ATb[MATRIX_DIM], CE[MATRIX_DIM][MATRIX_DIM], ce0[MATRIX_DIM], CI[MATRIX_DIM][MATRIX_DIM], ci0[MATRIX_DIM], x[MATRIX_DIM];
		
		
        bool isMC ( const PAtom* patm ) const;
		void AtomListNonSCopt(PProtein *lp, Atom_list &myatom_list, double occ,int dof, bool forw);
		void AtomListNonSCoptLastOnly(PProtein *lp, Atom_list &myatom_list, double occ,int dof, bool forw) const;
		bool TruncateSideChainAtom ( const PProteinResidue*, const PAtom* ) const;
		void AtomListSCopt(PProtein *lp, Atom_list &myatom_list, double occ,int dof, bool forw);
		void AtomListSCoptLastOnly(PProtein *lp, Atom_list &myatom_list, double occ,int dof, bool forw) const;
		float DensitySum ( PProtein* lp, const int dof, const bool forw, const bool WithSC );
		float DensityCC ( PProtein* lp, const int dof, const bool forw, const bool WithSC );
		bool SortByDensitySum ( vector<PProtein*>& lps, const int dof, const bool forw, const bool WithSC );
		void CreateGridBox(Xmap<float> & map, Coord_orth & p1, Coord_orth & p2);
		clipper::Xmap<float> GenerateMap(string infile);
		Xmap<float> CalculateMap(PProtein *lp,double occ,int dof, bool forw, bool WithSC);
		void CalculateMap(PProtein *lp,double occ,int dof, bool forw, bool WithSCi, Xmap<float>& );
		Xmap<float> CalculateMap(const Atom_list& atoms);
		void CalculateMap(const Atom_list& atoms, Xmap<float>&);
		Xmap<float> CalculateMap(PProtein *lp);
		void CreateMaskMap ( vector<PProtein*>& lps, const int dof, const bool forw, const bool WithSC );
                bool Include ( PProtein* lp, const int dof, const bool forw, const bool WithSC );
		void QuadLP ( const vector<double>& r_obs, vector<double> r_calc[], const int, double[] );
		void CGAL_QuadLP ( const vector<double>& r_obs, vector<double> r_calc[], const int, double[] );
		void UpdateObservedMap ( vector<PProtein*>&, const int, const bool, const bool );
		bool BondAngle ( vector<PProtein*>& in_Lps, const BondDirection dir, const int dof );
		bool OmegaAngle ( vector<PProtein*>& in_Lps, const BondDirection dir, const int dof );

		void DiscretizeOneDOF(vector<PProtein*>& in_Lps, int dof, int max_outLps, bool forw);
		void DiscretizeOneDOF(vector<PProtein*>& in_Lps, int dof, int max_outLps, double del_rot, bool forw);
		void DiscretizeOneDOFAndSolve(vector<PProtein*>& in_Lps, int dof, int max_outLps, double del_rot, Xmap<float> map_t, int num_solve_one_iter, double prob_clust, bool forw);
		void DiscretizeTwoDOFs(vector<PProtein*>& in_Lps, int dof1, int dof2, int max_outLps, double del_rot, bool forw);
		void DiscretizeTwoDOFs(vector<PProtein*>& in_Lps, int dof1, int dof2, int K1, int K2, double del_rot, bool forw);
		void DiscretizeTwoDOFsAndSolve(vector<PProtein*>& in_Lps, int dof1, int dof2, int K1, int K2, double del_rot, Xmap<float> map_t, int num_solve_one_iter, double prob_clust, bool forw);
		void DiscretizeTFacSC(vector<PProtein*>& in_Lps, int numres, vector<double> Tval, bool forw);
		void DiscretizeTFacSC_(vector<PProtein*>& in_Lps, int numres, vector<double> Tval, bool forw);
		void DiscretizeTFac(vector<PProtein*>& in_Lps, int dof, vector<double>& Tval, bool forw);
		void DiscretizeTFac(vector<PProtein*>& in_Lps, vector<string> atm_n, vector<int> res_indices, vector<double>& Tval);
		void GenerateRotamers(vector<PProtein*> &in_Lps, int numres);
		void FillAllAtomList (PProtein *lp );
		void GenerateRotamersAndSolve(vector<PProtein*> &in_Lps, int numres, Xmap<float> map_t, int num_solve_one_iter, double prob_clust, bool forw);
		
		double TorsionAngle(Vector3 u,Vector3 v,Vector3 w);
		void RandomPerturb(PProtein *lp);
		void RandomPerturb(PProtein *lp, vector<vector<CDof> > Dofs);
		void RotamerApply(PResidue *res, const vector<Real> &rotamer);
		
		Real GetFirstPhi ( PProteinResidue *pres );
		Real GetLastPsi ( PProteinResidue *pres );
		
		void InitKernels(ftype var, const Xmap<ftype32>& );
		void ConvolveSep(Xmap<ftype32> & map, Xmap<ftype32> & xmap);

		vector<string> getAllLines(string paramfile);
};
#endif

