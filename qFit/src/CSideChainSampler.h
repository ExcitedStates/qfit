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
#define CSIDECHAINSAMPLER_H_

#ifndef CSOLVER_H_
#include "CSolver.h"
#endif

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <numeric>

#include "PBasic.h"
#include "PChain.h"
#include <PExtension.h>
#include "PResources.h"
#include "PTools.h"
#include "PMyatomsf.h"
#include "PIKAlgorithms.h"

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-mmdb.h>
#include <mmdb/mmdb_manager.h>
#include "POptimize.h"
#ifndef __RAMPROB_H_
#include "ramprob.h"
#endif
#include "myPDBIO.h"
#include "PScaleXMap.h"

static const char* SCdhdrls[20][20] = {
	/* { Resname, #Didrhal angles, # atoms , atom Names */
	{"ALA", "1", "CB", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""},
	{"ARG", "18", "CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"},
	{"ASN", "8", "CB", "HB2", "HB3", "CG", "OD1", "ND2", "HD21", "HD22", "", "", "", "", "", "", "", "", "", ""},
	{"ASP", "6", "CB", "HB2", "HB3", "CG", "OD1", "OD2", "", "", "", "", "", "", "", "", "", "", "", ""},
	{"CYS", "2", "CB", "SG", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""},
	{"GLN", "11", "CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "NE2", "HE21", "HE22", "", "", "", "", "", "", ""},
	{"GLU", "9", "CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "OE2", "", "", "", "", "", "", "", "", ""},
	{"GLY", "0", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""},
	{"HIS", "10", "CB", "HB2", "HB3", "CG", "ND1", "CE1", "HE1", "NE2", "CD2", "HD2", "", "", "", "", "", "", "", ""},
	{"ILE", "13", "CB", "HB", "CG1", "HG12", "HG13", "CD1", "HD11", "HD12", "HD12", "CG2", "HG21", "HG22", "HG23", "", "", "", "", ""},
	{"LEU", "13", "CB", "HB2", "HB3", "CG", "HG", "CD1", "HD11", "HD12", "HD13", "CD2", "HD21", "HD22", "HD23", "", "", "", "", ""},
	{"LYS", "16", "CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3", "CE", "HE2", "HE3", "NZ", "HZ1", "HZ2", "HZ3", "", ""},
	{"MET", "11", "CB", "HB2", "HB3", "CG", "HG2", "HG3", "SE", "CE", "HE1", "HE2", "HE3", "", "", "", "", "", "", ""},
//	{"MSE", "11", "CB", "HB2", "HB3", "CG", "HG2", "HG3", "SE", "CE", "HE1", "HE2", "HE3", "", "", "", "", "", "", ""},
//	{"PHE", "14", "CB", "HB2", "HB3", "CG", "CD1", "HD1", "CE1", "HE1", "CZ", "HZ", "CE2", "HE2", "CD2", "HD2", "", "", "", ""},
    {"PHE", "14", "CB", "HB2", "HB3", "CG", "HZ", "CZ", "CD1", "HD1", "CE1", "HE1", "CE2", "HE2", "CD2", "HD2", "", "", "", ""},
	{"PRO", "3", "CB", "CG", "CD", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""},
	{"SER", "2", "CB", "OG", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""},
	{"THR", "8", "CB", "HB", "OG1", "HG1", "CG2", "HG21", "HG22", "HG23", "", "", "", "", "", "", "", "", "", ""},
	{"TRP", "18", "CB", "HB2", "HB3", "CG", "CD1", "HD1", "NE1", "HE1", "CE2", "CZ2", "HZ2", "CH2", "HH2", "CZ3", "HZ3", "CE3", "HE3", "CD2"},
//	{"TYR", "15", "CB", "HB2", "HB3", "CG", "CD1", "HD1", "CE1", "HE1", "CZ", "OH", "HH", "CE2", "HE2", "CD2", "HD2", "", "", ""},
	{"TYR", "15", "CB", "HB2", "HB3", "CG", "OH", "CZ", "CD1", "HD1", "CE1", "HE1", "HH", "CE2", "HE2", "CD2", "HD2", "", "", ""},
	{"VAL", "10", "CB", "HB", "CG1", "HG11", "HG12", "HG13", "CG2", "HG21", "HG22", "HG23", "", "", "", "", "", "", "", ""}
	//{"VAL", "3", "CB", "CG1", "CG2", "", "", "", "", "", "", "", ""}
};


static const double ReferencePeptide[3][3] = { {0.503, 24.271, 6.725} /*CA*/,
                                               {2.832,21.322,5.949} /*CA+1*/,
                                               {1.380,22.332,7.992} /*O*/ };
                                               /* N: 0.764 24.940 5.394 */


static const double FD0[4][4] = {
    {-0.30019477009773254, -0.9518356323242188, -0.062385961413383484, 24.323608279629187}, 
    {-0.9510586857795715, 0.2936420142650604, 0.09623823314905167, 16.892928474206315},
    {-0.0732838436961174, 0.0882229283452034, -0.9934013485908508, 11.063804149400173}, 
    {0.0, 0.0, 0.0, 1.0}
};


static const double FD1[4][4] = {{-0.05113440006971359, -0.6621357798576355,-0.7476372718811035, 21.63548926737983},
                                 {-0.8933226466178894 ,  0.3650185763835907,-0.26217570900917053, 17.632786008435822},
                                 { 0.4464974105358124 ,  0.6544750928878784,-0.6101658940315247, -5.37073280048649},
                                 { 0.0                ,  0.0,0.0,1.0}};


static const double FD2[4][4] = {{-0.022191500291228294,-0.8806278109550476,0.4732886850833893,18.69209070411491},
                                 {-0.6436168551445007,0.3748406171798706,0.6672719717025757,11.113267179170059},
                                 {-0.7650260925292969,-0.28980880975723267,-0.5751051902770996,17.900355430109244},
                                 {0.0,0.0,0.0,1.0}};


static const double FD3[4][4] = {{-0.28591281175613403,-0.942484974861145,-0.17313562333583832,24.673346294025002},
                                 {-0.9384804964065552,0.23888792097568512,0.2493729591369629,17.507788294344422},
                                 {-0.19367024302482605,0.23378334939479828,-0.9528049230575562,7.369434679233631},
                                 {0.0,0.0,0.0,1.0}};


class MoveFunctCalculator: public FunctFunctor {

  public:

  MoveFunctCalculator(
          vector<PProtein*> loop, 
          vector<vector<CDof> > Dofs, 
          vector<PAtom*> AtomToMove, 
          vector<Vector3> MovementDirn
          )
  {
    funct_value = 0;
    Tfunct_value = 0;
    lps = loop;
    DofsToUse = Dofs;
    Atom = AtomToMove;
    for(int i=0; i<lps.size(); i++) {
      ndofs.push_back(lps[i]->size()*2);
      Vector3 pos(MovementDirn[i].x+Atom[i]->getPos().x, 
                  MovementDirn[i].y+Atom[i]->getPos().y,
                  MovementDirn[i].z+Atom[i]->getPos().z
                  );
      FinalPos.push_back(pos);
      std::cout << "Goal Pos: " << pos << std::endl;
    }
  }


  double operator()(double p[])
  {
    int k = 1;
    Tfunct_value = 0;
    for(int i=0;i<lps.size();i++){
      IKSolutions solutions;
      ChainMove CMove;
      vector<ChainMove> AllMove;
      AllMove.clear();
      for(int j=0; j < DofsToUse[i].size();j++){
        CMove.blockType = (DofsToUse[i])[j].blockType;
        CMove.dir = (DofsToUse[i])[j].dir;
        CMove.DOF_index = (DofsToUse[i])[j].DOF_index;
        CMove.degrees = p[k]*rad2deg;
        AllMove.push_back(CMove);
        k++;
      }
      solutions.push_back(AllMove);
      lps[i]->MultiRotate(solutions[0]);
      double Delta[4];
      Delta[1] = -FinalPos[i].x + Atom[i]->getPos().x;
      Delta[2] = -FinalPos[i].y + Atom[i]->getPos().y;
      Delta[3] = -FinalPos[i].z + Atom[i]->getPos().z;
      Tfunct_value += (Delta[1]*Delta[1]+Delta[2]*Delta[2]+Delta[3]*Delta[3]);
      //std::cout << "Tfunct_value " << Tfunct_value << "\n";
      lps[i]->AntiMultiRotate(solutions[0]);
    }
    return Tfunct_value;
  }


  double operator()(double p[],int i) {
    IKSolutions solutions;
    ChainMove CMove;
    vector<ChainMove> AllMove;
    AllMove.clear();
    for(int j=0;j<DofsToUse[i].size();j++){
      CMove.blockType = (DofsToUse[i])[j].blockType;
      CMove.dir = (DofsToUse[i])[j].dir;
      CMove.DOF_index = (DofsToUse[i])[j].DOF_index;
      CMove.degrees = p[j+1]*rad2deg;
      AllMove.push_back(CMove);
    }
    solutions.push_back(AllMove);
    lps[i]->MultiRotate(solutions[0]);
    double Delta[4];
    Delta[1] = -FinalPos[i].x + Atom[i]->getPos().x;
    Delta[2] = -FinalPos[i].y + Atom[i]->getPos().y;
    Delta[3] = -FinalPos[i].z + Atom[i]->getPos().z;
    funct_value = (Delta[1]*Delta[1]+Delta[2]*Delta[2]+Delta[3]*Delta[3]);
    //std::cout << "funct_value " << funct_value << "\n";
    lps[i]->AntiMultiRotate(solutions[0]);
    return funct_value;
  }
    
  private:

  double funct_value, Tfunct_value;
  vector<PProtein*> lps;
  vector<Vector3> FinalPos;
  vector<PAtom*> Atom;
  vector<int> ndofs;
  vector<vector<CDof> > DofsToUse;
};

class CSideChainSampler
{
public:
  CSideChainSampler();
  virtual ~CSideChainSampler();

  void Init(std::string inpdb, std::string inmtz, 
            const std::string, const int seqnum, 
            const bool fp = false
            );

  // Unused function
  void CalcConfs(
          PProtein* pChn, 
          const int start_res_idx, 
          const int stop_res_idx );

  void DoFindConfs(
          const std::string ChnID, 
          const int lSeqNum, 
          const int hSeqNum
          );

private:
  CMMDBManager myCMMDBManager;
  CSolver mCSolver;
  std::vector<std::string> atm_ids;
  int resnum;
  clipper::Metric_tensor Metric_tensor_direct;
  clipper::Mat33<> Metric_tensor_direct2;
  clipper::Matrix<> Ueigen;
  clipper::Vec3<> EVectorSum;
  Real CATempFactor;
  bool isMSE;
  bool flip_peptide;
  double MC_AMPL;
  double MC_SIG;
  int gSAMPLE_SIZE;

  enum FlipDirection {
      FlipBackward, 
      FlipForward
  };

  void CalcBndvctrRotPlaneIntsct(
          std::vector<double>& intsct, 
          std::vector<double>& bndvctr, 
          const std::vector<double>& bndtail, 
          const std::vector<double>& trmnlatm
          ) const;

  double CalcMinimizingAngle(
          std::vector<double>& bndvctr, 
          const std::vector<double>& bndtail, 
          const std::vector<double>& trmnlatms, 
          const std::vector<double>& anchrAtms
          ) const;

  void GetAtomIDs(PResidue *res);

  Atom_list GetAtoms(PResidue *res, const bool SConly);

  Atom_list GetAtoms(PProtein *pProt, const int excl_res);

  void SetAtoms(
          PCResidue res, 
          const Atom_list& atm_lst, 
          const char altLoc, 
          const double occupancy
          );

  void SetAtoms(
          PResidue *res, 
          const Atom_list& atm_lst, 
          const int c_idx = 1
          );
  
  void SamplePsiPhi(
          std::vector<PProtein*> pChnTbl, 
          const int res_idx
          );

  void RotamerApply(
          PResidue *res, 
          const vector<Real> &rotamer
          );

  void SetTFactor(PResidue *res, const Real TempFactor );

  void SetTFactor(
          PResidue *res, 
          const Real TempFactor, 
          const int chi_idx
          );

  void SampleRotamerNeighborhood(
          const PResidue *res, 
          const vector<Real> &rotamer, 
          vector<vector<Real> >& rotamer_nbrhd
          );

  void GenTrialPositions(
          PProteinResidue *pPres, 
          std::vector<clipper::Atom_list>& Atm_lst
          );

  bool BondAngle(
          PProteinResidue* pPres, 
          std::vector<clipper::Atom_list>& Atm_lst
          );

  bool SampleCB(
          PProteinResidue* pPres, 
          std::vector<clipper::Atom_list>& Atm_lst
          );

  bool Sample_OnMC(
          PProteinResidue* pPres, 
          const std::vector<double>& oc,  
          std::vector<clipper::Atom_list>& Atm_lst
          );

  bool SampleCB_MC(
          PProteinResidue* pPres, 
          std::vector<clipper::Atom_list>& Atm_lst
          );

  bool Calc_GradientMove(
          PResidue *res, 
          std::vector<clipper::ftype>& grad, 
          const double occ
          );

  float Calc_CC(PResidue *res);

  void SetADPAxis_Direct(const PCAtom pAtm);

  void SetChi(
          PResidue *res, 
          const int chi_idx, 
          const Real chi
          );

  bool SampleCB_MC_AtChi(
          PProteinResidue* pPres, 
          std::vector<clipper::Atom_list>& Atm_lst
          );

  void GenTrialPositions_AtChi(
          PProteinResidue* pPres, 
          const int chi_idx, 
          std::vector<clipper::Atom_list>& Atm_lst
          );

  bool Sample_NextChi(
          PProteinResidue* pPres, 
          const int next_chi_idx, 
          const std::vector<double>& oc, 
          std::vector<clipper::Atom_list>& Atm_lst
          );

  bool FlipPeptide(
          PProtein* pChn, 
          const FlipDirection, 
          const float dAngle
          );

  bool FlipTransformPeptide(
          PProtein* pChn,  
          const double (*FD)[4], 
          const int resno
          );

  void CoordFrame(
          const Vector3& a,  
          const Vector3& b, 
          const Vector3& c, 
          Matrix3& T, 
          Vector3& org
          );

  bool PepCoordFrame(
          PProtein* pChn, 
          Matrix3& T, 
          Vector3& org, 
          const int resno
          );

  void MakeHTMatrix(
          const Matrix3& rFrame, 
          const Vector3& org, 
          Matrix4& htMatrix
          );
  
  bool isMC(const PCAtom pAtm) const
  {
    return  (pAtm->CheckID("N", "*", "*") || 
             pAtm->CheckID("CA", "*", "*") || 
             pAtm->CheckID("C", "*", "*") || 
             pAtm->CheckID("O", "*", "*") 
             );
  };
 
  bool isLongSC ( const PResidue *res ) const
  {
  	return (res->getName() == PID::MET || 
                res->getName() == PID::LYS || 
                res->getName() == PID::ARG || 
                res->getName() == PID::GLN || 
                res->getName() == PID::GLU
                );
  };
  
  bool isRing ( const PResidue *res ) const
  {
  	return ( res->getName() == PID::TYR || res->getName() == PID::PHE ||
  	         res->getName() == PID::HIS || res->getName() == PID::TRP || 
		 res->getName() == PID::GLN || res->getName() == PID::GLU );
  }; 


  vector<IKSolutions> MoveAtom(
          vector<PProtein*> loops, 
          vector<vector<CDof> > Dofs, 
          vector<PAtom*> AtomToMove, 
          vector<Vector3> MovementDirn
          )
  {
    MoveFunctCalculator *FunctToOptimize = 
        new MoveFunctCalculator(loops, Dofs, AtomToMove, MovementDirn);
    POptimize *opt = new POptimize(loops, Dofs);
    vector<IKSolutions> CombinedSoln = opt->OptimizeNullSpace(FunctToOptimize).LoopSol;
    for(int i=0;i<loops.size();i++){
      IKSolutions Soln = CombinedSoln[i];
      if ((Soln[1])[0].DOF_index!=-1){
        loops[i]->MultiRotate(Soln[0]);
        loops[i]->MultiRotate(Soln[1]);
          }
      else cout<<"No solution found"<<endl;
    }
    delete FunctToOptimize;
    delete opt;
    return CombinedSoln;
  }
  
  
  IKSolutions MoveAtom(
          PProtein *loop, 
          PAtom* AtomToMove, 
          Vector3 MovementDirn)
  {
    vector<PProtein *> loops;
    loops.push_back(loop);
    vector<vector<CDof> > Dofs;
    Dofs = PTools::GetBBDofs(loops);
    vector<PAtom*> Atom;
    Atom.push_back(AtomToMove);
    vector<Vector3> MovementDir;
    MovementDir.push_back(MovementDirn);
    vector<IKSolutions> CombinedSoln = MoveAtom(loops, Dofs, Atom, MovementDir);
    return CombinedSoln[0];
  }

};

#endif /*CSIDECHAINSAMPLER_H_*/
