#ifndef _LABELING_
#define _LABELING_

#include <string>
#include <cfloat>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <numeric>
#include <map>
#include <set>
#include <mmdb/mmdb_manager.h>
#include <mmdb/mmdb_tables.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

struct Node {
    int ResInd;
    PPCAtom atoms;
    int atomNo;
    AltLoc altLoc;
};

struct PotentialPair {
    std::string a;
    std::string b;
    double epsilon;
};

bool EqualPair(PotentialPair pair1, PotentialPair pair2) {
    if((pair1.a == pair2.a && pair1.b == pair2.b) 
    || (pair1.a == pair2.b && pair1.b == pair2.a)) {
        return true;
    }
    return false;
}

class Label {
public:
    Label() {}
    bool ReadFile(std::string);
    bool SetParameter(int, int, std::string, bool);
    bool MCMCLabeling(int);
    void setupEpsilonVals(const std::string &);
    void setupVanDerWaalsRadii();
    double vanDerWaalsEnergy(const PCAtom, const PCAtom);
private:
    std::vector<double> Metropolis(std::vector<std::vector<double> >&, 
                                   std::vector<std::vector<int> >&, 
                                   std::vector<std::vector<int> >&, 
                                   std::vector<int>&); 
    double energy(Node, Node);  
    bool isMChain(PCAtom);
    double vanDerWaalsRadius(PCAtom);
    PCMMDBManager myPCMMDBManager;
    std::vector<PCResidue> pcresidues;
    std::string file;
    std::string outputfile;
    int MolelNo;
    int NumOfSim;
    int NumOfChains;
    bool useOldH;
    std::vector<PotentialPair> atomPairs;
    std::map<std::string, double> myVdWaalsRadius;
    std::map<std::string, double> myNewVdWaalsRadius;
    #define NEW_RADIUS_H    1.22 // neither aromatic nor polar
    #define NEW_RADIUS_HAR  1.05 // aromatic
    #define NEW_RADIUS_HPOL 1.05 // polar
    #define NEW_RADIUS_HAP  1.05 // aromatic and polar
};
#endif
