#include "Labeling.hpp"
using namespace std;

double norm(double x1, double y1, double z1,
            double x2, double y2, double z2) {
    double a, b, c;
    a = x1 - x2;
    b = y1 - y2;
    c = z1 - z2;
    return sqrt(a * a + b * b + c * c);
}

void Label::setupEpsilonVals(const string &epsilonFile) {
    for(int i = 0; i < nElementNames; ++i) {
        myVdWaalsRadius.insert(make_pair(string(ElementName[i]), VdWaalsRadius[i]));
    }
    xmlDocPtr doc = xmlReadFile(epsilonFile.c_str(), NULL, 0);
    if (doc == NULL) {
        printf("Error: File not found.");
    }

    xmlNodePtr start_node = doc->children->children;
    for (xmlNodePtr cur_node = start_node; cur_node != NULL; cur_node = cur_node->next) {
        if (string(reinterpret_cast<const char *>(cur_node->name)) == string("atom_pair")) {
            string first_elem = "", second_elem = "", epsilon_val = "";

            for (xmlNodePtr cur_child_node = cur_node->children;
                     cur_child_node != NULL;
                     cur_child_node = cur_child_node->next) {

                if (string(reinterpret_cast<const char *>(cur_child_node->name)) == string("first")) {
                    first_elem = reinterpret_cast<const char *>(cur_child_node->children->content);
                } else if (string(reinterpret_cast<const char *>(cur_child_node->name)) == string("second")) {
                    second_elem = reinterpret_cast<const char *>(cur_child_node->children->content);
                } else if (string(reinterpret_cast<const char *>(cur_child_node->name)) == string("epsilon")) {
                    epsilon_val = reinterpret_cast<const char *>(cur_child_node->children->content);
                }   
            }

            //PResources::AddEpsilonValue(make_pair(first_elem, second_elem),
            //                     atof(epsilon_val.c_str()));
            PotentialPair tmpPair;
            tmpPair.a = string(" ") + first_elem;
            tmpPair.b = string(" ") + second_elem;
            tmpPair.epsilon = atof(epsilon_val.c_str());
            //printf("%d, %d, %lf\n", tmpPair.a.size(), tmpPair.b.size(), atof(epsilon_val.c_str()));
            bool flag = true;
            for(int i = 0; i < atomPairs.size(); ++i) {
                if(EqualPair(atomPairs[i], tmpPair)) {
                    flag = false;
                    break;
                }
            }
            if(flag) {
                atomPairs.push_back(tmpPair);
            }
        }
    }
    /*
    printf("%d\n", atomPairs.size());
    for(int i = 0; i < atomPairs.size(); ++i) {
        printf("%s, %s: %lf\n", atomPairs[i].a.c_str(), atomPairs[i].b.c_str(), atomPairs[i].epsilon);
    }
    */
    xmlFreeDoc(doc);
}


// Stores MMDB van der Waals radii for all atoms (heavy or hydrogen) in memory.
// If the user did not request old hydrogen van der Waals radii (based on nuclear 
// positions), also stores "new Reduce" van der Waals radii for all hydrogens.
// TODO: Make this read an XML file instead of hard-coding the numbers here!
void Label::setupVanDerWaalsRadii() {
    for(int i = 0; i < nElementNames; ++i) {
        myVdWaalsRadius.insert(make_pair(string(ElementName[i]), VdWaalsRadius[i]));
    }
    if(!useOldH) {
        myNewVdWaalsRadius[string("LYS, HZ3")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("LYS, HZ2")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("LYS, HZ1")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("LYS,3HZ ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("LYS,2HZ ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("LYS,1HZ ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("LYS, HE3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS, HE2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS, HD3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS, HD2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS, HG3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS, HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS,2HE ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS,1HE ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS,2HD ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS,1HD ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS,2HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS,1HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LYS,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLY, HA3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLY, HA2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLY,2HA ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLY,1HA ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLU, HG3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLU, HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLU, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLU, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLU,2HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLU,1HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLU,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLU,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THR,HG23")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THR,HG22")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THR,HG21")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THR,3HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THR,2HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THR,1HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THR, HG1")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("THR, HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ALA, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ALA, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ALA, HB1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ALA,3HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ALA,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ALA,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PHE, HZ ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("PHE, HE2")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("PHE, HE1")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("PHE, HD2")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("PHE, HD1")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("PHE, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PHE, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PHE,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PHE,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ARG,HH22")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ARG,HH21")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ARG,HH12")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ARG,HH11")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ARG,HH22")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ARG,HH21")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ARG,HH12")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ARG,HH11")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ARG,2HH2")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ARG,1HH2")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ARG,2HH1")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ARG,1HH1")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ARG, HE ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ARG, HD3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ARG, HD2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ARG, HG3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ARG, HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ARG, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ARG, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ARG,2HD ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ARG,1HD ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ARG,2HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ARG,1HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ARG,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ARG,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HIS, HE2")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string("HIS, HE1")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("HIS, HD2")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("HIS, HD1")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string("HIS, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HIS, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HIS,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HIS,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET, HE3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET, HE2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET, HE1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET, HG3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET, HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET,3HE ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET,2HE ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET,1HE ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET,2HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET,1HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MET,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ASP, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ASP, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ASP,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ASP,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("SER, HG ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("SER, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("SER, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("SER,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("SER,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ASN,HD22")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ASN,HD21")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ASN,HD22")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ASN,HD21")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ASN,2HD2")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ASN,1HD2")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ASN, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ASN, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ASN,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ASN,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("TYR, HH ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("TYR, HE2")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("TYR, HE1")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("TYR, HD2")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("TYR, HD1")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("TYR, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("TYR, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("TYR,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("TYR,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("CYS, HG ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("CYS, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("CYS, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("CYS,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("CYS,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLN,HE22")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("GLN,HE21")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("GLN,HE22")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("GLN,HE21")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("GLN,2HE2")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("GLN,1HE2")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("GLN, HG3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLN, HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLN, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLN, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLN,2HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLN,1HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLN,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLN,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,HD23")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,HD22")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,HD21")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,HD13")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,HD12")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,HD11")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,3HD2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,2HD2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,1HD2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,3HD1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,2HD1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,1HD1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU, HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("LEU,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PRO, HD3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PRO, HD2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PRO, HG3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PRO, HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PRO, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PRO, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PRO,2HD ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PRO,1HD ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PRO,2HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PRO,1HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PRO,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PRO,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("VAL,HG23")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("VAL,HG22")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("VAL,HG21")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("VAL,HG13")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("VAL,HG12")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("VAL,HG11")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("VAL,3HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("VAL,2HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("VAL,1HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("VAL,3HG1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("VAL,2HG1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("VAL,1HG1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("VAL, HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,HD13")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,HD12")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,HD11")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,HG23")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,HG22")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,HG21")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,HG13")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,HG12")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,3HD1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,2HD1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,1HD1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,3HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,2HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,1HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,2HG1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE,1HG1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ILE, HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("TRP, HH2")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("TRP, HZ3")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("TRP, HZ2")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("TRP, HE3")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("TRP, HE1")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string("TRP, HD1")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("TRP, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("TRP, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("TRP,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("TRP,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  U, H6 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("  U, H5 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("  U, H3 ")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string("  U, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  U, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  U, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  U, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("U  , H6 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("U  , H5 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("U  , H3 ")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string("U  , H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("U  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("U  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("U  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("URA, H6 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("URA, H5 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("URA, H3 ")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string("URA, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("URA, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("URA, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("URA, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H6 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("  T,3H5M")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T,2H5M")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T,1H5M")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H53")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H52")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H51")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T,3H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T,2H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T,1H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T,3H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T,2H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T,1H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H53")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H52")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H51")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H73")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H72")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H71")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H73")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H72")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H71")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H73")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H72")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H71")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H3 ")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string("  T, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  T, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H6 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("T  ,3H5M")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  ,2H5M")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  ,1H5M")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H53")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H52")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H51")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  ,3H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  ,2H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  ,1H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  ,3H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  ,2H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  ,1H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H53")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H52")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H51")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H73")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H72")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H71")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H73")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H72")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H71")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H73")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H72")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H71")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H3 ")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string("T  , H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("T  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H6 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("THY,3H5M")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY,2H5M")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY,1H5M")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H53")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H52")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H51")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY,3H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY,2H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY,1H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY,3H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY,2H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY,1H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H53")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H52")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H51")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H73")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H72")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H71")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H73")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H72")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H71")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H73")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H72")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H71")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H3 ")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string("THY, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("THY, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  A, H2 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("  A,2H6 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  A,1H6 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  A, H62")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  A, H61")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  A, H62")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  A, H61")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  A, H8 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("  A, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  A, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  A, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  A, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("A  , H2 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("A  ,2H6 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("A  ,1H6 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("A  , H62")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("A  , H61")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("A  , H62")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("A  , H61")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("A  , H8 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("A  , H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("A  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("A  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("A  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ADE, H62")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ADE, H61")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ADE, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ADE, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ADE, H2 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("ADE,2H6 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ADE,1H6 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ADE, H62")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ADE, H61")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("ADE, H8 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("ADE, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ADE, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  C, H6 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("  C, H5 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("  C,2H4 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  C,1H4 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  C, H42")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  C, H41")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  C, H42")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  C, H41")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  C, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  C, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  C, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  C, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("C  , H6 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("C  , H5 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("C  ,2H4 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("C  ,1H4 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("C  , H42")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("C  , H41")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("C  , H42")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("C  , H41")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("C  , H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("C  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("C  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("C  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("CYT, H6 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("CYT, H5 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("CYT,2H4 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("CYT,1H4 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("CYT, H42")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("CYT, H41")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("CYT, H42")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("CYT, H41")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("CYT, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("CYT, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("CYT, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("CYT, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  G,2H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  G,1H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  G, H22")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  G, H21")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  G, H22")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  G, H21")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("  G, H1 ")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string("  G, H8 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("  G, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  G, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  G, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("  G, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("G  ,2H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("G  ,1H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("G  , H22")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("G  , H21")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("G  , H22")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("G  , H21")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("G  , H1 ")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string("G  , H8 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("G  , H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("G  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("G  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("G  , H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GUA,2H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("GUA,1H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("GUA, H22")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("GUA, H21")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("GUA, H1 ")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string("GUA, H8 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("GUA, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GUA, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GUA, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GUA, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H6 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string(" DT,3H5M")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT,2H5M")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT,1H5M")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H53")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H52")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H51")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT,3H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT,2H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT,1H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT,3H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT,2H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT,1H5A")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H53")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H52")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H51")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H73")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H72")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H71")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H73")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H72")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H71")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H73")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H72")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H71")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H3 ")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string(" DT, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DT, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DA, H2 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string(" DA,2H6 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DA,1H6 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DA, H62")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DA, H61")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DA, H62")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DA, H61")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DA, H8 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string(" DA, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DA, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DA, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DA, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DC, H6 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string(" DC, H5 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string(" DC,2H4 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DC,1H4 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DC, H42")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DC, H41")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DC, H42")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DC, H41")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DC, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DC, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DC, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DC, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DG,2H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DG,1H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DG, H22")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DG, H21")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DG, H22")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DG, H21")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string(" DG, H1 ")] = NEW_RADIUS_HAP;
        myNewVdWaalsRadius[string(" DG, H8 ")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string(" DG, H1*")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DG, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DG, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string(" DG, H1'")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("AIB,3HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("AIB,2HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("AIB,1HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("AIB,3HB1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("AIB,2HB1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("AIB,1HB1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("AIB,HB23")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("AIB,HB22")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("AIB,HB21")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("AIB,HB13")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("AIB,HB12")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("AIB,HB11")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ABU,3HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ABU, HE2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ABU,2HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ABU,1HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ABU,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ABU,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ABU, HE2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ABU, HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ABU, HG1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ABU, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ABU, HB1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ACE,3HH3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ACE,2HH3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ACE,1HH3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ACE,3H  ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ACE,2H  ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ACE,1H  ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ACE, H3 ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ACE, H2 ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ACE, H1 ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ASX,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ASX,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ASX, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("ASX, HB1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLX,2HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLX,1HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLX,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLX,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLX, HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLX, HG1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLX, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("GLX, HB1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE,3HE ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE,2HE ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE,1HE ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE,2HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE,1HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE, HE3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE, HE2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE, HE1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE, HE3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE, HE2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE, HE1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE, HG3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE, HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE, HG3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE, HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("MSE, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PCA,2HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PCA,1HG ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PCA,2HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PCA,1HB ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PCA, HG3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PCA, HG2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PCA, HB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("PCA, HB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("NH2, H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("NH2, HN1")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("NH2, H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("NH2, HN1")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("NH2,2HN ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("NH2,1HN ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("NME,3HH3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("NME,2HH3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("NME,1HH3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("NME, H  ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("NME,3H  ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("NME,2H  ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("NME,1H  ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("NME,2HN ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("NME,1HN ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("NME, H3 ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("NME, H2 ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("NME, H1 ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("NME, HN2")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("NME, HN1")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("FOR, H  ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("FOR,2H  ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("FOR,1H  ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("FOR, H2 ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("FOR, H1 ")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,2HAD")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,1HAD")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,2HBD")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,1HBD")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,2HAA")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,1HAA")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,2HBA")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,1HBA")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,2HBC")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,1HBC")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM, HAC")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,2HBB")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,1HBB")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM, HAB")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,3HMD")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,2HMD")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,1HMD")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,3HMC")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,2HMC")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,1HMC")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,3HMB")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,2HMB")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,1HMB")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,3HMA")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,2HMA")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,1HMA")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM, HHD")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("HEM, HHC")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("HEM, HHB")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("HEM, HHA")] = NEW_RADIUS_HAR;
        myNewVdWaalsRadius[string("HEM,HAD2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HAD1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HBD2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HBD1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HAA2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HAA1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HBA2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HBA1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HBC2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HBC1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HBB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HBB1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HMD3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HMD2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HMD1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HMC3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HMC2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HMC1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HMB3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HMB2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HMB1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HMA3")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HMA2")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HEM,HMA1")] = NEW_RADIUS_H;
        myNewVdWaalsRadius[string("HOH, H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH, H1 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH, H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH, H1 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH, H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH, H1 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH, H2 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH, H1 ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH,2H  ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH,1H  ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH,2H  ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH,1H  ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH,2H  ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH,1H  ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH,2H  ")] = NEW_RADIUS_HPOL;
        myNewVdWaalsRadius[string("HOH,1H  ")] = NEW_RADIUS_HPOL;
    }
}

double Label::vanDerWaalsEnergy(const PCAtom a1, const PCAtom a2)
{   
    PotentialPair tmpPair;
    tmpPair.a = string(a1->element);
    tmpPair.b = string(a2->element);
    if(tmpPair.a == string("SE")) tmpPair.a = string(" S");
    if(tmpPair.b == string("SE")) tmpPair.b = string(" S");
    double e;
    for(int i = 0; i < atomPairs.size(); ++i) {
        if(EqualPair(atomPairs[i], tmpPair)) {
            e = atomPairs[i].epsilon;
            //printf("atomPairs[i].epsilon = %lf\n", atomPairs[i].epsilon);
            break;
        }
    }
    //printf("e = %lf\n", e);
    //double s = (myVdWaalsRadius[string(a1->element)] + myVdWaalsRadius[string(a2->element)])/1.122;
    double s = (vanDerWaalsRadius(a1) + 
                vanDerWaalsRadius(a2)) / 1.122;
    double r = norm(a1->x, a1->y, a1->z, 
                    a2->x, a2->y, a2->z); 
    return(4 * e * (pow(s/r, 12) - pow(s/r, 6)));
}

double Label::vanDerWaalsRadius(const PCAtom a)
{
    if(useOldH) {
        return myVdWaalsRadius[string(a->element)];
    }
    else {
        string resName = string(a->GetResName());
        string atomName = string(a->name);
        if(myNewVdWaalsRadius.find(resName+","+atomName) != myNewVdWaalsRadius.end()) {
            return myNewVdWaalsRadius[resName+","+atomName];
        }
        else {
            if(atomName == " C  ") {
                return 1.65;
            }
            else if(atomName == " H  ") {
                return 1.00;
            }
        }
    }
    return myVdWaalsRadius[string(a->element)];
}

double Label::energy(Node node1, Node node2) {
    double energy = 0.;
    if(node1.ResInd != node2.ResInd && abs(double(pcresidues[node1.ResInd]->GetSeqNum()) - double(pcresidues[node2.ResInd]->GetSeqNum())) <= 1.5) {
        //printf("enter\n");
        for(int i = 0; i < node1.atomNo; ++i) {
            for(int j = 0; j < node2.atomNo; ++j) {
                if(!isMChain(node1.atoms[i]) || !isMChain(node2.atoms[j])) {
                    energy += vanDerWaalsEnergy(node1.atoms[i], node2.atoms[j]);
                }
            }
        }
        return energy;
    }
    else {
        //printf("not enter\n");
        for(int i = 0; i < node1.atomNo; ++i) {
            for(int j = 0; j < node2.atomNo; ++j) {
                energy += vanDerWaalsEnergy(node1.atoms[i], node2.atoms[j]);
            }
        }
        return energy;
    }
    return 0;
}

int randInt(int a, int b) {
    return a + (rand() % (b - a + 1));
}

void randomPermute(vector<int>& arr) {
    for(int i = 0; i < arr.size(); ++i) {
        arr[i] = i;
    }
    for(int i = 0; i < arr.size() - 1; ++i) {
        swap(arr[i], arr[randInt(i, arr.size() - 1)]);
    }
}
bool Label::isMChain(PCAtom atom) {
    if(atom->CheckID("CA", "C", "*") || atom->CheckID("N", "N", "*") || atom->CheckID("H", "H", "*") || atom->CheckID("HA", "H", "*") || atom->CheckID("O", "O", "*") || atom->CheckID("C", "C", "*")) 
        return true;
    return false;
}

bool Label::ReadFile(string filename) {
	int RC,lcount;
	char S[500];
	InitMatType();
	PCMMDBManager myMMDBManager = NULL;
	myMMDBManager = new CMMDBManager();
	myMMDBManager->SetFlag(MMDBF_PrintCIFWarnings);
  	RC = myMMDBManager->ReadPDBASCII(filename.c_str());
	if(RC) {
		printf(" ***** ERROR #%i READ:\n\n %s\n\n",RC,GetErrorDescription(RC));
		myMMDBManager->GetInputBuffer(S,lcount);
    		if (lcount >= 0) 
      			printf("       LINE #%i:\n%s\n\n", lcount, S);
    		else if (lcount == -1)
      			printf("       CIF ITEM: %s\n\n", S);
    		delete myMMDBManager;	
    		return false;
	}
	myPCMMDBManager = myMMDBManager;
	return true;
}

vector<double> Label::Metropolis(std::vector<std::vector<double> >& metric, std::vector<std::vector<int> >& permutation, vector<vector<int> >& ResInd2NodeID, vector<int>& NumOfLoc) {
    printf("Enter MCMC\n");
    vector<double> energyList;
    vector<double> energies;
    int NumOfClusters = 0;
    vector<vector<int> > clusters;

    for(int i = 0; i < permutation.size(); ++i) {
        if(NumOfClusters < permutation[i].size()) {
            NumOfClusters = permutation[i].size();
        }
    }
    //printf("size = %d", permutation.size());
    for(int i = 0; i < permutation.size(); ++i) {
        for(int j = 0; j < permutation[i].size(); ++j) {
            //printf("%d ", permutation[i][j]);
        }
    }
    energies.resize(NumOfClusters, 0);
    clusters.resize(NumOfClusters);
    for(int ii = 0; ii < ResInd2NodeID.size(); ++ii) {
        //printf("%d\n", ResInd2NodeID[ii].size());
        if(ResInd2NodeID[ii].size() == 1) {
            for(int jj = 0; jj < NumOfClusters; ++jj) {
                clusters[jj].push_back(ResInd2NodeID[ii][0]);
            }
        }
        else {
            for(int jj = 0; jj < ResInd2NodeID[ii].size(); ++jj) {
                clusters[permutation[ii][jj]].push_back(ResInd2NodeID[ii][jj]);
            }
        }
    }
    for(int i = 0; i < NumOfClusters; ++i) {
        for(int j = 0; j < clusters[i].size(); ++j) {
            for(int k = j + 1; k < clusters[i].size(); ++k) {
                energies[i] += metric[clusters[i][j]][clusters[i][k]];
            }
        }
    }
    energyList.push_back(accumulate(energies.begin(), energies.end(), 0.));
    for(int i = 0; i < NumOfSim; ++i) {
        //printf("%d\n", i);
        vector<vector<int> > tmpPerm;
        tmpPerm.resize(permutation.size());
        tmpPerm[0].resize(permutation[0].size());
        randomPermute(tmpPerm[0]);
        for(int jj = 1; jj < permutation.size(); ++jj) {
            if(strcmp(pcresidues[jj - 1]->GetChainID(), pcresidues[jj]->GetChainID()) || permutation[jj - 1].size() != permutation[jj].size()) {
                tmpPerm[jj].resize(permutation[jj].size());
                randomPermute(tmpPerm[jj]);
            }
            else {
                tmpPerm[jj].resize(permutation[jj].size());
                tmpPerm[jj] = tmpPerm[jj - 1];
            }
        }
        for(int ii = 0; ii < tmpPerm.size(); ++ii) {
            for(int jj = 0; jj < tmpPerm[ii].size(); ++jj) {
                //printf("%d ", tmpPerm[ii][jj]);
            }
        }
        vector<vector<int> > tmpCluster;
        vector<double> tmpEnergies;
        tmpEnergies.resize(NumOfClusters, 0);
        tmpCluster.resize(NumOfClusters);
        for(int ii = 0; ii < ResInd2NodeID.size(); ++ii) {
            //printf("%d\n", ResInd2NodeID[ii].size());
            if(ResInd2NodeID[ii].size() == 1) {
                for(int jj = 0; jj < NumOfClusters; ++jj) {
                    tmpCluster[jj].push_back(ResInd2NodeID[ii][0]);
                }
            }
            else {
                for(int jj = 0; jj < ResInd2NodeID[ii].size(); ++jj) {
                    tmpCluster[tmpPerm[ii][jj]].push_back(ResInd2NodeID[ii][jj]);
                }
            }
        }
        
        //printf("metric.size = %d, clusters = %d\n", metric.size(), NumOfClusters);
        for(int ii = 0; ii < NumOfClusters; ++ii) {
            //printf("%d, %d\n", clusters[ii].size(), tmpCluster[ii].size());
            for(int jj = 0; jj < tmpCluster[ii].size(); ++jj) {
                //printf("%d ", tmpCluster[ii][jj]);
            }
        }
        for(int ii = 0; ii < NumOfClusters; ++ii) {
            for(int jj = 0; jj < tmpCluster[ii].size(); ++jj) {
                for(int kk = jj + 1; kk < tmpCluster[ii].size(); ++kk) {
                    //printf("(%d, %d)\n", tmpCluster[ii][jj], tmpCluster[ii][kk]);
                    tmpEnergies[ii] += metric[tmpCluster[ii][jj]][tmpCluster[ii][kk]];
                }
            }
        }
        //double u = drand();
        //printf("Difference = %lf\n", accumulate(energies.begin(), energies.end(), 0.) - accumulate(tmpEnergies.begin(), tmpEnergies.end(), 0.));
        if(0 <= 0.1 * accumulate(energies.begin(), energies.end(), 0.) - 0.1 * accumulate(tmpEnergies.begin(), tmpEnergies.end(), 0.)) {
            energies = tmpEnergies;
            clusters = tmpCluster;
            permutation = tmpPerm;
            energyList.push_back(accumulate(energies.begin(), energies.end(), 0.));
        }
        else {
            energyList.push_back(accumulate(energies.begin(), energies.end(), 0.));
        }
    }
    return energyList;
}

bool Label::MCMCLabeling(int MolNo) {
    if(MolNo < 1 || MolNo > myPCMMDBManager->GetNumberOfModels()) {
        printf("MolNo should be between 1 and %d\n", myPCMMDBManager->GetNumberOfModels());
        return false;
    }
    MolelNo = MolNo;
    PCModel pCModel = myPCMMDBManager->GetModel(MolelNo);
    int SelHnd = myPCMMDBManager->NewSelection();
    myPCMMDBManager->Select(
        SelHnd,
        STYPE_RESIDUE,
        MolelNo,
        "*",
        ANY_RES, "*",
        ANY_RES, "*",
        "!HOH",
        "*",
        "*",
        "*",
        SKEY_NEW
                    );
    PPCResidue pPCRes;
    int nPCRes;
    myPCMMDBManager->GetSelIndex(SelHnd, pPCRes, nPCRes);
    pcresidues.resize(nPCRes);
    for(int i = 0; i < nPCRes; ++i) {
        pcresidues[i] = pPCRes[i];
    }
    printf("Number of Residues = %d\n", nPCRes);
    
    std::vector<Node> nodes;
    std::vector<int> NumOfLoc;
    NumOfLoc.resize(nPCRes);
    for(int i = 0; i < nPCRes; ++i) {
        PPCAtom AtomTable;
        int NumOfAtoms;
        pcresidues[i]->GetAtomTable(AtomTable, NumOfAtoms);
        set<string> locs;
        for(int j = 0; j < NumOfAtoms; ++j) {
            string loc = string(AtomTable[j]->altLoc);
            locs.insert(loc);
        }
        if(locs.size() > 1) {
            locs.erase("");
        }
        NumOfLoc[i] = locs.size();
        for(set<string>::iterator ite = locs.begin(); ite != locs.end(); ++ite) {
            Node TmpNode;
            TmpNode.ResInd = i;
            strcpy(TmpNode.altLoc, (*ite).c_str());  
            nodes.push_back(TmpNode);          
        }
    }
    
    for(int i = 0; i < nodes.size(); ++i) {
        int SelHnd = myPCMMDBManager->NewSelection();
        //printf("%d, %s, %d\n", pcresidues[nodes[i].ResInd]->GetModelNum(), pcresidues[nodes[i].ResInd]->GetChainID(), pcresidues[nodes[i].ResInd]->GetSeqNum());
        myPCMMDBManager->Select(
            SelHnd,
            STYPE_ATOM,
            pcresidues[nodes[i].ResInd]->GetModelNum(),
            pcresidues[nodes[i].ResInd]->GetChainID(),
            pcresidues[nodes[i].ResInd]->GetSeqNum(), "*",
            pcresidues[nodes[i].ResInd]->GetSeqNum(), "*",
            "*",
            "*",
            "*",
            (string(",") + string(nodes[i].altLoc)).c_str(),
            SKEY_NEW
                    );
        PPCAtom ppCAtom = NULL;
        int atomNo;
        myPCMMDBManager->GetSelIndex(SelHnd, ppCAtom, atomNo);
        nodes[i].atoms = ppCAtom;
        nodes[i].atomNo = atomNo;
        
        //test part
        /*
        if(pcresidues[nodes[i].ResInd]->GetSeqNum() == 300)
            printf("%d\n", nodes[i].atomNo);
        */
    }
    
    vector<vector<int> > ResInd2NodeID;
    ResInd2NodeID.resize(pcresidues.size());
    for(int i = 0; i < nodes.size(); ++i) {
        ResInd2NodeID[nodes[i].ResInd].push_back(i);
    }
    
    vector<vector<int> > permutation;
    permutation.resize(pcresidues.size());
    for(int i = 0; i < permutation.size(); ++i) {
        permutation[i].resize(ResInd2NodeID[i].size());
        for(int j = 0; j < permutation[i].size(); ++j) {
            permutation[i][j] = j;
        }
    }
    
    vector<vector<double> > metric;
    metric.resize(nodes.size());
    for(int i = 0; i < metric.size(); ++i) {
        metric[i].resize(nodes.size());
        for(int j = 0; j < metric[i].size(); ++j) {
            metric[i][j] = DBL_MAX;
        }
    }
    for(int i = 0; i < metric.size(); ++i) {
        for(int j = i + 1; j < metric.size(); ++j) {
            //printf("%d, %d\n", i, j);
            if(nodes[i].ResInd != nodes[j].ResInd) {
            	metric[i][j] = energy(nodes[i], nodes[j]);
            	metric[j][i] = metric[i][j];
            }
            //printf("metric[%d][%d] = %lf\n", i, j, metric[i][j]);
        }
    }
    vector<vector<vector<int> > > perm;
    perm.resize(NumOfChains);
    vector<vector<double> > energyList;
    energyList.resize(NumOfChains);
    for(int i = 0; i < NumOfChains; ++i) {
        perm[i] = permutation;
        energyList[i] = Metropolis(metric, perm[i], ResInd2NodeID, NumOfLoc);
        //printf("%d\n", energyList[i].size());
        /*
        for(int i = 0; i < energyList.size(); ++i) {
            printf("%lf\n", energyList[i]);
        }
        */
    }
    int minInd;
    double minEnergy = DBL_MAX;
    for(int i = 0; i < NumOfChains; ++i) {
        if(minEnergy > energyList[i][energyList[i].size() - 1]) {
            minInd = i;
            minEnergy = energyList[i][energyList[i].size() - 1];
        }
    }
    for(int i = 0; i < energyList[minInd].size(); ++i) {
        printf("%lf\n", energyList[minInd][i]);
    }
    permutation = perm[minInd];
    
    map<int, string> Num2Str;
    map<string, int> Str2Num;
    for(int i = 0; i < 26; ++i) {
        Num2Str[i] = string(1, char(i + 65));
        Str2Num[string(1, char(i + 65))] = i;
    }
    
    for(int i = 0; i < permutation.size(); ++i) {
        for(int j = 0; j < permutation[i].size(); ++j) {
            printf("%d ", permutation[i][j]);
        }
        printf("\n");
    }
    
    for(int i = 0; i < pcresidues.size(); ++i) {
        if(ResInd2NodeID[i].size() > 1) {
            int SelHnd = myPCMMDBManager->NewSelection();
            myPCMMDBManager->Select(
                SelHnd,
                STYPE_ATOM,
                pcresidues[i]->GetModelNum(),
                pcresidues[i]->GetChainID(),
                pcresidues[i]->GetSeqNum(), pcresidues[i]->GetInsCode(),
                pcresidues[i]->GetSeqNum(), pcresidues[i]->GetInsCode(),
                "*",
                "*",
                "*",
                "*",
                SKEY_NEW
                        );
            PPCAtom ppCAtom = NULL;
            int atomNo;
            myPCMMDBManager->GetSelIndex(SelHnd, ppCAtom, atomNo);
            for(int j = 0; j < atomNo; ++j) {
                if(strcmp(ppCAtom[j]->altLoc, "")) {
                    strcpy(ppCAtom[j]->altLoc, Num2Str[permutation[i][Str2Num[string(ppCAtom[j]->altLoc)]]].c_str());
                }
            }
        }                   
    }
    myPCMMDBManager->WritePDBASCII(outputfile.c_str());
}

bool Label::SetParameter(int numsim, int numchain, string ofile, bool oldH) {
    NumOfSim = numsim;
    NumOfChains = numchain;
    outputfile = ofile;
    useOldH = oldH;
    return true;
}

int main(int argc, char** argv) {
    if(argc != 5 && argc != 6) {
        printf("Usage: Label in.pdb out.pdb num_sims num_markov_chains [use_old_h]\n");
        exit(1);
    }
    Label label;
    label.setupEpsilonVals("epsilon.xml");
    label.setupVanDerWaalsRadii();
    string PDB_filename, outfile;
    PDB_filename = string(argv[1]);
    outfile = string(argv[2]);
    int numsim;
    numsim = atoi(argv[3]);
    int numchain;
    numchain = atoi(argv[4]);
    bool oldH = false;
    if(argc == 6 && (string(argv[5]) == "T" || atoi(argv[5]) == 1)) {
        oldH = true;
        printf("Warning: Using hydrogen vdW radii from MMDB (~based on nuclear positions)");
        printf(" instead of from new Reduce (based on electron cloud centers)\n");
    }
    if(!oldH) {
        printf("Using hydrogen vdW radii from new Reduce (based on electron cloud centers)\n");
    }
    label.ReadFile(PDB_filename);  
    label.SetParameter(numsim, numchain, outfile, oldH);
    label.MCMCLabeling(1);
}

