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


/*
    LoopTK: Protein Loop Kinematic Toolkit
    Copyright (C) 2007 Stanford University

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "PConstants.h"
#include "myPDBIO.h"
#include "PExtension.h"
#include "PLibraries.h"
#include "PResources.h"
#include "PUtilities.h"

#include <iomanip>

/*
 * Public methods.
 */

PProtein* myPDBIO::readFromFile(const string &fileName)
{
  vector<string> fileLines = PUtilities::getLines(fileName),
                 atomLines = getAllAtomLines(fileLines);
  PProtein *result = readFromLines(atomLines);

  if (result->InAnyCollision()) {
    LoopTK::DisplayWarning("PDB file " + fileName + " contains one or more collisions.");
  }

  return result;
}
/*
void myPDBIO::writeToFile(PChain *chain, const string &fileName, const string &prefix)
{
  PResidue *curRes;
  HASH_MAP_STR(PAtom*) *curMap;

  int curAtomNum = 0;
  string curLine;
  ofstream outFile(fileName.c_str());

  if (chain == NULL) {
    PUtilities::AbortProgram("Error: Chain to be written to " + fileName + " is null.");
  }

  for(unsigned i = 0; i < chain->size(); i++) {
    curRes = chain->getResidue(i);
    curMap = curRes->getAtomMap();

    for(HASH_MAP_STR(PAtom *)::const_iterator it = curMap->begin(); it != curMap->end(); ++it) {
      curAtomNum++;
      curLine = "ATOM  ";

      addAtomNum(curLine, curAtomNum, prefix.length());
      addAtomName(curLine, it->first, prefix.length());
      addResName(curLine, curRes->getName(), prefix.length());
      addChainID(curLine, 'A', prefix.length());
      addResNum(curLine, i + 1, prefix.length());
      addInsertionCode(curLine, " ", prefix.length());
      addAtomPos(curLine, it->second->getPos(), prefix.length());
      addAuxData(curLine, it->second->getOccupancy(), it->second->getTempFactor(), "", prefix.length());
      addElemName(curLine, it->second->getName(), prefix.length());

      outFile << prefix << curLine << endl;
    }
  }
  outFile.close();
}
*/
void myPDBIO::writeToFile(PChain *chain, const string &fileName, const string& altLoc, const string &prefix)
{
  PResidue *curRes;
  HASH_MAP_STR(PAtom*) *curMap;

  int curAtomNum = 0;
  string curLine;
  ofstream outFile(fileName.c_str());

  if (chain == NULL) {
    PUtilities::AbortProgram("Error: Chain to be written to " + fileName + " is null.");
  }

  for(unsigned i = 0; i < chain->size(); i++) {
    curRes = chain->getResidue(i);
    curMap = curRes->getAtomMap();

    for(HASH_MAP_STR(PAtom *)::const_iterator it = curMap->begin(); it != curMap->end(); ++it) {
      curAtomNum++;
      curLine = "ATOM  ";

      addAtomNum(curLine, curAtomNum, prefix.length());
      addAtomName(curLine, it->first, prefix.length());
      addAltLoc (curLine, altLoc, prefix.length ());
      addResName(curLine, curRes->getName(), prefix.length());
      addChainID(curLine, 'A', prefix.length());
      addResNum(curLine, i + 1, prefix.length());
      addInsertionCode(curLine, " ", prefix.length());
      addAtomPos(curLine, it->second->getPos(), prefix.length());
      addAuxData(curLine, it->second->getOccupancy(), it->second->getTempFactor(), "", prefix.length());
      addElemName(curLine, it->second->getName(), prefix.length());

      outFile << prefix << curLine << endl;
    }
  }
  outFile.close();
}

vector<string> myPDBIO::getAllAtomLines(const vector<string> &pdbLines)
{
  vector<string> atomLines = PUtilities::getLinesStarting(pdbLines, "ATOM  ", false, "TER   ");
  vector<string> hetLines = PUtilities::getLinesStarting(pdbLines, "HETATM", false, "TER   ");

  atomLines.insert(atomLines.end(), hetLines.begin(), hetLines.end());
  return atomLines;
}

myPDBIO::pdbAtomLine myPDBIO::parseAtomLine(const string &line)
{
  pdbAtomLine curData;

  curData.atomNum = atoi(line.substr(6, 5).c_str());
  curData.atomName = line.substr(12, 4);
  PUtilities::trimSpaces(curData.atomName);

  curData.altLoc = line[16];
  curData.resName = line.substr(17, 3);

  curData.chainID = line[21];
  curData.resNum = atoi(line.substr(22, 4).c_str());
  curData.insCode = line[26];

  Real x = atof(line.substr(30, 8).c_str()),
       y = atof(line.substr(38, 8).c_str()), 
       z = atof(line.substr(46, 8).c_str());
  curData.pos = Vector3(x, y, z);
  
  if (line.length() > 54) {
    string occupancy = line.substr(54, 6);
    PUtilities::trimSpaces(occupancy);
    curData.occupancy = occupancy;

    if (line.length() > 60) {
      string tempFactor = line.substr(60, 6);
      PUtilities::trimSpaces(tempFactor);
      curData.tempFactor = tempFactor;

      if (line.length() > 76) {
        curData.elem = (line[76] == ' ' ? line.substr(77, 1) : line.substr(76, 2));
      }
    }
  }

  return curData;
}

/*
 * Private PDB-reading helper methods.
 */

bool myPDBIO::isCyclicPDB(const string &fileName)
{
  vector<string> keywordLines = PUtilities::getLinesStarting(fileName, "KEYWDS");

  for(unsigned i = 0; i < keywordLines.size(); i++) {
    if (keywordLines[i].find("CYCLIC BACKBONE") != string::npos ||
        keywordLines[i].find("CYCLIC PEPTIDE") != string::npos) {
      return true;
    }
  }

  return false;
}

PProtein* myPDBIO::readFromLines(const vector<string> &pdbLines)
{
  pdbAtomLine curData;
  hash_map<int, vector<pdbAtomLine> > atomMap;
  set<int> allResNums;

  if (pdbLines.size() == 0) {
    PUtilities::AbortProgram("Error: No valid \"ATOM\" lines to parse.");
  }

  for(unsigned i = 0; i < pdbLines.size(); i++) {
    curData = parseAtomLine(pdbLines[i]);
    if (PResources::ContainsAtomIDMapping(curData.atomName, curData.resName)) {
      curData.atomName = PResources::GetAtomIDMapping(curData.atomName, curData.resName);
    }
    if ((curData.altLoc == ' ' || curData.altLoc == 'A') && curData.insCode == ' ') {
      atomMap[curData.resNum].push_back(curData);
      allResNums.insert(curData.resNum);
    }
  }

  return readFromMap(atomMap, allResNums);  
}

PProtein* myPDBIO::readFromMap(const hash_map<int, vector<pdbAtomLine> > &atomMap, 
           const set<int> &allResNums)
{
  PProtein *protein = NULL;
  vector<int> resNums(allResNums.begin(), allResNums.end());

  sort(resNums.begin(), resNums.end());
  trimHeadAndTail(resNums, atomMap);

  for(unsigned i = 0; i < resNums.size(); i++) {
    vector<pdbAtomLine> curResIndex = atomMap.find(resNums[i])->second;
    PResidueSpec curResidueSpec = getSpec(curResIndex);

    updateProtein(protein, curResIndex[0].resName, curResidueSpec, resNums[i]);

    for (unsigned j = 0; j < curResIndex.size(); j++) {
      if (curResIndex[j].occupancy != "") {
        PAtom *curAtom = protein->getAtomAtRes(curResIndex[j].atomName, protein->size() - 1);
        if (curAtom != NULL) curAtom->setOccupancy(atof(curResIndex[j].occupancy.c_str()));
      }
      
      if (curResIndex[j].tempFactor != "") {
        PAtom *curAtom = protein->getAtomAtRes(curResIndex[j].atomName, protein->size() - 1);
        if (curAtom != NULL) curAtom->setTempFactor(atof(curResIndex[j].tempFactor.c_str()));
      }
    }
  }

  if (protein == NULL) {
    PUtilities::AbortProgram("Error: PDB file couldn't be loaded.");
  }
  
  protein->finalize();
  return protein;
}

void myPDBIO::updateProtein(PProtein* &protein, const string &resName, PResidueSpec &spec,
        int resNum)
{
  string shellToUse = resName;

  if (definesBackbone(spec)) {
    if (!PResources::ContainsResidueShell(resName)) {
      LoopTK::DisplayWarning("Residue \"" + resName + "\" (" + PUtilities::toStr(resNum) +
                             ") unknown; using only backbone atoms.");
      shellToUse = "backbone";
    }

    if (protein == NULL) protein = new PProtein(shellToUse, spec);
    else protein->AddResidue(shellToUse, spec);
  } else {
    PUtilities::AbortProgram("Error: Interior residue " + resName + " (" +
           PUtilities::toStr(resNum) + ") does not contain a backbone.");
  }
}

void myPDBIO::trimHeadAndTail(vector<int> &resIndices, const hash_map<int, vector<pdbAtomLine> > &atomMap)
{
  /* Trim head residues that don't have a backbone. */
  while(true) {
    if (resIndices.empty()) break;

    vector<pdbAtomLine> curLines = atomMap.find(resIndices[0])->second;
    PResidueSpec spec = getSpec(curLines);
    if (definesBackbone(spec)) break;

    LoopTK::DisplayWarning("Discarded head residue " + curLines[0].resName + " that lacks a backbone.");
    resIndices.erase(resIndices.begin());
  }

  /* Trim tail residues that don't have a backbone. */
  while(true) {
    if (resIndices.empty()) break;

    vector<pdbAtomLine> curLines = atomMap.find(resIndices[resIndices.size() - 1])->second;
    PResidueSpec spec = getSpec(curLines);
    if (definesBackbone(spec)) break;

    LoopTK::DisplayWarning("Discarded tail residue " + curLines[0].resName + " that lacks a backbone.");
    resIndices.pop_back();
  }
}

PResidueSpec myPDBIO::getSpec(const vector<pdbAtomLine> &atomLines)
{
  PResidueSpec spec;

  if (atomLines.size() > 0) {
    spec.pdb_id = atomLines[0].resNum;
  }

  for(unsigned i = 0; i < atomLines.size(); i++) {
    spec.Atom_Positions.addAtom(atomLines[i].atomName, atomLines[i].pos);
  }

  return spec;
}

bool myPDBIO::definesBackbone(PResidueSpec &spec) {
  return(spec.Atom_Positions.contains(PID::N) &&
         spec.Atom_Positions.contains(PID::C_ALPHA) &&
         spec.Atom_Positions.contains(PID::C));
}

/*
 * Private PDB-writing helper methods.
 */

void myPDBIO::addAtomNum(string &pdbOutLine, int atomNum, string::size_type prefixLen)
{
  /* Atom serial number occupies columns 7-11. */

  pdbOutLine += PUtilities::padLeft(PUtilities::toStr(atomNum), 5);
  assert(pdbOutLine.length() - prefixLen == 11);
  pdbOutLine += ' ';
}

void myPDBIO::addAtomName(string &pdbOutLine, const string &atomName, string::size_type prefixLen)
{
  /* Atom name, right-padded to length 3 then left-padded to length 4.  */

  pdbOutLine += PUtilities::padLeft(PUtilities::padRight(atomName, 3), 4);
  assert(pdbOutLine.length() - prefixLen == 16);
  //pdbOutLine += ' ';
}

void myPDBIO::addAltLoc(string &pdbOutLine, const string &altLoc, string::size_type prefixLen)
{
  /* Residue name, left-padded to have length 3. */

  pdbOutLine += altLoc.c_str ( );
  assert(pdbOutLine.length() - prefixLen == 17);
}

void myPDBIO::addResName(string &pdbOutLine, const string &resName, string::size_type prefixLen)
{
  /* Residue name, left-padded to have length 3. */

  pdbOutLine += PUtilities::padLeft(resName, 3);
  assert(pdbOutLine.length() - prefixLen == 20);
  pdbOutLine += ' ';
}

void myPDBIO::addChainID(string &pdbOutLine, char chainID, string::size_type prefixLen)
{
  /* Chain identifier. */

  pdbOutLine += chainID;
  assert(pdbOutLine.length() - prefixLen == 22);
}

void myPDBIO::addResNum(string &pdbOutLine, int resNum, string::size_type prefixLen)
{
  /* Residue sequence number. */

  pdbOutLine += PUtilities::padLeft(PUtilities::toStr(resNum), 4);
  assert(pdbOutLine.length() - prefixLen == 26);
}

void myPDBIO::addInsertionCode(string &pdbOutLine, const string &insCode, string::size_type prefixLen)
{
  /* Code for insertion of residues, column 27. */

  pdbOutLine += insCode;
  pdbOutLine += "   ";  /* Columns 28-30 unused. */
  assert(pdbOutLine.length() - prefixLen == 30);
}

void myPDBIO::addAtomPos(string &pdbOutLine, const Vector3 &atomPos, string::size_type prefixLen)
{
  /* Atom position (x, y, z) occupies columns 31-38, 39-46 and 47-54,
   * respectively.  Each component of the vector is left-padded. */
  char buf[64];

  sprintf(buf, "%.4f", atomPos.x);
  pdbOutLine += PUtilities::padLeft(string(buf), 8);
  assert(pdbOutLine.length() - prefixLen == 38);

  sprintf(buf, "%.4f", atomPos.y);
  pdbOutLine += PUtilities::padLeft(string(buf), 8);
  assert(pdbOutLine.length() - prefixLen == 46);

  sprintf(buf, "%.4f", atomPos.z);
  pdbOutLine += PUtilities::padLeft(string(buf), 8);
  assert(pdbOutLine.length() - prefixLen == 54);
}

void myPDBIO::addAuxData(string &pdbOutLine, Real occupancy, Real tempFactor,
      const string &segmentID, string::size_type prefixLen)
{
  stringstream realConverter;
  /* Codes for occupancy, temperature factor and segment ID. */

  realConverter << setiosflags(ios::fixed) << setprecision(2) << occupancy;
  pdbOutLine += PUtilities::padLeft(realConverter.str(), 6);
  assert(pdbOutLine.length() - prefixLen == 60);

  if (tempFactor > 0) {
    realConverter.str("");
    realConverter << setiosflags(ios::fixed) << setprecision(2) << tempFactor;
    pdbOutLine += PUtilities::padLeft(realConverter.str(), 6);
  } else {
    pdbOutLine += PUtilities::padLeft("", 6);
  }
  assert(pdbOutLine.length() - prefixLen == 66);

  pdbOutLine += PUtilities::padLeft("", 6);
  assert(pdbOutLine.length() - prefixLen == 72);

  pdbOutLine += PUtilities::padLeft(segmentID, 4);
  assert(pdbOutLine.length() - prefixLen == 76);
}

void myPDBIO::addElemName(string &pdbOutLine, const string &elemName, string::size_type prefixLen)
{
  /* Element name, left-padded to have length 2. */

  pdbOutLine += PUtilities::padLeft(elemName, 2);
  assert(pdbOutLine.length() - prefixLen == 78);
}
