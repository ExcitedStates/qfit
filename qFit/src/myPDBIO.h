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

/*
 * PExtension defines a variety of auxiliary
 * classes that are not considered part of the
 * "core" LoopTK class hierarchy. These include
 * file I/O classes, initialization classes,
 * test classes, and so forth.
 */

#ifndef MYPDBIO_H
#define MYPDBIO_H

#include "PBasic.h"
#include "PIKAlgorithms.h"
#include "PConstants.h"
// @package Conformation Analysis

/**
 *
 * Defines a space of protein loop conformations,
 * as well as the static portions of the protein
 * before and after the loop.
 */

class myPDBIO {
  friend class CS2IO;	/* To access readFromLines(). */

  public:

    /**
     * A convenience structure that encapsulates
     * the data on one line of a PDB "ATOM" descriptor.
     */

    struct pdbAtomLine {
      string atomName, resName, elem;
      char altLoc, insCode, chainID;
      int atomNum, resNum;
      string occupancy, tempFactor;
      Vector3 pos;
    };

    /**
     * Reads the specified PDB file and returns the
     * contents as a PProtein pointer, or NULL if an
     * I/O error is encountered.
     */

    static PProtein* readFromFile(const string &fileName);

    /**
     * Writes the specified PChain <code>chain</code> to a file named
     * <code>fileName</code>, using the standard PDB format.  Optionally,
     * a <code>prefix</code> may be appended to each line of output.
     */

    //static void writeToFile(PChain *chain, const string &fileName, const string &prefix = "");
    static void writeToFile(PChain *chain, const string &fileName, const string& altLoc = " ", const string &prefix = "");

    /**
     * Returns <code>true</code> if the PDB file specified by
     * <code>fileName</code> appears to represent a cyclic protein,
     * <code>false</code> otherwise.
     */

    static bool isCyclicPDB(const string &fileName);

    /**
     * Returns a <code>vector</code> containing all the lines in
     * <code>pdbLines</code> that begin with <code>"ATOM  "</code>
     * or <code>"HETATM"</code>, up to the <code>"TER"</code>
     * terminating line (if one is present).
     */

    static vector<string> getAllAtomLines(const vector<string> &pdbLines);

    /**
     * Returns a pdbAtomLine structure containing
     * data from the specified line in a PDB file.
     */

    static pdbAtomLine parseAtomLine(const string &line);

  private:

    /* Helper functions for readFromFile(). */

    static PProtein* readFromLines(const vector<string> &pdbLines);
    static PProtein* readFromMap(const hash_map<int, vector<pdbAtomLine> > &atomMap, 
				 const set<int> &allResNums);
    static void updateProtein(PProtein* &protein, const string &resName,
				PResidueSpec &spec, int resNum);

    static void trimHeadAndTail(vector<int> &resIndices, const hash_map<int,
				vector<pdbAtomLine> > &atomMap);
    static PResidueSpec getSpec(const vector<pdbAtomLine> &atomLines);
    static bool definesBackbone(PResidueSpec &spec);

    /* Helper functions for writeToFile() that do some error checking. */

    static void addAtomNum(string &pdbOutLine, int atomNum, string::size_type prefixLen);
    static void addAtomName(string &pdbOutLine, const string &atomName, string::size_type prefixLen);
    static void addResName(string &pdbOutLine, const string &resName, string::size_type prefixLen);
    static void addAltLoc(string &pdbOutLine, const string &altLoc, string::size_type prefixLen);
    static void addChainID(string &pdbOutLine, char chainID, string::size_type prefixLen);
    static void addResNum(string &pdbOutLine, int resNum, string::size_type prefixLen);
    static void addInsertionCode(string &pdbOutLine, const string &insCode, string::size_type prefixLen);
    static void addAtomPos(string &pdbOutLine, const Vector3 &atomPos, string::size_type prefixLen);
    static void addAuxData(string &pdbOutLine, Real occupancy, Real tempFactor,
				const string &segmentID, string::size_type prefixLen);
    static void addElemName(string &pdbOutLine, const string &elemName, string::size_type prefixLen);

};


#endif
