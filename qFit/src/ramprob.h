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


#ifndef __RAMPROB_H_
#define __RAMPROB_H_

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
#endif //__RAMPROB_H_

class RamProb
{
	public:
		static void loadRamPlot(ifstream &fin,int ramplot);
		static void loadAllRamPlots();
		static double getRamProb( double phi, double psi,const char* res,const char* nres);
	private:
		static const char* ramFiles[5];
		//static vector<vector<double> > ramplots[5];
		static double ramplots[5][180][180];
};
