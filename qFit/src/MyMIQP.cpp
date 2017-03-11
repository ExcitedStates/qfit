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


#include <fstream>
#include <iostream>
#include <cstring>
#include <ilcplex/ilocplex.h>
#include "MyMIQP.hpp"

ILOSTLBEGIN

void MyMIQP(double** A, double* B, int n, double* s, double threshold) {
	IloEnv env;
	try {
		IloModel model(env);
		IloNumVarArray var(env);
		IloRangeArray con(env);
		for(int i = 0; i < n - 1; ++i) {
			var.add(IloNumVar(env, 0.0, 1.0));
		}
		var.add(IloNumVar(env, -1.0, 1.0)); /*Changed here from -40, 40*/
		for(int i = 0; i < n - 1; ++i) {
			var.add(IloNumVar(env, 0.0, 1.0, ILOINT));
		}
		IloExpr tmp(env);
		for(int i = 0; i < n; ++i) {
			tmp += B[i] * var[i];
		}
		for(int i = 0; i < n; ++i) {
			for(int j = 0; j < n; ++j) {
				tmp += A[i][j] * var[i] * var[j] / 2;
			}
		}
		/*
		for(int i = 0; i < n; ++i) {
			tmp += A[i][i] * var[i] * var[i] / 2;
		}
		*/
		model.add(IloMinimize(env, tmp));
		for(int i = 0; i < n - 1; ++i) {
			con.add(threshold * var[n + i] - var[i] <= 0);
			con.add(var[i] - var[n + i] <= 0);
		}
		IloExpr tmmm(env);
		for(int i = 0; i < n - 1; ++i) {
			tmmm += var[i];
		}
		con.add(tmmm <= 1);
		model.add(con);
		
		IloCplex cplex(model);
		cplex.solve();
		
	      	env.out() << "Solution status = " << cplex.getStatus() << endl;
	      	env.out() << "Solution value  = " << cplex.getObjValue() << endl;

	      	IloNumArray vals(env);
	      	cplex.getValues(vals, var);
	      	env.out() << "Values        = " << vals << endl;
	      	for(int i = 0; i < 2 * n - 1; ++i) {
	      		s[i] = vals[i];
	      	}
	      	cplex.getSlacks(vals, con);
	      	env.out() << "Slacks        = " << vals << endl;

	     	cplex.exportModel("miqpex1.lp");
   	}
   	catch (IloException& e) {
    	  cerr << "Concert exception caught: " << e << endl;
   	}
   	catch (...) {
    	  cerr << "Unknown exception caught" << endl;
   	}

   	env.end();
}
/*
int main() {
   	double** A = NULL;
   	int col;
   	std::cin >> col;
   	double* B = new double[col];
   	for(int i = 0; i < col; ++i) {
  		 A = new double* [col];
   	}
   	std::cout << "begins\n";
   	std::ifstream out("data.txt");
   	for(int i = 0; i < col; ++i) {
  		 A[i] = new double[col];
  		 for(int j = 0; j < col; ++j) {
  		   out >> A[i][j];
  		   std::cout << A[i][j] << " ";
  		 }
  		 std::cout << "\n";
   	}
   	for(int i = 0; i < col; ++i) {
  		 out >> B[i];
                 B[i] = 3 * B[i] - 3;
  		 std::cout << "B[i]" << B[i] << " ";
   	}
   	for(int i = 0; i < col; ++i) {
   		std::cout << B[i] << std::endl;
   	}
   	double* s = new double[col * 2 - 1];
   	MyMIQP(A, B, col, s, .1);
   	for(int i = 0; i < 2 * col - 1; ++i) {
   		std::cout << s[i] << " ";
   	}
   	return 0;
}
*/
