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


#ifndef PSCALEXMAP_H_
#define PSCALEXMAP_H_

#include <cstdlib>
#include <vector>
#include <list>
#include <algorithm>
#include <numeric>
#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

/* Sorting index table adapted from:
*  C. Hicks. "Creating an Index Table in STL" (C/C++ Users Journal, 18(8), August 2000).       
*/

template<class T, class U> 
struct ComparePair1stAsc
{
	bool operator()( const std::pair<T,U>& a, const std::pair<T,U>& b ) const 
	{
		return *a.first < *b.first; 
	}
};


template<class IterIn, class IterOut> 
void sort_idxtbl( IterIn first, IterIn last, IterOut out )
{
	std::vector< std::pair<IterIn,std::size_t> > s( last-first );
	for( std::size_t i=0; i < s.size(); ++i )
		s[i] = std::make_pair( first+i, i );
	std::stable_sort( s.begin(), s.end(), ComparePair1stAsc<IterIn,std::size_t>() );
	for( std::size_t i=0; i < s.size(); ++i, ++out )
		*out = s[i].second;
}


class PScaleXMap
{
public:

    PScaleXMap( );

    virtual ~PScaleXMap();
	
    bool DoScaleXMap(
            clipper::Xmap<clipper::ftype32>& xm, 
            clipper::HKL_info& hkl_info, 
            const std::string inmtz, 
            const clipper::Atom_list& atmList, 
            clipper::ftype32& m, 
            clipper::ftype32& va
            );

    bool L1Fit(
            const std::vector<double>&, 
            const std::vector<double>&, 
            clipper::ftype32&, 
            clipper::ftype32&
            ) const;

    bool L2Fit(
            const std::vector<double>&, 
            const std::vector<double>&, 
            clipper::ftype32&, 
            clipper::ftype32&
            ) const;

private:

    clipper::ftype32 RADIUS;
    clipper::ftype32 BSMEAR;

    std::vector<clipper::Xmap<clipper::ftype32>::Map_reference_coord> 
        CreateCarve(
                const clipper::Xmap<clipper::ftype32>&, 
                const clipper::Atom_list&
                ) const;

    bool Stats(
            const std::vector<double>& xmap_o, 
            clipper::ftype32& m, 
            clipper::ftype32& var
            );
};

#endif /*CSCALEXMAP_H_*/
