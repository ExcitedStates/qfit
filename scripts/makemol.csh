#!/bin/tcsh

set pdb = $1
set start = $2
set stop = $3


#pdbset xyzin $pdb xyzout noheader.pdb << eopdbset > pdbset.log
#exclude headers
#eopdbset

egrep '^ATOM|^HETATM' $pdb > noheader.pdb

grep CRYST1 $pdb > final.pdb
awk '{ if(substr($0,22,1)=="A"&&substr($0,24,3)+0<'$start') print $0}' noheader.pdb >> final.pdb
grep -v END ALL_idmulti.pdb >> final.pdb
awk '{ if(substr($0,22,1)=="A"&&substr($0,24,3)+0>'$stop') print $0}' noheader.pdb >> final.pdb
awk '{ if(substr($0,22,1)!="A") print $0}' noheader.pdb >> final.pdb
