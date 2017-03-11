#!/bin/tcsh

set pdb = $1
set chain = $2
set start = $3
set stop = $4

set i = $start

while ( $i <= $stop )

  cd $i
  if (! -f Res$i.pdb && ! -f Res.pdb ) then
    awk '{ if((substr($0,1,4)=="ATOM"||substr($0,1,6)=="HETATM")&&substr($0,22,1)=="'$chain'"&&substr($0,24,3)+0=='$i') print $0}' $pdb >> ../ALL.pdb
  else
    cat Res*.pdb >> ../ALL.pdb
  endif
  cd ..
  @ i = $i + 1
end
