#!/bin/tcsh

set bn = `basename $2 .pdb`

if ( -f tmp.txt )then
  rm tmp.txt
endif

/data/vdbedem/LX_ACONIO << eor > aconio1.log
1
${bn}.pdb
m2
eor

set m2s = m2_[1-6].pdb

foreach m ($m2s) 
  echo $m >> tmp.txt
end

/data/vdbedem/LX_ACONIO << eor > aconio2.log
3
`cat tmp.txt`

ABCDEF
N
1.02
-1.000
m2.pdb
eor

source /home/sw/rhel4/x86_64/phenix/phenix-1.3-final/phenix_env
phenix.refine $1 "m2.pdb" "main.number_of_macro_cycles=0" $3 "strategy=occupancies" "--overwrite"

mv m2_refine_001.pdb ${bn}_preprune.pdb
