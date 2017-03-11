#!/bin/tcsh -f

set fls = [1-9]*_[1-9]*

set i = 1

while ( $i <= ${#fls} )
if ( -d ${fls[$i]} ) then
echo ${fls[$i]}
pdbset xyzin ${fls[$i]}/frag.pdb xyzout ${fls[$i]}_frag.pdb << eof_1
chain $2
eof_1
endif
@ i = $i + 1
end

set bn = `basename ${1} .pdb`

set i = 1
while ( $i <= ${#fls} )
if ( -d ${fls[$i]} ) then
/home/sdcsoftware/linux/CNSNSS/CombineFrags_v0.2 $1 ${fls[$i]}_frag.pdb ${fls[$i]}_frag.pdb
cp ${bn}_${fls[$i]}_frag.pdb ${bn}.pdb
endif
@ i = $i + 1
end
