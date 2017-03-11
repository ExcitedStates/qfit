#!/bin/tcsh

source /home/sw/rhel5/x86_64/phenix/phenix-1.7.1-743/phenix_env

set cryst1 = `grep CRYST1 $1`
echo $cryst1

set i = 8
set cell = ""
while ( $i < ${#cryst1} )
	set c = $cryst1[$i]
	set cell = $cell$c
	@ i = $i + 1
end
phenix.cif_as_mtz $2 --unit-cell=$cryst1[2],$cryst1[3],$cryst1[4],$cryst1[5],$cryst1[6],$cryst1[7] --space-group=$cell --show-details-if-error
