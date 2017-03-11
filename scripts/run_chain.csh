#!/bin/tcsh

setenv DOKKUMDIR /data/vdbedem/qFitRun/DOKKUMDIR
set scriptsdir = /data/vdbedem/qFitRun
source /home/sdcsoftware/linux/ccp4-6.0.2/include/ccp4.setup

if ( ! -f chains.txt ) then
  echo "File chains.txt not found."
  exit
endif

set chaininfo = `awk '{if (NR=='$LSB_JOBINDEX') print $0 }' chains.txt`

set pdb = "$cwd/`basename $1`"
set mtz = "$cwd/`basename $2`"
set chain = $chaininfo[1]
set start = $chaininfo[2]
set stop = $chaininfo[3]
set F = $3
set SF = $4
set elbow = $5
echo $pdb
@ start = $start + 3
@ stop = $stop - 3

if ( $start <= 0 ) then
  set start = 1
endif
echo $start $stop
if ( ! -d $chain ) then
  mkdir $chain
endif

cd $chain

# RUN DOKKUM ON RESIDUES

set tmpstr  = `mktemp`
set bs = `basename $tmpstr`
set job_base = `echo $bs | cut -c 5-10`

bsub -J "${job_base}[$start-$stop]" -qsdcq -o/tmp/${job_base}.%I ${scriptsdir}/dokkum.csh ${pdb} ${mtz} $chain refinement.input.xray_data.labels=$F,$SF $elbow

# WAIT FOR JOBS TO FINISH AND BUILD MODEL FROM RESIDUES

rm -f ALL.pdb
bsub -K -qsdcq -o/tmp/LSF.%J -w "ended(${job_base})" ${scriptsdir}/get_Res.csh $pdb $chain $start $stop

grep -v END ALL.pdb > ALL_.pdb

pdbset xyzin ALL_.pdb xyzout ALL_.pdb << eopdb > pdbset.log
exclude hydrogens
chain $chain
end
eopdb

${DOKKUMDIR}/TRUNCATE/idmulti 0.1 ALL_.pdb ALL_idmulti.pdb $chain > fragids.txt

awk '{sub(/ SE /," SD ");print}' ALL_idmulti.pdb > ALL_idmulti_SD.pdb

fragdokkum:

set tmpstr  = `mktemp`
set bs = `basename $tmpstr`
set job_base = `echo $bs | cut -c 5-10`
set basepdb = `basename $pdb .pdb`
set NL = `wc fragids.txt`

bsub -J "${job_base}_frag[1-$NL[1]]" -qsdcq -o/tmp/${job_base}_frag.%I ${scriptsdir}/fragDokkum.csh ../../${basepdb}_refine_001.pdb ../../${basepdb}_refine_001.mtz $chain ../ALL_idmulti_SD.pdb

bsub -K -qsdcq -o/tmp/LSF.%J -w "ended(${job_base}_frag)" ${scriptsdir}/insert_frags.csh ALL_idmulti.pdb $chain

# MAKE COMPLETE CHAIN
makechain:
awk '{ if(substr($0,22,1)=="'$chain'"&&substr($0,23,4)+0<'$start') print $0}' ../noheader.pdb > final.pdb
grep -v END ALL_idmulti.pdb >> final.pdb
awk '{ if(substr($0,22,1)=="'$chain'"&&substr($0,23,4)+0>'$stop') print $0}' ../noheader.pdb >> final.pdb
