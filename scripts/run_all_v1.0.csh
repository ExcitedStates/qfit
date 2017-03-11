#!/bin/tcsh

source /home/sdcsoftware/linux/ccp4-6.0.2/include/ccp4.setup
set scriptsdir = /data/vdbedem/qFitRun
source /home/sw/rhel5/x86_64/phenix/phenix-1.8-1069/phenix_env


set pdb = $2
set mtz = $3
set F = $4
set SF = $5

cd $1

echo $2 $3 $4 $5 > qFit.log

mv $pdb $pdb.org

phenix.reduce -Trim -NOFLIP $pdb.org | grep -v USER > $pdb

#pdbset xyzin $pdb.org xyzout $pdb << eopdbset > pdbset.log
#exclude hydrogen
#eopdbset

${scriptsdir}/DOKKUMDIR/TRUNCATE/renumberInsertions $pdb

phenix.elbow --do-all $pdb
set basepdb = `basename $2 .pdb`
echo $basepdb
set elbow = ""
echo "elbow:" $elbow
if ( -r "elbow.${basepdb}_pdb.all.001.cif" ) then
  set elbow = "$cwd/elbow.${basepdb}_pdb.*001.cif"
endif
echo "elbow:" $elbow

$scriptsdir/DOKKUMDIR/TRUNCATE/chains $pdb > chains.txt 

# CREATE PDB AND MTZ FOR SCALING FRAGDOKKUM LATER 

${scriptsdir}/fragDokkum_mtz.csh ${pdb} ${mtz} $F $SF $elbow > fragDokkum_mtz.log

if ( ! -f ${basepdb}_refine_001.pdb ) then
  goto makelog
endif

egrep '^ATOM|^HETATM' $pdb > noheader.pdb

#pdbset xyzin $pdb xyzout noheader.pdb << eopdbset >> pdbset.log
#exclude headers
#eopdbset

set nchns = `wc chains.txt`

set tmpstr  = `mktemp`
set bs = `basename $tmpstr`
set job_base = `echo $bs | cut -c 5-10`
set njobs = $nchns[1]

bsub -J "${job_base}[1-$njobs]" -qsdcq -Rlinux64 -o/tmp/${job_base}.%I ${scriptsdir}/run_chain.csh ${pdb} ${mtz} $F $SF $elbow

# PUT THE MOLECULE TOGETHER

bsub -K -qsdcq -o/tmp/LSF.%J -w "ended(${job_base})" "grep CRYST1 $pdb > final.pdb"
set chns = `awk '{printf "%s ", $1 }' chains.txt`
foreach c ( $chns )
  grep -v END $c/final.pdb >> final.pdb
end
set mtch = `echo $chns | awk '{gsub(" ", "")}1'`
awk '{chains = "['$mtch']"} { if(substr($0,22,1) !~ chains ) print $0}' noheader.pdb >> final.pdb

tst:
mv final.pdb tmp.pdb
awk '{if(substr($0,18,3) ~ "MET" && substr($0,11,4) ~ "[0-9] SE" ) {gsub("   SE","    S") gsub("SE  "," SD ") } print $0}' tmp.pdb > final.pdb

${scriptsdir}/runlabel.csh final.pdb final_relabel.pdb

# REFINE
refine:
${scriptsdir}/refine_master.csh final_relabel.pdb ${mtz} $F $SF $elbow

makelog:
grep -i sorry fragDokkum_mtz.log >> qFit.log
