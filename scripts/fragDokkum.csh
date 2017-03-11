#!/bin/tcsh

setenv DOKKUMDIR /data/vdbedem/qFitRun/DOKKUMDIR
setenv LOOPTK_RESOURCEDIR ${DOKKUMDIR}/LOOPTK/resources/
source /home/sw/rhel5/x86_64/phenix/phenix-1.8-1069/phenix_env
setenv LD_LIBRARY_PATH ${DOKKUMDIR}/coin-Cbc/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /home/vdbedem//ILOG/CPLEX_Studio_AcademicResearch122/cplex/lib/x86-64_sles10_4.1/static_pic:/home/vdbedem//ILOG/CPLEX_Studio_AcademicResearch122/concert/lib/x86-64_sles10_4.1/static_pic:$LD_LIBRARY_PATH
source /home/sdcsoftware/linux/ccp4-6.0.2/include/ccp4.setup
setenv ILOG_LICENSE_FILE /home/vdbedem/DEV/ModelBuildTool/MBT/QFIT/ILOG/CPLEX_Studio_AcademicResearch122/licenses/access.ilm

set pdb_scale = $1
set mtz = $2
set chain = $3
set pdb = $4

set fragstart = `awk '{if (NR=='$LSB_JOBINDEX') print $1 }' fragids.txt`
set fragend = `awk '{if (NR=='$LSB_JOBINDEX') print $NF }' fragids.txt`

set jobdir = `awk '{if (NR=='$LSB_JOBINDEX') { print $1 "_" $NF }}' fragids.txt`

if ( ! -d $jobdir ) then
  mkdir $jobdir
endif

cd $jobdir

${DOKKUMDIR}/reDokkum_frag $pdb_scale $mtz $chain $fragstart $fragend $pdb > reDokkum_frag.log

#awk BEGIN'{ q="\047"; }  {if (NR==2) {gsub(" ","_"); cmd = "mkdir -p "q$0q; system(cmd) }}' fragids.txt
