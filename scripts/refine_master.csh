#!/bin/tcsh

source /home/sw/rhel5/x86_64/phenix/phenix-1.8-1069/phenix_env
#source /home/sdcsoftware/linux/ccp4-6.0.2/include/ccp4.setup
setenv DOKKUMDIR /data/vdbedem/qFitRun/DOKKUMDIR
set scriptsdir = /data/vdbedem/qFitRun

set pdb = $1
set mtz = $2
set F = $3
set SF = $4
set elbow = $5

set bspdb = `basename $pdb .pdb`

# DETERMINE RESOLUTION AND (AN)ISOTROPIC REFINEMENT
set resrange = `phenix.mtz.dump $mtz | grep "Resolution range:"`
set res = `echo $resrange[4] | cut -c 1-5`
set res1000 = `echo $res | awk '{tot = $1*1000}{print tot }'`

if ( $res1000 < 1550 ) then
  set adp = 'adp.individual.anisotropic="not (water or element H)"' 
else
  set adp = 'adp.individual.isotropic=all'
endif

phenix.refine $mtz $pdb ${scriptsdir}/1.params output.prefix=${bspdb} output.serial=2 --overwrite refinement.input.xray_data.labels=$F,$SF $elbow write_maps=false

#phenix.refine $mtz ${bspdb}_001.pdb refine.strategy=individual_sites refine.sites.torsion_angles=all output.prefix=${bspdb} output.serial=2 --overwrite  refinement.input.xray_data.labels=$F,$SF $elbow write_maps=false

set greptmp = `grep COORD ${bspdb}_002.pdb`
set dist = $greptmp[8]
set dist = `echo $dist | awk '{print $1*2.0}'`
set chns = `awk '{printf "%s ", $1 }' chains.txt`

cp ${bspdb}_002.pdb ${bspdb}_002_cluster.pdb

foreach c ( $chns )
  ln -s $c/fragids.txt
  ${DOKKUMDIR}/TRUNCATE/cluster ${bspdb}_002_cluster.pdb ${bspdb}_002_cluster.pdb $c $dist > cluster_$c.log
  rm -f fragids.txt
end
cp ${bspdb}_002_cluster.pdb ${bspdb}_002_cluster_unit.pdb

foreach c ( $chns )
	${DOKKUMDIR}/TRUNCATE/occ2unit ${bspdb}_002_cluster_unit.pdb ${bspdb}_002_cluster_unit.pdb $c
end

final:
echo phenix.refine $mtz ${bspdb}_002_cluster_unit.pdb ${scriptsdir}/3.params $adp --overwrite refinement.input.xray_data.labels=$F,$SF write_maps=false

phenix.refine $mtz ${bspdb}_002_cluster_unit.pdb ${scriptsdir}/3.params "$adp" --overwrite refinement.input.xray_data.labels=$F,$SF $elbow write_maps=false 
final2:

awk '{ if(substr($0,1,6)=="CRYST1"||(substr($0,1,4)=="ATOM"||substr($0,1,6)=="HETATM")&&substr($0,57,4)+0>0.19) print $0}' ${bspdb}_002_cluster_unit_refine_001.pdb > ${bspdb}_002_cluster_unit_refine_001_OCC.pdb

#pdbset xyzin ${bspdb}_002_cluster_unit_refine_001.pdb xyzout ${bspdb}_002_cluster_unit_refine_001_OCC.pdb << eopdbset > pdbsetOCC.log
#select occupancy 0.19
#eopdbset

foreach c ( $chns )
	${DOKKUMDIR}/TRUNCATE/occ2unit ${bspdb}_002_cluster_unit_refine_001_OCC.pdb ${bspdb}_002_cluster_unit_refine_001_OCC.pdb $c
end

phenix.refine $mtz ${bspdb}_002_cluster_unit_refine_001_OCC.pdb ${scriptsdir}/3.params "$adp" --overwrite refinement.input.xray_data.labels=$F,$SF $elbow write_maps=false 

cp ${bspdb}_002_cluster_unit_refine_001_OCC_refine_001.pdb qFit.pdb
cp ${bspdb}_002_cluster_unit_refine_001_OCC_refine_001.mtz qFit.mtz
