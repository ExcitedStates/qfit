#!/bin/bash

pdb=$1
mtz=$2

grep CRYST1 $pdb > final.pdb # only place $pdb is used here!

chns=`awk '{printf "%s ", $1}' chains.txt`
for c in $chns; do
    grep -v END $c/final.pdb | grep -v "USER  MOD" >> final.pdb
    #rm -rf $c  # clean up
done
mtch=`echo $chns | awk '{gsub(" ", "")}1'`
awk '{chains="['$mtch']"} {if(substr($0,22,1) !~ chains) print $0}' noheader.pdb >> final.pdb
#rm noheader.pdb

mv final.pdb tmp.pdb
awk '{if(substr($0,18,3) ~ "MET" && substr($0,11,4) ~ "[0-9] SE" ) {gsub("   SE","    S") gsub("SE  "," SD ") } print $0}' tmp.pdb > final.pdb
rm tmp.pdb

phenix.python ${SCRIPTSDIR}/fix_1of2rsd_flips.py final.pdb
mv final_fixedPepFlipGeom.pdb final.pdb

phenix.reduce -q final.pdb | grep -v USER > final_eH.pdb
cp final_eH.pdb preLabel.pdb
python $SCRIPTSDIR/fix_dupe_waters.py final_eH.pdb > final_eH_fixHOH.pdb
ln -s $QFITDIR/Label/epsilon.xml
$QFITDIR/Label/Label final_eH_fixHOH.pdb final_eH_relabel.pdb 10000 10
rm -f epsilon.xml  # remove the soft link we just made
cp final_eH_relabel.pdb postLabel.pdb
phenix.reduce -q -trim -noflip final_eH_relabel.pdb | grep -v USER > final_relabel.pdb

# old way: jump out to different script: ${SCRIPTSDIR}/final_refine.sh final_relabel.pdb $mtz

resolrange=`phenix.mtz.dump $mtz | grep "Resolution range:"`
resol=`echo $resrange | awk '{print $4}' | cut -c 1-5`
resol1000=`echo $resol | awk '{tot=$1*1000}{print tot }'`

if [ $resol1000 -lt 1450 ]; then
    adp='adp.individual.anisotropic="not (water or element H)"' 
else
    adp='adp.individual.isotropic=all'
fi

#phenix.refine $mtz final_relabel.pdb refinement.refine.strategy=individual_sites_real_space \
phenix.refine $mtz final_relabel.pdb refinement.refine.strategy=occupancies+individual_sites+individual_adp \
    refinement.main.number_of_macro_cycles=3 output.prefix=final_relabel output.serial=2 --overwrite \
    `cat phenix_flags.txt` write_maps=false write_def_file=false write_eff_file=false write_geo_file=false

greptmp=`grep COORD final_relabel_002.pdb`
dist=`echo $greptmp | awk '{print $8*2.0}'`
chns=`awk '{printf "%s ", $1}' chains.txt`

cp final_relabel_002.pdb final_relabel_002_cluster.pdb
for c in $chns; do
    ln -s $c/fragids.txt
    $QFITDIR/TRUNCATE/cluster final_relabel_002_cluster.pdb final_relabel_002_cluster.pdb $c $dist > cluster_$c.log
    #rm -f fragids.txt
done

cp final_relabel_002_cluster.pdb final_relabel_002_cluster_unit.pdb
for c in $chns; do
    $QFITDIR/TRUNCATE/occ2unit final_relabel_002_cluster_unit.pdb final_relabel_002_cluster_unit.pdb $c
done

phenix.refine $mtz final_relabel_002_cluster_unit.pdb \
    refinement.refine.strategy=occupancies+individual_sites+individual_adp \
    refinement.main.number_of_macro_cycles=5 "$adp" --overwrite `cat phenix_flags.txt` \
    write_maps=false write_def_file=false write_eff_file=false write_geo_file=false

awk '{if(substr($0,1,6)=="CRYST1" || (substr($0,1,4)=="ATOM" || substr($0,1,6)=="HETATM") && substr($0,57,4)+0>0.19) print $0}' \
    final_relabel_002_cluster_unit_refine_001.pdb > final_relabel_002_cluster_unit_refine_001_occ.pdb

for c in $chns; do
    $QFITDIR/TRUNCATE/occ2unit final_relabel_002_cluster_unit_refine_001_occ.pdb final_relabel_002_cluster_unit_refine_001_occ.pdb $c
done

phenix.refine $mtz final_relabel_002_cluster_unit_refine_001_occ.pdb \
    refinement.refine.strategy=occupancies+individual_sites+individual_adp \
    refinement.main.number_of_macro_cycles=5 "$adp" --overwrite `cat phenix_flags.txt` \
    write_maps=false write_def_file=false write_eff_file=false write_geo_file=false

grep "Final R" final_relabel_002_cluster_unit_refine_001_occ_refine_001.log | tail -n 1 >> qFit.log

cp final_relabel_002_cluster_unit_refine_001_occ_refine_001.pdb qFit.pdb
cp final_relabel_002_cluster_unit_refine_001_occ_refine_001.mtz qFit.mtz

#rm -f elbow* minocc.txt *.params phenix_flags.txt #final* *_refine*
