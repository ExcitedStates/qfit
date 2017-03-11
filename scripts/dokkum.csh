#!/bin/tcsh

setenv DOKKUMDIR /data/vdbedem/qFitRun_DVLP2/DOKKUMDIR
setenv LOOPTK_RESOURCEDIR ${DOKKUMDIR}/LOOPTK/resources/
source /data/vdbedem/qFitRun_DVLP2/DOKKUMDIR/phenix-1.9-1692
#setenv LD_LIBRARY_PATH ${DOKKUMDIR}/coin-Cbc/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /home/vdbedem/DEV/ModelBuildTool/MBT/QFIT/ILOG/CPLEX_Studio_AcademicResearch122/cplex/lib/x86-64_sles10_4.1/static_pic/:/home/vdbedem/DEV/ModelBuildTool/MBT/QFIT/ILOG/CPLEX_Studio_AcademicResearch122/concert/lib/x86-64_sles10_4.1/static_pic/:$LD_LIBRARY_PATH
#source /home/sdcsoftware/linux/ccp4-6.0.2/include/ccp4.setup
setenv ILOG_LICENSE_FILE /home/vdbedem/DEV/ModelBuildTool/MBT/QFIT/ILOG/CPLEX_Studio_AcademicResearch122/licenses/access.ilm

setenv FLIP_PEPTIDE "--flip-peptide"
#set LSB_JOBINDEX = $4

if ( ! -d $LSB_JOBINDEX ) then
  mkdir $LSB_JOBINDEX
endif

set CURRLOC = $PWD/$LSB_JOBINDEX

if ( ! -d /scratch/$user ) then
  mkdir -p /scratch/$user
endif

set TMPDIR=`mktemp -d -p /scratch/$user` || exit 1

cp $1 ${TMPDIR}/in.pdb
cp $2 ${TMPDIR}/in.hkl
cp ../minocc.txt ${TMPDIR}

cd $TMPDIR

phenix.pdbtools in.pdb modify.rename_chain_id.old_id=$3 modify.rename_chain_id.new_id=$3 $5 > pdbtools_rename.log
#pdbset xyzin in.pdb xyzout in.pdb << eopdb > pdbset.log
#$chain $3 $3
#end
#eopdb
@ i = $LSB_JOBINDEX
@ i -= 1
set PRVS_RSD = $i
@ i += 2
set NEXT_RSD = $i


# DETERMINE RESOLUTION AND (AN)ISOTROPIC REFINEMENT
set resrange = `phenix.mtz.dump in.hkl | grep "Resolution range:"`
set res = `echo $resrange[4] | cut -c 1-5`
set res1000 = `echo $res | awk '{tot = $1*1000}{print tot }'`

if ( $res1000 < 1450 ) then
  set adp = "not (water or element H)"
else
  set adp = "chain $3 and resid $LSB_JOBINDEX"
endif

set paramfile = {$3}{$LSB_JOBINDEX}adp.params

echo "refinement {" > $paramfile
echo "  electron_density_maps {" >> $paramfile
echo "    map_coefficients {" >> $paramfile
echo "      mtz_label_amplitudes = FWT" >> $paramfile
echo "      mtz_label_phases = PHWT" >> $paramfile
echo "      map_type = 2mFo-DFc" >> $paramfile
echo "    }" >> $paramfile
echo "  }" >> $paramfile
echo "  refine {" >> $paramfile
echo "    strategy = *individual_sites *individual_adp" >> $paramfile
echo "    adp {" >> $paramfile
echo "      individual {" >> $paramfile
echo "        anisotropic = $adp" >> $paramfile
echo "      }" >> $paramfile
echo "    }" >> $paramfile
echo "  }" >> $paramfile
echo "}" >> $paramfile


#Set occupancies to 0 for original input to obtain omit map
cp in.pdb out_{$3}{$LSB_JOBINDEX}.pdb
phenix.pdbtools \
	modify.selection="chain $3 and ( resseq $LSB_JOBINDEX and not ( name n or name ca or name c or name o or name cb ) or ( resseq $PRVS_RSD and name n) )" \
	modify.occupancies.set=0 \
	stop_for_unknowns=False \
	out_{$3}{$LSB_JOBINDEX}.pdb \
	output.file_name=out_{$3}{$LSB_JOBINDEX}_modified.pdb

phenix.reduce out_{$3}{$LSB_JOBINDEX}_modified.pdb > out_{$3}{$LSB_JOBINDEX}_modified_H.pdb

phenix.elbow --do-all out_{$3}{$LSB_JOBINDEX}_modified_H.pdb
set elbowbasepdb = out_{$3}{$LSB_JOBINDEX}_modified_H
set elbow = ""
if ( -r elbow.{$elbowbasepdb}_pdb.all.001.cif ) then
   set elbow = "$cwd/elbow.{$elbowbasepdb}_pdb.all.001.cif"
endif
echo "elbow:" $elbow

phenix.refine in.hkl out_{$3}{$LSB_JOBINDEX}_modified_H.pdb --overwrite $paramfile $4 $elbow
$DOKKUMDIR/TRUNCATE/truncateSC  \
	out_{$3}{$LSB_JOBINDEX}_modified_H_refine_001.pdb \
	out_{$3}{$LSB_JOBINDEX}_modified_H_refine_001_trunc.pdb \
	$3 \
	$LSB_JOBINDEX
phenix.reduce -Trim out_{$3}{$LSB_JOBINDEX}_modified_H_refine_001_trunc.pdb > out_{$3}{$LSB_JOBINDEX}_modified_H_refine_001_trunc_trim.pdb

# Decide whether to try flipping this peptide or not based on the secondary structure
# NB: The CalcSecStructure method in MMDB fails if there are any alt confs!
phenix.pdbtools \
    remove_alt_confs=True \
    stop_for_unknowns=False \
    in.pdb \
    output.file_name=in_1conf.pdb
    
set flip = `$DOKKUMDIR/TRUNCATE/flip_or_not in_1conf.pdb $3 $LSB_JOBINDEX`

uname -a > {$3}{$LSB_JOBINDEX}qFit.log
echo `pwd` >> {$3}{$LSB_JOBINDEX}qFit.log

echo "$DOKKUMDIR/qFit ${flip} out_{$3}{$LSB_JOBINDEX}_H_refine_001.pdb out_{$3}{$LSB_JOBINDEX}_modified_H_refine_001.mtz $3 $LSB_JOBINDEX" >> {$3}{$LSB_JOBINDEX}qFit.log

$DOKKUMDIR/qFit \
	$flip \
	out_{$3}{$LSB_JOBINDEX}_modified_H_refine_001_trunc_trim.pdb \
	out_{$3}{$LSB_JOBINDEX}_modified_H_refine_001.mtz \
	$3 \
	$LSB_JOBINDEX \
	>> {$3}{$LSB_JOBINDEX}qFit.log

rm *.map *.geo elbow*

rm in.hkl in.pdb

cp * ${CURRLOC}
