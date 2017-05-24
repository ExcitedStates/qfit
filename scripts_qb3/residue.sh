#!/bin/bash

pdb=$1
mtz=$2
chain=$3
first=$4
last=$5
flip_or_noflip=$6
res_nums_commas=$7

res_nums_spaces=`echo $res_nums_commas | tr , " "`
tasks=(0 $res_nums_spaces)
RES_NUM=${tasks[$SGE_TASK_ID]}

log=$PWD/$RES_NUM.log

######################  DECIDE IF WE'LL RUN QFIT FOR THIS RESIDUE  ######################

# First, build array of target residue numbers in this chain
# we're supposed to run qFit on (will get reused later).
# Yes, this means we have to do again for every residue job,
# but I couldn't figure out a nice way to make bash pass arrays as arguments,
# and this should add essentially negligibly to each residue job's runtime.
# (Reminder: we're in chain directory, and residues.txt (if it exists)
# should be in main work directory, i.e. one level up.)
residues=()
if [ -e ../residues.txt ]; then
# format: "A 1"
    for resnum in `grep $chain ../residues.txt | awk '{print $2}'`; do
        if [[ $resnum -lt $first || $resnum -gt $last ]]; then
            echo "$chain $resnum from residues.txt OUTSIDE $first-$last range, so skipping" #>> $log
        continue
        else
            echo "$chain $resnum from residues.txt INSIDE first-last range, so running qFit" #>> $log
            residues+=($resnum)
        fi
    done
    echo "final CUSTOM target residue list: ${residues[@]}" #>> $log
else
    for (( resnum=$first; resnum<=$last; resnum++ )); do
        echo "residues.txt not provided, so continuing with $chain $resnum" #>> $log
        residues+=($resnum)
    done
    echo "final NON-custom target residue list: ${residues[@]}" #>> $log
fi

# Second, do the check for this particular residue:
# make sure this residue number (RES_NUM) is in the whitelist.
this_res_in_target_list=false
for resnum in ${residues[@]}; do
    if [ $resnum -eq $RES_NUM ]; then
        echo "$chain $RES_NUM IS in final target residue list, so running residue.sh" #>> $log
        this_res_in_target_list=true
        break
    fi
done
if ! $this_res_in_target_list ; then
    echo "$chain $RES_NUM NOT in final target residue list, so skipping residue.sh" #>> $log
    sleep 5  # give a little time for that ^ echo statement to reach the log file (?)
    exit
fi

################################  PREP FOR AND RUN QFIT  ################################

if [ ! -d $RES_NUM ]; then
    mkdir $RES_NUM
fi

cp $pdb $RES_NUM/in.pdb
cp $mtz $RES_NUM/in.hkl
cp ../minocc.txt $RES_NUM
cp ../phenix_flags.txt $RES_NUM

cd $RES_NUM

CURRLOC=$PWD # in a ##/ residue directory

TMPDIR=/scratch/$user
if [ ! -d $TMPDIR ]; then
    mkdir -p $TMPDIR
fi
TMPLOC=`mktemp -d -p $TMPDIR` || exit 1

cp in.pdb $TMPLOC
cp in.hkl $TMPLOC
cp minocc.txt $TMPLOC
cp phenix_flags.txt $TMPLOC

cd $TMPLOC

phenix.pdbtools in.pdb stop_for_unknowns=False \
    modify.rename_chain_id.old_id=$chain \
    modify.rename_chain_id.new_id=$chain \
    > pdbtools_rename.log
mv in.pdb_modified.pdb in.pdb

resolrange=`phenix.mtz.dump in.hkl | grep "Resolution range:"`
resol=`echo $resolrange | awk '{print $4}' | cut -c 1-5`
resol1000=`echo $resol | awk '{tot=$1*1000}{print tot }'`

if [ $resol1000 -lt 1450 ]; then
    adp="all"
else
    adp="chain $chain and resid $RES_NUM"
fi

paramfile=${chain}${RES_NUM}adp.params

echo "refinement {" > $paramfile
echo "  electron_density_maps {" >> $paramfile
echo "    map_coefficients {" >> $paramfile
echo "      mtz_label_amplitudes=FWT" >> $paramfile
echo "      mtz_label_phases=PHWT" >> $paramfile
echo "      map_type=2mFo-DFc" >> $paramfile
echo "    }" >> $paramfile
echo "  }" >> $paramfile
echo "  refine {" >> $paramfile
echo "    strategy=*individual_sites *individual_adp" >> $paramfile
echo "    adp {" >> $paramfile
echo "      individual {" >> $paramfile
echo "        anisotropic=$adp" >> $paramfile
echo "      }" >> $paramfile
echo "    }" >> $paramfile
echo "  }" >> $paramfile
echo "}" >> $paramfile

# Get omit map and starting model
prevresid=$((RES_NUM-1)) # for a distal dummy selection, so pdbtools works for Gly/Ala
nextresid=$((RES_NUM+1)) # not used!
cp in.pdb in_${chain}${RES_NUM}.pdb
phenix.pdbtools \
    in_${chain}${RES_NUM}.pdb \
    modify.selection="chain $chain and ( (resseq $RES_NUM and not ( name n or name ca or name c or name o or name cb ) ) or ( resseq $prevresid and ( name n ) ) )" \
    modify.occupancies.set=0 \
    stop_for_unknowns=False
cp in_${chain}${RES_NUM}.pdb_modified.pdb in_${chain}${RES_NUM}_modified.pdb
phenix.refine in.hkl in_${chain}${RES_NUM}_modified.pdb --overwrite $paramfile \
    write_def_file=false write_eff_file=false write_geo_file=false `cat phenix_flags.txt`
$QFITDIR/TRUNCATE/truncateSC \
    in_${chain}${RES_NUM}_modified_refine_001.pdb \
    in_${chain}${RES_NUM}_modified_refine_001_trunc.pdb \
    $chain \
    $RES_NUM

flip="DON'T KNOW YET"
if [ $flip_or_noflip == "flip" ]; then
    flip="--flip-peptide"
elif [ $flip_or_noflip == "noflip" ]; then
    flip=""
else
    echo "WARNING: flip_or_noflip variable for residue $RES_NUM set to something weird (not 'flip' or 'noflip'): $flip_or_noflip"
    echo "Not flipping, just to be safe..."
    flip=""
fi

uname -a > ${chain}${RES_NUM}resolve.log
echo `pwd` >> ${chain}${RES_NUM}resolve.log

echo "Memory usage before running qFit:"
qstat -f -j $JOB_ID | grep mem | grep "${RES_NUM}:" >> ${chain}${RES_NUM}qFit.log

$QFITDIR/qFit/bin/qFit \
    $flip \
    in_${chain}${RES_NUM}_modified_refine_001_trunc.pdb \
    in_${chain}${RES_NUM}_modified_refine_001.mtz \
    $chain \
    $RES_NUM \
    >> ${chain}${RES_NUM}qFit.log

echo "Memory usage after running qFit:"
qstat -f -j $JOB_ID | grep mem | grep "${RES_NUM}:" >> ${chain}${RES_NUM}qFit.log

rm *.map *.geo
rm in.hkl in.pdb

cp Res*.pdb *.out *.err *.log $CURRLOC

cd $CURRLOC  # back from TMPLOC to residue (RES_NUM)
cd ..  # back from residue (RES_NUM) to chain 

echo "done with residue.sh for residue $chain $RES_NUM at `date`" #>> $log

