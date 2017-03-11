#!/bin/bash

pdb=$1
mtz=$2
chain=$3
first=$4
last=$5

log=$PWD/$SGE_TASK_ID.log

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
# make sure this residue number (SGE_TASK_ID) is in the whitelist.
this_res_in_target_list=false
for resnum in ${residues[@]}; do
    if [ $resnum -eq $SGE_TASK_ID ]; then
        echo "$chain $SGE_TASK_ID IS in final target residue list, so running residue.sh" #>> $log
        this_res_in_target_list=true
        break
    fi
done
if ! $this_res_in_target_list ; then
    echo "$chain $SGE_TASK_ID NOT in final target residue list, so skipping residue.sh" #>> $log
    sleep 5  # give a little time for that ^ echo statement to reach the log file (?)
    exit
fi

################################  PREP FOR AND RUN QFIT  ################################

if [ ! -d $SGE_TASK_ID ]; then
    mkdir $SGE_TASK_ID
fi

cp $pdb $SGE_TASK_ID/in.pdb
cp $mtz $SGE_TASK_ID/in.hkl
cp ../minocc.txt $SGE_TASK_ID
cp ../phenix_flags.txt $SGE_TASK_ID

cd $SGE_TASK_ID

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
    adp="chain $chain and resid $SGE_TASK_ID"
fi

paramfile=${chain}${SGE_TASK_ID}adp.params

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
prevresid=$((SGE_TASK_ID-1)) # for a distal dummy selection, so pdbtools works for Gly/Ala
nextresid=$((SGE_TASK_ID+1)) # not used!
cp in.pdb in_${chain}${SGE_TASK_ID}.pdb
phenix.pdbtools \
    in_${chain}${SGE_TASK_ID}.pdb \
    modify.selection="chain $chain and ( (resseq $SGE_TASK_ID and not ( name n or name ca or name c or name o or name cb ) ) or ( resseq $prevresid and ( name n ) ) )" \
    modify.occupancies.set=0 \
    stop_for_unknowns=False
cp in_${chain}${SGE_TASK_ID}.pdb_modified.pdb in_${chain}${SGE_TASK_ID}_modified.pdb
phenix.refine in.hkl in_${chain}${SGE_TASK_ID}_modified.pdb --overwrite $paramfile \
    write_def_file=false write_eff_file=false write_geo_file=false `cat phenix_flags.txt`
$QFITDIR/TRUNCATE/truncateSC \
    in_${chain}${SGE_TASK_ID}_modified_refine_001.pdb \
    in_${chain}${SGE_TASK_ID}_modified_refine_001_trunc.pdb \
    $chain \
    $SGE_TASK_ID

# Decide whether to try flipping this peptide or not based on the secondary structure
# NB: The CalcSecStructure method in MMDB fails if there are any alt confs!
phenix.pdbtools \
    remove_alt_confs=True \
    stop_for_unknowns=False \
    in.pdb \
    output.file_name=in_1conf.pdb
flip=`$QFITDIR/TRUNCATE/flipornot in_1conf.pdb $chain $SGE_TASK_ID`
# Provide a placeholder "do not flip" string
# so that we avoid trying to pass an empty string to another bash script,
# which would mess up the $1, $2, ... argument interpretation
if [[ $flip == "" ]]; then
    flip="DO-NOT-FLIP-PLACEHOLDER-STRING"
fi

uname -a > ${chain}${SGE_TASK_ID}resolve.log
echo `pwd` >> ${chain}${SGE_TASK_ID}resolve.log

# $QFITDIR/qFit/bin/qFit \
#     $flip \
#     in_${chain}${SGE_TASK_ID}_modified_refine_001_trunc.pdb \
#     in_${chain}${SGE_TASK_ID}_modified_refine_001.mtz \
#     $chain \
#     $SGE_TASK_ID \
#     >> ${chain}${SGE_TASK_ID}qFit.log

basepdb=`basename $pdb .pdb`

tmpstr=`mktemp`
basestr=`basename $tmpstr`
jobbase=`echo $basestr | cut -c 1-10`

memfree="1G" # we're fine with nodes with less memory
if [[ $flip == "--flip-peptide" ]]; then
    memfree="5G" # request a node with > 4 GB
fi

#cp * $CURRLOC
cp in_${chain}${SGE_TASK_ID}_modified_refine_001_trunc.pdb \
   in_${chain}${SGE_TASK_ID}_modified_refine_001.mtz \
   $CURRLOC
cd $CURRLOC

## If we got the memory requirement ~right, 24 hours should be PLENTY of time for this
#    -l h_rt=24:0:0 \

# Let's try the max time allocation just to be safe
qsub \
    -S /bin/bash \
    -l arch=linux-x64 \
    -l netapp=1G \
    -l scratch=1G \
    -l mem_free=${memfree} \
    -l h_rt=168:0:0 \
    -cwd \
    -o qfit_${chain}${SGE_TASK_ID}.out \
    -e qfit_${chain}${SGE_TASK_ID}.err \
    -j n \
    -N q_${basepdb}_${chain}${SGE_TASK_ID}_${jobbase} \
    -V \
    -sync y \
    $SCRIPTSDIR/qfit.sh \
        $flip \
        in_${chain}${SGE_TASK_ID}_modified_refine_001_trunc.pdb \
        in_${chain}${SGE_TASK_ID}_modified_refine_001.mtz \
        $chain \
        $SGE_TASK_ID \
        ${chain}${SGE_TASK_ID}qFit.log

# wait 'til that's done ("-sync y")...
# We have to do that ^ because the finish_chain job is waiting for this residue job array
# to finish -- but it isn't *efffectively* complete until the qfit jobs finish...

#rm in.hkl in.pdb *.map

#cp * $CURRLOC
#cp Res*.pdb after*.pdb *.log $CURRLOC

#cd $CURRLOC  # back from TMPLOC to residue (SGE_TASK_ID)

cd ..  # back from residue (SGE_TASK_ID) to chain 

echo "done with residue.sh for residue $chain $SGE_TASK_ID at `date`" #>> $log

##########################  CHECK STATUS OF DAISY CHAIN  #########################
#
#all_residues_done=true
#for resnum in ${residues[@]}; do
#    # qFit makes Res##.pdb for each residue (including Gly) when it's done, 
#    # but I'm not totally convinced that successfully happens in all edge cases 
#    # (e.g. weird occupancy combinations), so let's just check that residue.sh 
#    # has *reported* as being complete for all other residues.  That will be false 
#    # for residue.sh jobs that are still in progress or failed due to some bug.
#    done_1_or_0=`grep "done with residue.sh" $resnum.log | wc | awk '{print $1}'`
#    if [ $done_1_or_0 -eq 0 ]; then
#        echo "at least one residue ($resnum) not done yet" >> $log
#        all_residues_done=false
#        break
#    fi
#done
#
#if $all_residues_done ; then
#    biglog=$first-$last.log
#    for (( resnum=$first; resnum<=$last; resnum++ )); do
#        cat $resnum.log >> $biglog
#        rm $resnum.log
#    done
#    echo "all residues done! ($SGE_TASK_ID was the last one)" >> $biglog
#
#    tmpstr=`mktemp`
#    basestr=`basename $tmpstr`
#    jobbase=`echo $basestr | cut -c 1-10`
#
#    basepdb=`basename $pdb .pdb`
#
#    qsub \
#        -S /bin/bash \
#        -l arch=linux-x64,netapp=1G,scratch=1G,mem_free=4G,h_rt=168:0:0 \
#        -cwd -o finish_chain_${chain}.out -e finish_chain_${chain}.err -j n \
#        -N fc_${basepdb}${chain}_${jobbase} \
#        -V \
#        $SCRIPTSDIR/finish_chain.sh $pdb $mtz $chain $first $last
#    echo "qsubbed finish_chain.sh for chain $chain at `date`" >> $biglog
#fi
#