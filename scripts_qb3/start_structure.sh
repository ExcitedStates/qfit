#!/bin/bash

export LD_LIBRARY_PATH=/netapp/home/dkeedy/qFit/mmdb-1.23.2.2/lib:/netapp/home/dkeedy/qFit/cplex/ILOG/CPLEX_Studio124/concert/lib/x86-64_sles10_4.1/static_pic:/netapp/home/dkeedy/qFit/cplex/ILOG/CPLEX_Studio124/cplex/lib/x86-64_sles10_4.1/static_pic:/netapp/home/dkeedy/qFit/clipper-2.1/lib
export ILOG_LICENSE_FILE=/netapp/home/dkeedy/qFit/cplex/ILOG/CPLEX_Studio124/licenses/LI_en.txt
export LOOPTK_RESOURCEDIR=/netapp/home/dkeedy/qFit/looptk-2.0.1/resources

export QFITDIR=/netapp/home/dkeedy/qFit/QFIT_2.1/branches/dev-2.1-qb3
export SCRIPTSDIR=$QFITDIR/scripts_qb3

workdir=""
pdb=""
mtz=""
sflab=""
rfreelab=""
genfree=""
residues=""

while getopts :w:p:m:l:f:g:o:r: opt; do
    case $opt in 
        w)
            workdir=$OPTARG ;;
        p)
            pdb=$OPTARG ;;
        m)
            mtz=$OPTARG ;;
        l)
            sflab=$OPTARG ;;
        f)
            rfreelab=$OPTARG ;;
        g)
            genrfree=$OPTARG ;;
        o)
            minocc=$OPTARG ;;
        r)
            residues=$OPTARG ;;
        \?)
            echo "Invalid option: -$OPTARG"; exit 1 ;;
        :)
            echo "Option -$OPTARG requires an argument"; exit 1 ;;
    esac
done

if [ "$workdir" = "" ]; then
    echo "Must supply work directory name with -w dirname"
    exit 1
fi
if [ "$pdb" = "" ]; then
    echo "Must supply PDB file with -p in.pdb"
    exit 1
fi
if [ "$mtz" = "" ]; then
    echo "Must supply MTZ file with -m in.mtz"
    exit 1
fi
if [ "$sflab" = "" ]; then
    echo "Must supply structure factor labels with (e.g.) -l \"FOBS,SIGFOBS\""
    exit 1
fi

if [ "$rfreelab" != "" ]; then
    echo "R-free labels were provided"
fi
if [ "$genfree" != "" ]; then
    echo "Request to have PHENIX generate R-free labels was provided"
fi
if [ "$residues" != "" ]; then
    echo "Target residue list file was provided, so renaming \"residues.txt\" and using"
    cp $residues $workdir/residues.txt
fi

echo $pdb $mtz > qFit.log
cat $minocc >> qFit.log
cat $residues >> qFit.log

source /netapp/home/dkeedy/phenix/phenix-1.9-1692/phenix_env.sh

mv $pdb $pdb.org

phenix.reduce -q -trim -noflip $pdb.org | grep -v USER > $pdb.trim
mv $pdb.trim $pdb

$SCRIPTSDIR/../TRUNCATE/renumberInsertions $pdb

python $SCRIPTSDIR/fix_dupe_waters.py $pdb > $pdb.fixhoh
mv $pdb.fixhoh $pdb

basepdb=`basename $pdb .pdb`

phenix.elbow --do-all $pdb
rm elbow*.pdb elbow*.pickle
elbow=""
elbowbasepdb=`echo $basepdb | tr . _`
if [ -r "elbow.${elbowbasepdb}_pdb.all.001.cif" ]; then
    elbow="$PWD/elbow.${elbowbasepdb}_pdb.all.001.cif"
fi
echo "restraints from elbow: $elbow"

echo $sflab > phenix_flags.txt
echo $elbow >> phenix_flags.txt
echo $rfreelab >> phenix_flags.txt
echo $genrfree >> phenix_flags.txt

$SCRIPTSDIR/../TRUNCATE/chains $pdb > chains.txt

#$SCRIPTSDIR/initial_refine.sh $pdb $mtz > initial_refine.log

paramfile=makemtz.params

echo "refinement {" > $paramfile
echo "  electron_density_maps {" >> $paramfile
echo "    map_coefficients {" >> $paramfile
echo "      mtz_label_amplitudes=FWT" >> $paramfile
echo "      mtz_label_phases=PHWT" >> $paramfile
echo "      map_type=2mFo-DFc" >> $paramfile
echo "    }" >> $paramfile
echo "  }" >> $paramfile
echo "}" >> $paramfile

phenix.refine $mtz $pdb adp.individual.isotropic=all \
    write_def_file=false write_eff_file=false write_geo_file=false \
    `cat phenix_flags.txt` $paramfile --overwrite

if [ ! -f ${basepdb}_refine_001.pdb ]; then
    grep -i sorry ${basepdb}_refine_001.log >> qFit.log
    sorries=`grep -i sorry ${basepdb}_refine_001.log | wc | awk '{print $1}'`
    if [ $((sorries)) -gt 0 ]; then
        echo "initial refinement failed -- exiting now!" >> qFit.log
        exit 1
    fi
fi

egrep '^ATOM|^HETATM' $pdb > noheader.pdb

TMPDIR=/scratch/$user
if [ ! -d $TMPDIR ]; then
    mkdir -p $TMPDIR
fi

tmpstr=`mktemp`
basestr=`basename $tmpstr`
jobbase=`echo $basestr | cut -c 1-10`
nchains=`wc chains.txt | awk '{print $1}'`

qsub \
    -S /bin/bash \
    -l arch=linux-x64,netapp=1G,scratch=1G,mem_free=1G,h_rt=0:29:0 \
    -cwd -o start_chain_1-$nchains.out -e start_chain_1-$nchains.err -j n \
    -N sc_${basepdb}_1-${nchains}_${jobbase} \
    -t 1-$nchains \
    -V \
    $SCRIPTSDIR/start_chain.sh $pdb $mtz
