#!/bin/bash

# This exists as a separate script entirely so we can call it with different mem_free requests.

flip=$1
inpdb=$2
inmtz=$3
chain=$4
resnum=$5
log=$6

# Provide a placeholder "do not flip" string
# so that we avoid trying to pass an empty string to another bash script,
# which would mess up the $1, $2, ... argument interpretation
if [[ $flip == "DO-NOT-FLIP-PLACEHOLDER-STRING" ]]; then
    flip=""
fi

CURRLOC=$PWD # in a ##/ residue directory

TMPDIR=/scratch/$user
if [ ! -d $TMPDIR ]; then
    mkdir -p $TMPDIR
fi
TMPLOC=`mktemp -d -p $TMPDIR` || exit 1

cp $inpdb $TMPLOC
cp $inmtz $TMPLOC
cp minocc.txt $TMPLOC

cd $TMPLOC

$QFITDIR/qFit/bin/qFit \
    $flip \
    $inpdb \
    $inmtz \
    $chain \
    $resnum \
    >> $log

#cp * $CURRLOC
cp Res*.pdb qfit*.out qfit*.err *qFit.log $CURRLOC
cd $CURRLOC
