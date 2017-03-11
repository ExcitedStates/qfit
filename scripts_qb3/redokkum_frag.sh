#!/bin/bash

# This exists as a separate script entirely so we can call it with different mem_free requests.

pdb_scale=$1
mtz=$2
chain=$3
fragstart=$4
fragend=$5
pdb=$6
log=$7

$QFITDIR/qFit/bin/reDokkum_frag \
    $pdb_scale \
    $mtz \
    $chain \
    $fragstart \
    $fragend \
    $pdb \
    > reDokkum_frag.log
