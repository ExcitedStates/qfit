#!/bin/bash -f

pdb=$1
chain=$2

# The "reDokkum_frag" executable (in bin/ with the "qFit" executable)
# should have already made the "##_##" directories

#for frag in [1-9]*_[1-9]*; do
#    if [ -d $frag ]; then
#        echo "working on fragment $frag in chain $chain"
#        awk 'BEGIN{OFS=""} {if(($0~"ATOM  ")||($0~"HETATM")) print substr($0,0,21),'$chain',substr($0,23); else print $0}' $frag/frag.pdb > ${frag}_frag.pdb
#    fi
#    i=$((i+1))
#done

#bn=`basename $pdb .pdb`

for frag in [1-9]*_[1-9]*; do
    first=${frag%_*}
    last=${frag#*_}
    if [ -d $frag ]; then
        
        # Old way with CombineFrags:
        #/netapp/home/dkeedy/qFit/CombineFrags/bin/CombineFrags \
        #$QFITDIR/CombineFrags/bin/CombineFrags \
        #    $pdb ${frag}_frag.pdb ${frag}_frag.pdb
        #cp ${bn}_${frag}_frag.pdb $bn.pdb
        
        # New way with just awk:
        if [ -e $frag/frag.pdb ]; then
            # Use the frag.pdb file made by reDokkum_frag
            awk '{if(substr($0,22,1)=="'$chain'"&&substr($0,23,4)+0<'$first') print $0}' $pdb > tmp.pdb
            awk -v chain="$chain" 'BEGIN{OFS=""} {if(($0~"ATOM  ")||($0~"HETATM")) print substr($0,0,21),chain,substr($0,23)}' $frag/frag.pdb >> tmp.pdb
            awk '{if(substr($0,22,1)=="'$chain'"&&substr($0,23,4)+0>'$last') print $0}' $pdb >> tmp.pdb
            mv tmp.pdb $pdb
        fi
        # else reDokkum_frag failed for some reason, so we keep what's in the input for this residue range instead
    fi
  i=$((i+1))
done
