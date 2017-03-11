#!/bin/tcsh

source /home/sw/rhel5/x86_64/phenix/phenix-1.8-1069/phenix_env

set pdb = $1
set mtz = $2
set F = $3
set SF = $4
set elbow=$5

# DETERMINE RESOLUTION AND (AN)ISOTROPIC REFINEMENT
#set resrange = `phenix.mtz.dump $mtz | grep "Resolution range:"`
#set res = `echo $resrange[4] | cut -c 1-5`
#set res1000 = `echo $res | awk '{tot = $1*1000}{print tot }'`

set paramfile = makemtz.params

echo "refinement {" > $paramfile
echo "  electron_density_maps {" >> $paramfile
echo "    map_coefficients {" >> $paramfile
echo "      mtz_label_amplitudes = FWT" >> $paramfile
echo "      mtz_label_phases = PHWT" >> $paramfile
echo "      map_type = 2mFo-DFc" >> $paramfile
echo "    }" >> $paramfile
echo "  }" >> $paramfile
echo "}" >> $paramfile

#if ( $res1000 < 1450 ) then
#  phenix.refine $mtz $pdb adp.individual.anisotropic="not (water or element H)" adp.individual.isotropic="water or element H" refinement.input.xray_data.labels=$F,$SF $paramfile --overwrite
#else
phenix.refine $mtz $pdb adp.individual.isotropic=all refinement.input.xray_data.labels=$F,$SF $elbow $paramfile --overwrite
 #endif

