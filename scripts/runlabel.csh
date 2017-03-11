#!/bin/tcsh -f

source /home/sw/rhel5/x86_64/phenix/phenix-1.8-1069/phenix_env
ln -s /data/vdbedem/qFitRun/DOKKUMDIR/epsilon.xml
phenix.reduce $1 | grep -v USER > finalH.pdb
/data/vdbedem/qFitRun/DOKKUMDIR/Label finalH.pdb finalH_relabel.pdb 10000 10
phenix.reduce -Trim -NOFLIP finalH_relabel.pdb | grep -v USER > $2
rm -f epsilon.xml
