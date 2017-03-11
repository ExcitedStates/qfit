#!/bin/tcsh

../DOKKUM/get_Res.csh

grep -v END ALL.pdb > ALL_.pdb

pdbset xyzin ALL_.pdb xyzout ALL_.pdb << eopdb > pdbset.log
exclude hydrogens
chain A
end
eopdb

/home/vdbedem/DEV/ModelBuildTool/MBT/LPEnsemble/CODE/TRUNCATE/idmulti ALL_.pdb ALL_idmulti.pdb A > fragids.txt


