#!/bin/tcsh

source /home/sw/rhel5/x86_64/phenix/phenix-1.6.4-486/phenix_env

phenix.refine ../../../L2_L123_eden-unique_no_phase.mtz ../1_001_cluster_unit.pdb ../3.params --overwrite refinement.main.use_experimental_phases=False refinement.input.xray_data.labels=F_l2,SIGF_l2 optimize_wxc=true optimize_wxu=true

../preprune.csh ../../../L2_L123_eden-unique_no_phase.mtz 1_001_cluster_unit_refine_001.pdb refinement.input.xray_data.labels=F_l2,SIGF_l2 > preprune.log
source /home/sdcsoftware/linux/phenix-1.3b/phenix_env
phenix.get_cc_mtz_pdb 1_001_cluster_unit_refine_001_map_coeffs.mtz 1_001_cluster_unit_refine_001_preprune.pdb "FP=2FOFCWT PHIB=PH2FOFCWT"
/home/vdbedem/DEV/ModelBuildTool/MBT/LPEnsemble/CODE/TRUNCATE/prune 1_001_cluster_unit_refine_001.pdb 1_001_cluster_unit_refine_001_prune.pdb A > prune.log

phenix.refine ../../../L2_L123_eden-unique_no_phase.mtz 1_001_cluster_unit_refine_001_prune.pdb ../3.params --overwrite refinement.main.use_experimental_phases=False refinement.input.xray_data.labels=F_l2,SIGF_l2 output.prefix=1_001_cluster_unit_refine_001_prune_wxoptimize optimize_wxc=true optimize_wxu=true

rm *.map
