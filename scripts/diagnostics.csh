#!/bin/tcsh

foreach h ( sdc35 sdc36 sdc37 sdc38 sdc39 sdc40 sdc41 sdc42 sdc43 sdc44 sdc45 sdc46 sdc47 sdc48 sdc49 sdc50 sdc51 sdc52 sdc53 sdc54 sdc55 sdc56 sdc57 sdc58 sdc59 sdc60 sdc61 sdc62 sdc63 sdc64 sdc89 sdc90 sdc91 sdc92 sdc93 sdc94 sdc95 sdc96 sdc97 sdc98 )
  bsub -qsdcq -m "$h" -oout.%J "ls /scratch/vdbedem/"
end
