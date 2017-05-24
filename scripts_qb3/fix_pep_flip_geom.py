from __future__ import division
import iotbx.pdb
import sys
import os
import itertools

def get_prev_rsd_flip_occs(rg):
  
  co_angles = []
  o_occs = {}

  altlocs = set()
  for ag in rg.atom_groups():
    altlocs.add(ag.altloc)

  for altloc_pair in itertools.combinations(altlocs, 2):
    if " " in altloc_pair or "" in altloc_pair:
      continue
    
    c1_atom = None
    o1_atom = None
    c2_atom = None
    o2_atom = None
    for ag in rg.atom_groups():
      if ag.altloc == altloc_pair[0]:
        for atom in ag.atoms():
          if atom.name == " C  ":
            c1_atom = atom
          elif atom.name == " O  ":
            o1_atom = atom
            o_occs[ag.altloc] = atom.occ
    for ag in rg.atom_groups():
      if ag.altloc == altloc_pair[1]:
        for atom in ag.atoms():
          if atom.name == " C  ":
            c2_atom = atom
          elif atom.name == " O  ":
            o2_atom = atom
            o_occs[ag.altloc] = atom.occ
    
    if not (c1_atom and o1_atom and c2_atom and o2_atom):
      return None

    # Translate second C=O so carbons are aligned
    new_o2_xyz = []
    for i in range(3):
      diff = c2_atom.xyz[i] - c1_atom.xyz[i]
      new_o2_xyz.append(o2_atom.xyz[i] - diff)
    new_o2_atom = o2_atom.detached_copy()
    new_o2_atom.xyz = new_o2_xyz

    co_angle = c1_atom.angle(atom_1=o1_atom, atom_3=new_o2_atom, deg=True)
    # ^ c1_atom is "atom_2", so to speak, in the angle
    co_angles.append(co_angle)

  if len(co_angles) >= 1 and max(co_angles) > 90:
    return o_occs
  return None

if __name__ == "__main__":

  if len(sys.argv) != 2:
    print >> sys.stderr, \
      'Usage: python %s in.pdb' % os.path.basename(__file__)
    sys.exit(1)

  file_name = sys.argv[1]
  pdb_obj = iotbx.pdb.hierarchy.input(file_name=file_name)
  for model in pdb_obj.hierarchy.models():
    for chain in model.chains():
      
      prev_rsd_flip_occs = None

      for rg in chain.residue_groups():
        
        if prev_rsd_flip_occs:
          
          # See if we need to split the N
          N_altlocs = set()
          #N_atom = None
          for ag in rg.atom_groups():
            for atom in ag.atoms():
              if atom.name.strip() == "N":
                N_altlocs.add(ag.altloc)
                #N_atom = atom
          
          #if len(N_altlocs) == 1:
          if len(N_altlocs) != len(prev_rsd_flip_occs):
            # We *DO* need to split the N
            print chain.id, rg.resseq, rg.atom_groups()[0].resname, \
              "needs to be split to fix the i to i-1 peptide geometry"
            
            new_ags_to_append = []
            
            for altloc in prev_rsd_flip_occs:
              print "altloc", altloc, prev_rsd_flip_occs[altloc]
              
              matching_ag = None
              matching_N_atom = None
              for ag in rg.atom_groups():
                if ag.altloc == altloc:
                  matching_ag = ag
                  for atom in ag.atoms():
                    if atom.name.strip() == "N":
                      matching_N_atom = atom

              # OPTION 1: atom_group exists and already has N
              if matching_ag and matching_N_atom:
                #print 'OPTION 1'
                matching_N_atom.occ = prev_rsd_flip_occs[altloc]
                # This ^ makes the N's occ match that of the preceding CO,
                # but its occ could now not match its own CA's occ!

              # OPTION 2: atom_group exists, but N does not 
              # (e.g. atom_group is only for sidechain)
              elif matching_ag and (not matching_N_atom):
                #print 'OPTION 2'
                new_ag = matching_ag.detached_copy()
                # Possible TODO: change coordinates?
                new_ag.occupancy = prev_rsd_flip_occs[altloc]
                #print new_ag.occupancy
                for atom in new_ag.atoms():
                  if atom.name.strip() != "N":
                    new_ag.remove_atom(atom)
                  else:
                     atom.occ = prev_rsd_flip_occs[altloc]
                rg.merge_atom_groups(matching_ag, new_ag)

              # OPTION 3: atom_group does not exist at all
              else:
                #print 'OPTION 3'
                new_ag = rg.atom_groups()[0].detached_copy()
                # Possible TODO: pick altloc that is closest / has best geometry
                # to use for template atom_group, instead of just default index 0?
                # Possible TODO: change coordinates?
                new_ag.altloc = altloc
                new_ag.occupancy = prev_rsd_flip_occs[altloc]
                #print new_ag.occupancy
                for atom in new_ag.atoms():
                  if atom.name.strip() != "N":
                    new_ag.remove_atom(atom)
                  else:
                     atom.occ = prev_rsd_flip_occs[altloc]
                new_ags_to_append.append(new_ag)

            for new_ag in new_ags_to_append:
              rg.append_atom_group(new_ag)
            #rg.atom_groups()[0].remove_atom(N_atom)

        # Regardless of what we did for this residue,
        # plan ahead for the next residue
        prev_rsd_flip_occs = get_prev_rsd_flip_occs(rg)
        if prev_rsd_flip_occs:
          print chain.id, rg.resseq, rg.atom_groups()[0].resname, \
            "has a peptide flip (CO rotates by >90 degrees)"

  suffix = "_fixedPepFlipGeom"
  output_pdb = os.path.basename(file_name).split(".pdb")[0]+suffix+".pdb"
  pdb_obj.hierarchy.write_pdb_file(file_name=output_pdb,
    crystal_symmetry=pdb_obj.input.crystal_symmetry(), append_end=True)
  print "Wrote", output_pdb
