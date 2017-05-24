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
          N_atom = None
          for ag in rg.atom_groups():
            for atom in ag.atoms():
              if atom.name.strip() == "N":
                N_altlocs.add(ag.altloc)
                N_atom = atom
          
          if len(N_altlocs) == 1:
            # We *DO* need to split the N
            print chain.id, rg.resseq, rg.atom_groups()[0].resname, \
              "needs to be split to fix the i to i-1 peptide geometry"
            new_ags = []
            for altloc in prev_rsd_flip_occs:
              #new_atom = N_atom.detached_copy()
              #new_atom.occ = prev_rsd_flip_occs[altloc]
              # possible TODO: change coordinates?
              new_ag = rg.atom_groups()[0].detached_copy()
              new_ag.altloc = altloc
              new_ag.occupancy = prev_rsd_flip_occs[altloc]
              print new_ag.occupancy
              for atom in new_ag.atoms():
                if atom.name.strip() != "N":
                  new_ag.remove_atom(atom)
                else:
                   atom.occ = prev_rsd_flip_occs[altloc]
              new_ags.append(new_ag)
            for new_ag in new_ags:
              rg.append_atom_group(new_ag)
            rg.atom_groups()[0].remove_atom(N_atom)
            prev_rsd_flip_occs = None

        else:
          # Previous residue does not have peptide flip, so we 
          # don't need to worry about splitting this residue's N
          prev_rsd_flip_occs = get_prev_rsd_flip_occs(rg)
          if prev_rsd_flip_occs:
            print chain.id, rg.resseq, rg.atom_groups()[0].resname, \
              "has a peptide flip (CO rotates by >90 degrees)"
          else:
            prev_rsd_flip_occs = None

  suffix = "_fixedPepFlipGeom"
  output_pdb = os.path.basename(file_name).split(".pdb")[0]+suffix+".pdb"
  pdb_obj.hierarchy.write_pdb_file(file_name=output_pdb,
    crystal_symmetry=pdb_obj.input.crystal_symmetry(), append_end=True)
  print "Wrote", output_pdb
