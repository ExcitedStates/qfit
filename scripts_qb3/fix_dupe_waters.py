# This script trys to avoid the following error message from Label:
# "Duplicate sequence number and insertion code."
# It only addresses waters for now since that's where I've seen
# the problem thus far and since numbering for waters is more arbitrary
# than for other molecules -- but it could easily be made more general. 
# Actually, as of 11/25/13, I can't remember which cases actually had
# this problem -- but I do remember such cases existing...

import sys

def parse_line(line):
  altconf = line[16:17]
  restype = line[17:20]
  chain = line[21:22]
  resseq = int(line[22:26].strip())
  icode = line[26:27]
  return chain, resseq, icode, restype, altconf

# Figure out the max residue number per chain so we can know 
# where to start incrementing from in the next step
max_resseqs = {} # chain --> int
f = open(sys.argv[1], "r")
lines = f.readlines()
f.close()
for line in lines:
  if line.startswith("ATOM  ") or line.startswith("HETATM"):
    chain, resseq, icode, restype, altconf = parse_line(line)
    if chain not in max_resseqs:
      max_resseqs[chain] = resseq
    elif resseq > max_resseqs[chain]:
      max_resseqs[chain] = resseq

# Find combinations of residue number + insertion code that are duplicated
# within a chain, and immediately renumber them
chain_resseq_icodes = set()
prev_line_resseq = None
for line in lines:
  if line.startswith("ATOM  ") \
  or line.startswith("HETATM") \
  or line.startswith("ANISOU"):
    chain, resseq, icode, restype, altconf = parse_line(line)
    cri = (chain, resseq, icode)
    if cri in chain_resseq_icodes and altconf == ' ' and restype == 'HOH':
      # Duplicate water!
      if line.startswith("ATOM  ") or line.startswith("HETATM"):
        new_max_resseq = max_resseqs[chain] + 1
        max_resseqs[chain] = new_max_resseq
        new_resseq = new_max_resseq
        prev_line_resseq = new_max_resseq
      elif line.startswith("ANISOU"):
        # Assume previous line was ATOM or HETATM for same residue...
        #print 'HOH %d' % prev_line_resseq
        new_resseq = prev_line_resseq
      new_line = line[0:22] + str(new_resseq).rjust(4) + line[26:]
      print new_line ,
    else:
      # Unique (so far) ATOM/HETATM/ANISOU, or alt conf
      print line ,
      prev_line_resseq = resseq
    if line.startswith("ATOM  ") or line.startswith("HETATM"):
      chain_resseq_icodes.add(cri)
  else:
    # Header or something
    print line ,
