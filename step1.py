import shutil

import numpy as np
import MDAnalysis as mda
import gromacs as gmx

"""Create initial conformation for simulation of protein in Octanol wall"""

# Copy needed files
shutil.copy("topol/prot.gro", "prot.gro")

# NOTE: files for development purposes
# shutil.copy("topol/box.gro", "conf.gro")
# shutil.copy("topol/box.top", "topol.top")

# Insert protein into slab
u1 = mda.Universe("conf.gro")
u2 = mda.Universe("prot.gro")
u = mda.Merge(u1.atoms, u2.atoms)

# Set new universe dimensions
u.dimensions = u1.dimensions
u.dimensions[0] = 200.0 # A

# Move antibody so that Fv region is on the edge of the octanol slab
protein = u.select_atoms("protein")
com = protein.select_atoms("resid 1-100").center_of_mass()
x, y, z = u1.dimensions[:3]
delta = np.array([x, y/2, z/2]) - com 
protein.atoms.positions += delta

# Remove octanols that overlap protein
new = u.select_atoms("not (resname OCOH and same resid as (around 2.5 protein))")
nmol = new.select_atoms("resname OCOH").n_residues
new.atoms.write("conf.gro")

# Alter topology to include new number of octanols and protein
with open("topol.top", "r") as f:
    lines = f.read().splitlines()
lines.pop()
lines.append(f"{'OCOH':<10s} {nmol:10d}")
lines.append(f"{'Protein':<10s} {1:10d}")

# Write new topology
with open("topol.top", "w") as f:
    for line in lines:
        f.write(line + "\n")

# Solvate system
gmx.solvate(cp="conf.gro", cs="spc216.gro", p="topol.top", o="conf.gro", backup=False)

# Neutralize system
gmx.grompp(f="mdp/min.mdp", c="conf.gro", p="topol.top", o="ions", maxwarn=5, backup=False)
gmx.genion(s="ions.tpr", p="topol.top", o="conf.gro", neutral=True, conc=0.05,
           pname="NA", pq=1, nname="CL", nq=-1, backup=False, input=("SOL"))

# Generate reference file for constraints
# Reference position is in the middle of Octanol slab
ref = u1.dimensions[0] 
u = mda.Universe("conf.gro")
octanol = u.select_atoms("resname OCOH")
positions = octanol.atoms.positions
positions[:, 0] = ref / 2
octanol.atoms.positions = positions
u.atoms.write("ref.gro")

# Generate restraint file for all heteroatoms in Octanol
id_list = list()
for atom in octanol.residues[0].atoms:
    if not atom.name.startswith("H"):
        id_list.append(atom.id)

# Generate file parameters
ftype = 2  # Flat-bottom potential
shape = 3  # In x-axis
width = ref / 2 / 10 * 1.05  # Half of Octanol slab width in [nm] + 5%
const = 5000  # kJ/mol K
with open("posres.itp", "w") as f:
    f.write("[ position_restraints ]\n")
    for id in id_list:
        f.write(f"{id:5d}{ftype:5d}{shape:5d}{width:6.2f}{const:6d}\n")

# Run minimization again
gmx.grompp(f="mdp/min.mdp", c="conf.gro", p="topol.top", o="min", nobackup=True)
gmx.mdrun(deffnm="min", v=True, nobackup=True)

# Run MD with flat-bottom potential
gmx.grompp(f="mdp/fbp.mdp", c="min.gro", p="topol.top",
           r="ref.gro", o="run", nobackup=True)
gmx.mdrun(deffnm="run", v=True, nobackup=True)
