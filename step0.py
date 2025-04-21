import shutil
import subprocess

import gromacs as gmx
import numpy as np

# Copy needed files
shutil.copy("topol/OCOH.gro", "conf.gro")
shutil.copy("topol/topol.top", "topol.top")

# Octanol properties
density = 824 # kg/m^3
mass = 0.130231 # kg/mol
box = np.array((7.0, 9.0, 9.0))  # nm
Na = 6.022E+23 # 1/mol

# Calculate number of molecules needed
nmol = round(np.prod(box * 1E-9) * density * Na / mass)

# Construct simulation box
gmx.insert_molecules(ci="conf.gro", o="conf.gro", box=(box * 1.1).tolist(),
                     nmol=nmol, _try=1000, backup=False)
with open("topol.top", "a") as f:
    f.write(f"{'OCOH':<10s} {nmol:10d}\n")

# Run minimization
gmx.grompp(f="mdp/min.mdp", c="conf.gro", p="topol.top", o="min", nobackup=True)
gmx.mdrun(deffnm="min", v=True, nobackup=True)

# Run equilibration
gmx.grompp(f="mdp/npt.mdp", c="min.gro", p="topol.top", o="run", nobackup=True)
gmx.mdrun(deffnm="run", v=True, nobackup=True)

# Remove periodic boundary
gmx.trjconv(f="run.gro", s="run.tpr", o="conf.gro", pbc="mol", nobackup=True, input=("System"))
