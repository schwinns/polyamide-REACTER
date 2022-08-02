# Script to move molecules broken across z boundary, so we can introduce walls

import matplotlib.pyplot as plt
from lammps_topology import LAMMPSTopology

# Read in topology
top = LAMMPSTopology()
top.read('pre_polym.data', skip_bonded=True, skip_velocities=True)

# Move molecules crossing the z boundary
zlo = top.box[2,0]
zhi = top.box[2,1]
z_box = zhi - zlo

max_z = -100
min_z = 100

for mol in top.Molecules:
    atom0 = mol.atoms[0]
    z0 = atom0.xyz[2]
    for atom in mol.atoms[1:]:
        z = atom.xyz[2]

        if (z0 - z) >= z_box / 2:
            atom.xyz[2] = z + z_box
            if (z + z_box) > max_z:
                max_z = z + z_box

        elif (z0 - z) <= -z_box / 2: 
            atom.xyz[2] = z - z_box
            if (z - z_box) < min_z:
                min_z = z - z_box

# Write a new file with new atom coordinates
original = top._original
out = open('polym_in.data', 'w')

    # Copy all header information except z limits of box
mass_split =  original.split('Masses')
for line in mass_split[0].split('\n'):
    if 'zlo' in line:
        out.write('{:.15f} {:.15f} zlo zhi\n'.format(min_z, max_z))

    else:
        out.write(line + '\n')

out.write('Masses')
atom_split = mass_split[1].split('Atoms')
out.write(atom_split[0])
out.write('Atoms\n\n')

    # Write Atoms section manually
for atom in top.Atoms:

    xyz = atom.xyz
    line = '{a_id} {mol_id} {atype} {charge} {x} {y} {z}\n'.format(a_id=atom.id, mol_id=atom.mol_id,
                                                                   atype=atom.atype, charge=atom.charge,
                                                                   x=xyz[0], y=xyz[1], z=xyz[2])
    out.write(line)

    # Copy all velocities and bonded information
velocity_split = atom_split[1].split('Velocities')
out.write('\nVelocities')
out.write(velocity_split[1])
out.close()
