# Script to write moltemplate input files

# Imports
from topology_class import Topology
import mdtraj as md

cwd = '/home/nate/Projects/polymer_membrane/polyamide-REACTER/examples/MPD_TMC_moltemplate/'
gro = 'MPD-L.gro'
top = 'MPD-L.top'

output = 'MPD-L.lt'

coordinate = md.load(cwd + gro, top=cwd + gro) # mdtraj atoms are 0-indexed
topology = Topology(cwd + top) # Topology atoms are 1-indexed
coords = coordinate.xyz[0,:,:]

# Write the name of Moltemplate object
mol = output.split('.')[0]
out = open(cwd + output, 'w')
out.write(mol + ' {\n\n')

# Data Atoms section
out.write('\twrite("Data Atoms") {\n')
out.write('\t\t# atomID\tmolID\tatomType\tcharge\tx\ty\tz\n')

i = 0
for atom in coordinate.top.atoms:

    top_atom = topology.atoms[i+1]

    atomID = top_atom['atom_name']
    molID = top_atom['residue']
    atomType = top_atom['type']
    charge = top_atom['charge']
    x = coords[i,0]*10
    y = coords[i,1]*10
    z = coords[i,2]*10

    line = '\t\t$atom:{:<2}\t$mol:{}\t@atom:{}\t{}\t{}\t{}\t{}\n'.format(atomID, molID, atomType, charge, x, y, z)
    out.write(line)
    i += 1

out.write('\t}\n\n')

# Data Masses section
out.write('\twrite_once("Data Masses") {\n')

i = 0
for atom in coordinate.top.atoms:
    
    i += 1
    top_atom = topology.atoms[i]

    atomType = top_atom['type']
    mass = top_atom['mass']

    line = '\t\t@atom:{:<2}\t{}\n'.format(atomType, mass)
    out.write(line)

out.write('\t}\n\n')

# Data Bonds section
out.write('\twrite("Data Bonds") {\n')
out.write('\t\t# bondID\tbondType\tatomID1\tatomID2\n')

for i in topology.bonds:

    bond = topology.bonds[i]

    a1 = bond['a1']
    a2 = bond['a2']
    
    a1_name = topology.atoms[a1]['atom_name']
    a2_name = topology.atoms[a2]['atom_name']
    
    element1 = a1_name.strip('0123456789')
    element2 = a2_name.strip('0123456789')
    
    bond_name = element1 + element2

    line = '\t\t$bond:b{}\t@bond:{}\t$atom:{:<2}\t$atom:{:<2}\n'.format(i, bond_name, a1_name, a2_name)
    out.write(line)

out.write('\t}\n\n')

# Data Angles section
out.write('\twrite("Data Angles") {\n')
out.write('\t\t# angleID\tangleType\tatomID1\tatomID2\tatomID3\n')

for i in topology.angles:

    angle = topology.angles[i]

    a1 = angle['a1']
    a2 = angle['a2']
    a3 = angle['a3']

    a1_name = topology.atoms[a1]['atom_name']
    a2_name = topology.atoms[a2]['atom_name']
    a3_name = topology.atoms[a3]['atom_name']
    
    element1 = a1_name.strip('0123456789')
    element2 = a2_name.strip('0123456789')
    element3 = a3_name.strip('0123456789')
    
    angle_name = element1 + element2 + element3

    line = '\t\t$angle:a{}\t@angle:{}\t$atom:{:<2}\t$atom:{:<2}\t$atom:{:<2}\n'.format(i, angle_name, a1_name, a2_name, a3_name)
    out.write(line)

out.write('\t}\n\n')

# Data Dihedrals section
out.write('\twrite("Data Dihedrals") {\n')
out.write('\t\t# dihID\tdihType\tatomID1\tatomID2\tatomID3\tatomID4\n')

for i in topology.dihedrals:

    dih = topology.dihedrals[i]

    a1 = dih['a1']
    a2 = dih['a2']
    a3 = dih['a3']
    a4 = dih['a4']
    
    a1_name = topology.atoms[a1]['atom_name']
    a2_name = topology.atoms[a2]['atom_name']
    a3_name = topology.atoms[a3]['atom_name']
    a4_name = topology.atoms[a4]['atom_name']
    
    element1 = a1_name.strip('0123456789')
    element2 = a2_name.strip('0123456789')
    element3 = a3_name.strip('0123456789')
    element4 = a4_name.strip('0123456789')

    dih_name = element1 + element2 + element3 + element4
    if dih['func'] in [2,4]: # mark improper dihedrals
        dih_name += '_imp'
    
    line = '\t\t$dihedral:a{}\t@dihedral:{}\t$atom:{:<2}\t$atom:{:<2}\t$atom:{:<2}\t$atom:{:<2}\n'.format(i, dih_name, a1_name, a2_name, a3_name, a4_name)
    out.write(line)

out.write('\t}\n\n')

# In Settings section
out.write('\twrite_once("In Settings") {\n\n')

#    Pair coefficients
out.write('\t\t# Pair Coeffs\n')
out.write('\t\t# \tatomType1\tatomType2\tparams (epsilon, sigma)\n')
for i in topology.atomtypes:

    atype = topology.atomtypes[i]

    eps = atype['epsilon']/4.184 # units need to be kcal/mol
    sig = atype['sigma']*10      # units need to be Angstroms

    line = '\t\tpair_coeff\t@atom:{}\t@atom:{}\t{}\t{}\n'.format(i, i, eps, sig)
    out.write(line)

#   Bond Coefficients (assuming all bonds are harmonic)
out.write('\n\t\t# Bond Coeffs\n')
out.write('\t\t# \tbondType\tparams (k_bond, r0)\n')
for i in topology.bonds:

    bond = topology.bonds[i]

    a1 = bond['a1']
    a2 = bond['a2']
    
    a1_name = topology.atoms[a1]['atom_name']
    a2_name = topology.atoms[a2]['atom_name']
    
    element1 = a1_name.strip('0123456789')
    element2 = a2_name.strip('0123456789')
    
    bond_name = element1 + element2

    k_bond = bond['params'][1]/4.184/10**2 # units need to be kcal/mol/Ang^2
    r0 = bond['params'][0]*10              # units need to be Angstroms

    line = '\t\tbond_coeff\t@bond:{}\t{}\t{}\n'.format(bond_name, k_bond, r0)
    out.write(line)






out.write('\t}\n\n')