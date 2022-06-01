# Script to write moltemplate input files

# Imports
import numpy as np

def remove_comments(line):
    return line.split(';')[0]

class Topology():

    def __init__(self, topology_file):

        f = open(topology_file, 'r')

        # HEADER INFORMATION #
        for line in f:
            # Go to defaults directive
            if line.startswith('[ defaults ]'):
                break

        # DEFAULTS DIRECTIVE #
        for line in f:
            # Go to atomtypes directive
            if line.startswith('[ atomtypes ]'):
                break
            
            # Ignore comments and blank lines
            elif not line.startswith(';') and len(line.split()) != 0:

                line = remove_comments(line)
                if line.split()[2] == 'yes':
                    gen_pairs = True
                else:
                    gen_pairs = False

                self.defaults_directive = {
                    'nbfunc'    : int(line.split()[0]),
                    'comb-rule' : int(line.split()[1]),
                    'gen-pairs' : gen_pairs,
                    'fudgeLJ'   : float(line.split()[3]),
                    'fudgeQQ'   : float(line.split()[4])
                }

        # ATOMTYPES DIRECTIVE #
        self.atomtypes = {}
        for line in f:
            # Go to moleculetype directive
            if line.startswith('[ moleculetype ]'):
                break

            # Ignore comments and blank lines
            elif not line.startswith(';') and len(line.split()) != 0:

                line = remove_comments(line)
                atom = line.split()
                self.atomtypes[atom[0]] = {
                    'at_num'  : int(atom[1]),
                    'mass'    : float(atom[2]),
                    'charge'  : float(atom[3]),
                    'ptype'   : atom[4],
                    'sigma'   : float(atom[5]),
                    'epsilon' : float(atom[6]) 
                }

        # MOLECULETYPE DIRECTIVE #
        self.moleculetype = {}
        for line in f:
            # Go to atoms directive
            if line.startswith('[ atoms ]'):
                break

            # Ignore comments and blank lines
            elif not line.startswith(';') and len(line.split()) != 0:

                line = remove_comments(line)
                self.moleculetype[line.split()[0]] = {
                    'nrexcl'  : int(line.split()[1])
                }

        # ATOMS DIRECTIVE #
        self.atoms = {}
        for line in f:
            # Go to bonds directive
            if line.startswith('[ bonds ]'):
                break

            # Ignore comments and blank lines
            elif not line.startswith(';') and len(line.split()) != 0:

                line = remove_comments(line)
                atom = line.split()
                atom_id = int(atom[0])
                self.atoms[atom_id] = {
                    'type'      : atom[1],
                    'resnr'     : int(atom[2]),
                    'residue'   : atom[3],
                    'atom_name' : atom[4],
                    'cgnr'      : int(atom[5]),
                    'charge'    : float(atom[6]),
                    'mass'      : float(atom[7])
                }

        # BONDS DIRECTIVE #
        self.bonds = {}
        bond_id = 0
        for line in f:
            # Go to angles directive
            if line.startswith('[ pairs ]'):
                break

            # Ignore comments and blank lines
            elif not line.startswith(';') and len(line.split()) != 0:

                line = remove_comments(line)
                bond = line.split()
                bond_id += 1
                self.bonds[bond_id] = {
                    'a1'     : int(bond[0]),
                    'a2'     : int(bond[1]),
                    'func'   : int(bond[2]),
                    'params' : [float(b) for b in bond[3:]]
                }


################# TESTING #################
topology_file = './MPD_TMC_moltemplate/MPD-L.top'
MPD = Topology(topology_file)

print('self.defaults_directive')
print('\t', MPD.defaults_directive)

print('self.atomtypes')
for a in MPD.atomtypes:
    print('\t', a, MPD.atomtypes[a])

print('self.moleculetype')
for m in MPD.moleculetype:
    print('\t', m, MPD.moleculetype[m])

print('self.atoms')
for a in MPD.atoms:
    print('\t', a, MPD.atoms[a])

print('self.bonds')
for b in MPD.bonds:
    print('\t', b, MPD.bonds[b])