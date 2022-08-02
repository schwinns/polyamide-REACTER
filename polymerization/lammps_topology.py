# Creates a lammps topology class with all necessary information for a lammps data topology

import numpy as np
from time import time

def remove_comments(line):
    return line.split('#')[0].lstrip()


def isfloat(string):
    '''
    Test if a string can be expressed as a float. If it can, return the float.
    '''
    try:
        float(string)
        return float(string)
    except ValueError:
        return string


class LAMMPSTopology():

    def __init__(self):

        # Save filename and all original lines
        self._filename = None
        self._original = None

        # Save header information
        self.box = np.array([[0.0,0.0],[0.0,0.0],[0.0,0.0]])
        self.masses = {}
        self.pair_coeffs = {}
        self.bond_coeffs = {}
        self.angle_coeffs = {}
        self.dihedral_coeffs = {}
        self.improper_coeffs = {}

        # Define lists of Atom objects and Molecule objects
        self.Atoms = []
        self.Molecules = []

        # Lists of all molecules and atoms
        self._atoms = []
        self._molecules = []


    def atom(self, index):
        '''
        Get Atom object from index in Atoms list
        '''
        return self.Atoms[index]


    def atom_by_id(self, a_id):
        '''
        Get Atom object from atom id in data file
        '''
        idx = self._atoms.index(a_id)  
        return self.Atoms[idx]
                
    
    # TO DO
    # def atom_by_name(self, name):
    #     '''
    #     Get Atom object from name defined as 'molecule id'-'atom type'-'number in group'
    #     '''
    #     for a in self.Atoms:
    #         if a.name == name:
    #             return a

    
    def molecule(self, index):
        '''
        Get Molecule object from index in Atoms list
        '''        
        return self.Molecules[index]


    def mol_by_id(self, mol_id):
        '''
        Get Molecule object from molecule id in data file
        '''
        idx = self._molecules.index(mol_id)
        return self.Molecules[idx]


    class Atom():
        '''
        Class to store all atom information as identified by id
        '''
 
        def __init__(self, Topology, a_id, mol_id, atype, charge, xyz):
            
            self.Topology = Topology # access outer class LAMMPSTopology

            self.id = a_id
            self.mol_id = mol_id
            self.atype = atype
            self.charge = charge
            self.xyz = xyz
            self.mass = self.Topology.masses[atype]
            self.velocity = np.zeros(3)
            self.bonds = {}


        def molecule(self):
            '''
            Get molecule of the Atom
            '''
            return self.Topology.mol_by_id(self.mol_id)



    class Molecule():
        '''
        Class to store all molecule information as identified id
        '''
 
        def __init__(self, Topology, atom_s, mol_id):

            self.Topology = Topology # access outer class LAMMPSTopology
            self.id = mol_id

            self.atoms = []
            self.atoms_by_id = []
            self.mass = 0

            self.bonds = {}
            self.angles = {}
            self.dihedrals = {}
            self.impropers = {}

            if atom_s is list:
                for A in atom_s:
                    self.atoms.append(A)
                    self.atoms_by_id.append(A.id)
                    self.mass += A.mass
            else:
                self.atoms.append(atom_s)
                self.atoms_by_id.append(atom_s.id)
                self.mass += atom_s.mass


        def add_atom(self, atom):
            '''
            Add an Atom object to the Molecule atoms list
            '''
            self.atoms.append(atom)
            self.atoms_by_id.append(atom.id)
            self.mass += atom.mass

        
        def calculate_com(self):
            '''
            Calculate center of mass of molecule and save as Attribute
            '''
            self.com = np.zeros(3)

            for A in self.atoms:
                self.com += A.xyz * A.mass

            self.com = self.com / self.mass

            return self.com


    def read(self, filename, skip_bonded=True, skip_velocities=True):
        '''
        Method to read an input data file

        Inputs:
            filename (string)      -- name of data file to read in
            skip_bonded (bool)     -- whether to read in bonds, angle, dihedrals, and impropers : default = True
            skip_velocities (bool) -- whether to read in velocities : default = True
        '''

        start = time()
        self.filename = filename
        f = open(filename, 'r')
        self._original = f.read()
        f.close()
        f = open(filename, 'r')

    # HEADER INFORMATION
        self.header = f.readline()
        for line in f:

            line = remove_comments(line)
            l = line.split()
            
    # Number of atoms, bonds, angles, dihedrals, impropers
            if len(l) == 2:

                n = int(l[0])
                obj = l[1]

                if obj == 'atoms':
                    self.n_atoms = n
                elif obj == 'bonds':
                    self.n_bonds = n
                elif obj == 'angles':
                    self.n_angles = n
                elif obj == 'dihedrals':
                    self.n_dihedrals = n
                elif obj == 'impropers':
                    self.n_impropers = n
                else:
                    raise Exception('Unrecognized header information: {}'.format(line))

    # Number of atom, bond, etc. types
            elif len(l) == 3:

                n = int(l[0])
                obj = l[1]

                if obj == 'atom':
                    self.atom_types = [0]*n
                    n_atypes = 0
                elif obj == 'bond':
                    self.bond_types = [0]*n
                    n_btypes = 0
                elif obj == 'angle':
                    self.angle_types = [0]*n
                    n_angtypes = 0
                elif obj == 'dihedral':
                    self.dihedrals_types = [0]*n
                    n_dihtypes = 0
                elif obj == 'improper':
                    self.impropers_types = [0]*n
                    n_imptypes = 0
                else:
                    raise Exception('Unrecognized header information: {}'.format(line))

    # Box vectors
            elif len(l) == 4:

                lo = float(l[0])
                hi = float(l[1])
                obj = l[2]

                if obj == 'xlo':
                    self.box[0,0] = lo
                    self.box[0,1] = hi
                elif obj == 'ylo':
                    self.box[1,0] = lo
                    self.box[1,1] = hi
                elif obj == 'zlo':
                    self.box[2,0] = lo
                    self.box[2,1] = hi
                else:
                    raise Exception('Unrecognized box definition: {}'.format(line))

            # Go to Masses section
            elif line.startswith('Masses'):
                break

            # Skip blank lines
            elif len(l) == 0:
                pass

            # Print unused lines
            else:
                print('Unused line:\n\t{}'.format(line))

    # Masses
        for line in f:

            line = remove_comments(line)
            l = line.split()

            # Go to Pair Coeffs
            if line.startswith('Pair'):
                atom = False
                break

            # Go to Atoms section
            elif line.startswith('Atoms'):
                atom = True
                break

            # Masses
            elif len(l) == 2:

                atype = int(l[0])
                m = float(l[1])

                self.masses[atype] = m

    # Pair Coeffs
        if not atom:
            for line in f:

                line = remove_comments(line)
                l = line.split()

                # Go to Bond Coeffs
                if line.startswith('Bond'):
                    atom = False
                    break

                # Go to Atoms section
                elif line.startswith('Atoms'):
                    atom = True
                    break

                # Pair Coeffs
                elif len(l) == 3:

                    at = int(l[0])
                    eps = float(l[1])
                    sig = float(l[2])

                    self.pair_coeffs[at] = [eps, sig]

                # Skip blank lines
                elif len(l) == 0:
                    pass

                else:
                    raise Exception('Unrecognized pair coeffient form: {}'.format(line))

    # Bond Coeffs
        if not atom:
            for line in f:

                line = remove_comments(line)
                l = line.split()

                # Go to Angle Coeffs
                if line.startswith('Angle'):
                    atom = False
                    break

                # Go to Atoms section
                elif line.startswith('Atoms'):
                    atom = True
                    break

                # Cannot read other coefficients now
                elif len(l) == 2:
                    atom = True
                    print('WARNING: Cannot read {} parameters. Continuing to Atoms section.'.format(line.split('\n')[0]))
                    break

                # Skip blank lines and save coefficients
                elif len(l) > 1:

                    at = int(l[0])
                    self.bond_coeffs[at] = []

                    for coeff in l[1:]:
                        self.bond_coeffs[at].append(isfloat(coeff))

    # Angle Coeffs
        if not atom:
            for line in f:

                line = remove_comments(line)
                l = line.split()

                # Go to Dihedral Coeffs
                if line.startswith('Dihedral'):
                    atom = False
                    break

                # Go to Atoms section
                elif line.startswith('Atoms'):
                    atom = True
                    break

                # Cannot read other coefficients now
                elif len(l) == 2:
                    atom = True
                    print('WARNING: Cannot read {} parameters. Continuing to Atoms section.'.format(line.split('\n')[0]))
                    break

                # Skip blank lines and save coefficients
                elif len(line.split()) > 1:

                    at = int(l[0])
                    self.angle_coeffs[at] = []

                    for coeff in l[1:]:
                        self.angle_coeffs[at].append(isfloat(coeff))

    # Dihedral Coeffs
        if not atom:
            for line in f:

                line = remove_comments(line)
                l = line.split()

                # Go to Dihedral Coeffs
                if line.startswith('Improper'):
                    atom = False
                    break

                # Go to Atoms section
                elif line.startswith('Atoms'):
                    atom = True
                    break

                # Cannot read other coefficients now
                elif len(l) == 2:
                    atom = True
                    print('WARNING: Cannot read {} parameters. Continuing to Atoms section.'.format(line.split('\n')[0]))
                    break

                # Skip blank lines and save coefficients
                elif len(l) > 1:

                    at = int(l[0])
                    self.dihedral_coeffs[at] = []

                    for coeff in l[1:]:
                        self.dihedral_coeffs[at].append(isfloat(coeff))

    # Improper Coeffs
        if not atom:
            for line in f:

                line = remove_comments(line)
                l = line.split()

                # Go to Atoms section
                if line.startswith('Atoms'):
                    atom = True
                    break

                # Cannot read other coefficients now
                elif len(l) == 2:
                    atom = True
                    print('WARNING: Cannot read {} parameters. Continuing to Atoms section.'.format(line.split('\n')[0]))
                    break

                # Skip blank lines and save coefficients
                elif len(l) > 1:

                    at = int(l[0])
                    self.dihedral_coeffs[at] = []

                    for coeff in l[1:]:
                        self.dihedral_coeffs[at].append(isfloat(coeff))

    # Other Coeffs (will not be read)
        if not atom:
            print('WARNING: Cannot read further parameters. Continuing to Atoms section.')

        header_time = time()
    # Atoms section
        if atom:
            for line in f:

                line = remove_comments(line)
                l = line.split()

                # Go to Velocities section
                if line.startswith('Velocities'):
                    bond = False
                    break

                # Go to Bonds section
                elif line.startswith('Bonds'):
                    bond = True
                    break

                # Cannot read other sections
                elif len(l) > 0 and len(l) < 7:
                    bond = True
                    print('WARNING: Cannot read {} section. Continuing to Bonds section.'.format(line.split('\n')[0]))
                    break

                # Skip blank lines and save atom information
                elif len(l) >= 7:

                    a_id = int(l[0])
                    mol_id = int(l[1])
                    atype = int(l[2])
                    charge = float(l[3])
                    x = float(l[4])
                    y = float(l[5])
                    z = float(l[6])
                    # name = mol_id + '-' + atype + '-' + at_num # TO DO add name to Atom class

                    if atype not in self.atom_types:
                        self.atom_types[n_atypes] = atype
                        n_atypes += 1

                    self._atoms.append(a_id)
                    A = self.Atom(self, a_id, mol_id, atype, charge, np.array([x,y,z]))
                    self.Atoms.append(A)

                    if mol_id not in self._molecules: # create Molecule object
                        self._molecules.append(mol_id)
                        M = self.Molecule(self, A, mol_id)
                        self.Molecules.append(M)

                    else:                             # add atom to existing Molecule object
                        M = self.mol_by_id(mol_id)
                        M.add_atom(A)

        atom_time = time()
        print('Finished reading atoms...')

        if not skip_velocities:
    # Velocities section
            if not bond:
                for line in f:

                    line = remove_comments(line)
                    l = line.split()

                    # Go to Bonds section
                    if line.startswith('Bonds'):
                        bond = True
                        break

                    # Cannot read other sections
                    elif len(l) > 0 and len(l) < 4:
                        bond = True
                        print('WARNING: Cannot read {} section. Continuing to Bonds section.'.format(line.split('\n')[0]))
                        break

                    # Velocity section
                    elif len(l) == 4:
                        
                        a_id = int(l[0])
                        vx = float(l[1])
                        vy = float(l[2])
                        vz = float(l[3])

                        A = self.atom_by_id(a_id)
                        A.velocity = np.array([vx,vy,vz])

            velocity_time = time()
            print('Finished reading velocitites...')

        if not skip_bonded:
    # Bonds section
            if bond:
                for line in f:

                    line = remove_comments(line)
                    l = line.split()

                    # Go to Angles section
                    if line.startswith('Angles'):
                        angle = True
                        break

                    # Cannot read other sections
                    elif len(l) > 0 and len(l) < 4:
                        angle = True
                        print('WARNING: Cannot read {} section. Continuing to Angles section.'.format(line.split('\n')[0]))
                        break

                    # Bonds section
                    elif len(l) == 4:

                        b_id = int(l[0])
                        btype = int(l[1])
                        a1 = int(l[2])
                        a2 = int(l[3])

                        A1 = self.atom_by_id(a1)
                        A1.bonds[b_id] = {
                            'type'  : btype,
                            'atoms' : [a1, a2] 
                        }

                        A2 = self.atom_by_id(a2)
                        A2.bonds[b_id] = {
                            'type'  : btype,
                            'atoms' : [a1, a2] 
                        }

                        M = self.mol_by_id(A1.mol_id)
                        M.bonds[b_id] = {
                            'type'  : btype,
                            'atoms' : [a1, a2] 
                        }

            bond_time = time()
            print('Finished reading bonds...')
        # Angles section
            if angle:
                for line in f:

                    line = remove_comments(line)
                    l = line.split()

                    # Go to Dihedrals section
                    if line.startswith('Dihedrals'):
                        dih = True
                        break

                    # Cannot read other sections
                    elif len(l) > 0 and len(l) < 5:
                        dih = True
                        print('WARNING: Cannot read {} section. Continuing to Dihedrals section.'.format(line.split('\n')[0]))
                        break

                    # Angles section
                    elif len(l) == 5:

                        ang_id = int(l[0])
                        angtype = int(l[1])
                        a1 = int(l[2])
                        a2 = int(l[3])
                        a3 = int(l[4])

                        A1 = self.atom_by_id(a1)
                        M = self.mol_by_id(A1.mol_id)
                        M.angles[ang_id] = {
                            'type'  : angtype,
                            'atoms' : [a1, a2, a3]
                        }

            angle_time = time()
            print('Finished reading angles...')
        # Dihedrals section
            if dih:
                for line in f:

                    line = remove_comments(line)
                    l = line.split()

                    # Go to Impropers section
                    if line.startswith('Impropers'):
                        imp = True
                        break

                    # Cannot read other sections
                    elif len(l) > 0 and len(l) < 6:
                        imp = True
                        print('WARNING: Cannot read {} section. Continuing to Impropers section.'.format(line.split('\n')[0]))
                        break

                    # Dihedrals section
                    elif len(l) == 6:

                        dih_id = int(l[0])
                        dihtype = int(l[1])
                        a1 = int(l[2])
                        a2 = int(l[3])
                        a3 = int(l[4])
                        a4 = int(l[5])

                        A1 = self.atom_by_id(a1)
                        M = self.mol_by_id(A1.mol_id)
                        M.dihedrals[dih_id] = {
                            'type'  : dihtype,
                            'atoms' : [a1, a2, a3, a4]
                        }

            dihedral_time = time()
            print('Finished reading dihedrals...')
        # Impropers section
            if imp:
                for line in f:

                    line = remove_comments(line)
                    l = line.split()

                    # Cannot read other sections
                    if len(l) > 0 and len(l) < 6:
                        print('WARNING: Cannot read {} section. Gathered all possible information'.format(line.split('\n')[0]))
                        break

                    # Impropers section
                    elif len(l) == 6:

                        imp_id = int(l[0])
                        imptype = int(l[1])
                        a1 = int(l[2])
                        a2 = int(l[3])
                        a3 = int(l[4])
                        a4 = int(l[5])

                        A1 = self.atom_by_id(a1)
                        M = self.mol_by_id(A1.mol_id)
                        M.dihedrals[imp_id] = {
                            'type'  : imptype,
                            'atoms' : [a1, a2, a3, a4]
                        }

            improper_time = time()
            print('Finished reading impropers...')

        end_time = time()
        f.close()
        print('\n---------- FILE READING TIMING ---------')
        print('Section\t\t\t\tTime (s)\n')
        print('Header\t\t\t\t{:>.4f}'.format(header_time - start))
        print('Atoms\t\t\t\t{:>.4f}'.format(atom_time - header_time))
        if not skip_velocities:
            print('Velocities\t\t\t{:>.4f}'.format(velocity_time - atom_time))
        if not skip_bonded:
            print('Bonds\t\t\t\t{:>.4f}'.format(bond_time - velocity_time))
            print('Angles\t\t\t\t{:>.4f}'.format(angle_time - bond_time))
            print('Dihedrals\t\t\t{:>.4f}'.format(dihedral_time - angle_time))
            print('Impropers\t\t\t{:>.4f}'.format(improper_time - dihedral_time))
        print('\nTotal time\t\t\t{:>.4f}'.format(end_time - start))
        print('-'*41 + '\n')


    def write(self, filename, keep_header=True, keep_bonded=False):
        '''
        Method to write a data file from LAMMPSTopology object
        '''

        original = self._original

        out = open(filename, 'w')
        out.close()


















        
