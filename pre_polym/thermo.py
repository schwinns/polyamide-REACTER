# Script to extract thermodynamic data from LAMMPS log file

import pandas as pd
import matplotlib.pyplot as plt

def get_thermo(filename):
    # Read file
    reading_file = True
    f = open(filename, 'r')

    data = []
    while reading_file:
        # Get thermodynamic inputs and number of steps for run command
        for line in f:
            if len(line.split()) > 0 and not line.startswith('#'):

                if line.startswith('thermo_modify'):
                    print('Thermodynamic data with modified output: {}'.format(line.split()[1]))
                    print('See LAMMPS documentation about command thermo_modify:')
                    print('\thttps://docs.lammps.org/thermo_modify.html\n\n')

                elif line.startswith('thermo_style'):
                    if line.split()[1] == 'one':
                        thermo_style = 'one'
                        cols = ['step', 'temp', 'epair', 'emol', 'etotal', 'press']

                    elif line.split()[1] == 'multi':
                        thermo_style = 'multi'
                        cols = ['etotal', 'ke', 'temp', 'pe', 'ebond', 'eangle', 'edihed', 'eimp', 'evdwl', 'ecoul', 'elong', 'press']

                    elif line.split()[1] == 'yaml':
                        print('Thermodynamic data should be extracted as a yaml file. Use the following command line syntax:')
                        print("\tegrep  '^(keywords:|data:$|---$|\.\.\.$|  - \[)' log.lammps > log.yaml")
                        exit()

                    elif line.split()[1] == 'custom':
                        thermo_style = 'custom'
                        cols = line.split()[2:]

                elif line.startswith('thermo'):
                    save_every = int(line.split()[1])

                elif line.startswith('minimize'):
                    break

                elif line.startswith('run'):
                    if len(line.split()) > 2:
                        print('Extra run keywords... check your input script')
                        exit()
                    else:
                        n_steps = int(line.split()[1])
                        break

        # Find where thermodynamic data starts
        for line in f:
            if len(line.split()) > 0 and not line.startswith('#'):

                if thermo_style == 'one':
                    if line.startswith('Step'):
                        break
                    
                elif thermo_style == 'multi':
                    if line.startswith('TotEng'):
                        break

                elif thermo_style == 'custom':
                    if line.startswith(cols[0].capitalize()): # this is not very general (will only work with Step, Temp, Press as first column)
                        break

        # Read thermodynamic data
        for line in f:
            vals = line.split()
            
            if len(vals) == len(cols):
                vals = [float(x) for x in vals]
                data.append(vals)

            else:
                break

        next_line = f.readline()
        if not next_line:
            reading_file = False
            print('Finished reading file...')

    # Save data as a pandas dataframe
    df = pd.DataFrame(data, columns=cols)
    return df
