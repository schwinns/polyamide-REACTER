# Script to call packmol and then convert to lammps with moltemplate

./packmol < initial.inp
moltemplate.sh -pdb system.pdb system.lt