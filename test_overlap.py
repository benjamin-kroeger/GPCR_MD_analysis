# Paths to the PDB files of the two structures
pdb_file1 = '/home/benjaminkroeger/Documents/Master/Master_2_Semester/Internship/New Folder/3ODU_318.pdb'
pdb_file2 = '/home/benjaminkroeger/Documents/Master/Master_2_Semester/Internship/New Folder/3ODU_503.pdb'

import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Superimposer

# Function to calculate the van der Waals radius of an atom
def van_der_waals_radius(atom_name):
    vdw_radii = {
        'H': 1.20,
        'C': 1.70,
        'N': 1.55,
        'O': 1.52,
        'S': 1.80,
        # Add more elements as needed
    }
    return vdw_radii.get(atom_name.strip()[0], 1.70)  # Default radius for unknown elements is set to C

# Function to calculate the volume of a structure using van der Waals radii
def calculate_structure_volume(structure):
    atoms = [atom for atom in structure.get_atoms()]
    volume = sum(4/3 * np.pi * (van_der_waals_radius(atom.get_name()) ** 3) for atom in atoms)
    return volume

# Function to align two structures using Superimposer
def align_structures(structure1, structure2):
    atoms1 = [atom for atom in structure1.get_atoms()]
    atoms2 = [atom for atom in structure2.get_atoms()]

    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    sup.apply(atoms2)


# Parse the PDB files and get the structures
parser = PDBParser()
structure1 = parser.get_structure('structure1', pdb_file1)
structure2 = parser.get_structure('structure2', pdb_file2)

# Align the structures
align_structures(structure1, structure2)

# Calculate the volumes of each aligned structure using van der Waals radii
volume1 = calculate_structure_volume(structure1)
volume2 = calculate_structure_volume(structure2)

# Calculate the volume overlap (common volume) of the two aligned structures
overlap_volume = volume1 + volume2

# Calculate the total volume of both structures
total_volume = max(volume1, volume2)

# Calculate the volume overlap percentage
overlap_percentage = (overlap_volume / total_volume) * 100

print("Volume of Structure 1:", volume1)
print("Volume of Structure 2:", volume2)
print("Volume Overlap of Structures:", overlap_volume)
print("Overlap Percentage:", overlap_percentage, "%")
