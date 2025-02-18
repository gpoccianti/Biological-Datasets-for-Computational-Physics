# -*- coding: utf-8 -*-
"""ARC_TAU_new.ipynb

"""# Here begins the simulations proper"""

from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import *
from sys import stdout

import numpy as np
import mdtraj as mdt

import random

# Set random seed for reproducibility
random.seed(42)

"""## Download, fix missing atoms, solvate

Can also be done on the command line with the `pdbfixer` executable.
"""

###############################################
#Take the file as input
import sys

if len(sys.argv) < 2:
    print("Usage: python script.py <filename>")
    sys.exit(1)

filename = sys.argv[1]  # Get the filename from the command-line argument

print(f"Processing file: {filename}")
# Add your file processing logic here


# Load the PDB file with PDBFixer
fixer = PDBFixer(filename=filename+".pdb")

# Add missing (unresolved) residues. We don't want to model anything.
fixer.findMissingResidues()
fixer.missingResidues = {}
#fixer.addMissingResidues()

# Add missing (unresolved) atoms
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# Protonate (roughly) at chosen pH
fixer.addMissingHydrogens(pH=7.0)

positions = fixer.positions

# Save the file so it can be inspected
PDBFile.writeFile(fixer.topology, fixer.positions, open("fixed_"+filename+".pdb", "w"))

# Load your fixed complex PD
pdb = PDBFile("fixed_"+filename+".pdb")

# Load the force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Create a modeller object
modeller = Modeller(pdb.topology, pdb.positions)

# Create a system from the Modeller topology
system = forcefield.createSystem(modeller.topology)

# Step 3: Calculate the net charge of the system
nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
charges = []
for i in range(system.getNumParticles()):
    charge, sigma, epsilon = nonbonded.getParticleParameters(i)
    charges.append(charge)

print("Net Charge of the system: ", sum(charges))
####################################################

"""## Solvation
We use a dodecahedron solvation box while rotating the protein in the optimal conformation using Simulaid


"""## Packmol part

#We use this application to add the ions at a certain distance from the function


#Compute the Octahedron box
pdb = app.PDBFile("fixed_" + filename + ".pdb")
modeller = app.Modeller(pdb.topology, pdb.positions)

# Define padding
geompadding = 6.0 * nanometer

# Convert positions to NumPy array
positions = np.array([[p.x, p.y, p.z] for p in pdb.positions])

# Define unit vectors of the truncated octahedron box
vectors = np.array([
    [1, 0, 0],
    [1/3, 2*sqrt(2)/3, 0],
    [-1/3, sqrt(2)/3, sqrt(6)/3]
])

# Project positions onto each of the box vectors
projected_sizes = [
    np.max(np.dot(positions, vec)) - np.min(np.dot(positions, vec))
    for vec in vectors
]

# Select the largest projection as the defining dimension
max_projection = max(projected_sizes)

# Compute edge length of the truncated octahedron
edge_length = max_projection + geompadding.value_in_unit(nanometer)

print("Optimized truncated octahedron edge length:", edge_length)

# Convert box vectors to OpenMM format and scale them
boxVectors = [edge_length * Vec3(*vec) for vec in vectors]

import numpy as np
from Bio.PDB import PDBParser
from math import sqrt


# Define the padding and geometry parameters
#nanometer = 10  # Assuming this is how we convert to angstroms

# Function to calculate the maximum radius of the protein from the center of mass
def calculate_maximum_radius(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    # Get atom coordinates
    atom_coords = np.array([atom.coord for atom in structure.get_atoms()])

    # Compute the center of mass
    center_of_mass = np.mean(atom_coords, axis=0)

    # Compute the distance of each atom from the center of mass
    distances = np.linalg.norm(atom_coords - center_of_mass, axis=1)

    # The radius is the maximum distance from the center of mass
    protein_radius = np.max(distances)

    return protein_radius


# The radius for the inside and outside spheres
protein_radius_value = calculate_maximum_radius('fixed_complex.pdb')  # Replace with your protein PDB file
protein_radius_in_A = protein_radius_value  # Radius of the protein in Å

# Define the radius for the spheres (inside and outside)
padding = 10  # Padding around the protein in Å
inside_sphere_radius = protein_radius_in_A + padding  # Ions inside this sphere






# Define the dimensions of the truncated octahedron box
boxVectors = np.array([
    [1, 0, 0],
    [1/3, 2*sqrt(2)/3, 0],
    [-1/3, sqrt(2)/3, sqrt(6)/3]
])

# Compute the edge length of the box
positions = np.array([[p.x, p.y, p.z] for p in pdb.positions])  # Assuming pdb.positions exists
projected_sizes = [np.max(np.dot(positions, vec)) - np.min(np.dot(positions, vec)) for vec in boxVectors]
max_projection = max(projected_sizes)
edge_length = max_projection + geompadding.value_in_unit(nanometer) # Convert nanometer to Ångström
half_box_size_value = [0.5 * np.linalg.norm(vector) * edge_length for vector in boxVectors]
half_box_size = list(map(int, half_box_size_value))[0]

print("Edge length of the truncated octahedron:", edge_length)
print("Half box sizes (in Å):", half_box_size_value)

# Adjust the outside sphere radius to fit within the truncated octahedron
max_sphere_radius = half_box_size * sqrt(3) / 2  # Largest inscribed sphere in the box
outside_sphere_radius = min(max_sphere_radius, half_box_size + padding)


# Print calculated values for verification
print(f"Protein radius: {protein_radius_in_A:.2f} Å")
print(f"Inside sphere radius: {inside_sphere_radius:.2f} Å")
print(f"Outside sphere radius: {outside_sphere_radius:.2f} Å")

# Write the Packmol input file
inp_filename = "packmol_input.inp"
output_pdb = "ion_" + filename + ".pdb"  # Example filename, replace as needed
input_pdb = "fixed_" + filename + ".pdb"  # Example filename, replace as needed

with open(inp_filename, 'w') as inp_file:
    inp_file.write(f"""tolerance 2.0
filetype pdb
output {output_pdb}

structure {input_pdb}
  number 1
  fixed 0.0 0.0 0.0 0.0 0.0 0.0
end structure

structure na.pdb
  number 519
  inside sphere 0.0 0.0 0.0 {inside_sphere_radius}
  outside sphere 0.0 0.0 0.0 {outside_sphere_radius}
end structure

structure cl.pdb
  number 514
  inside sphere 0.0 0.0 0.0 {inside_sphere_radius}
  outside sphere 0.0 0.0 0.0 {outside_sphere_radius}
end structure
""")


print(f"Packmol input file saved as {inp_filename}")


print(f"Packmol input file saved as {inp_filename}")

#Generate na.pdb and cl.pdb files

# Define a simple PDB format for single ions
na_pdb_content = """\
HETATM    1 NA   NA     1      0.000   0.000   0.000  1.00  0.00           NA
END
"""

cl_pdb_content = """\
HETATM    1 CL   CL     1      0.000   0.000   0.000  1.00  0.00           CL
END
"""

# Save the files
with open("na.pdb", "w") as na_file:
    na_file.write(na_pdb_content)

with open("cl.pdb", "w") as cl_file:
    cl_file.write(cl_pdb_content)

print("Generated na.pdb and cl.pdb successfully!")



import subprocess

# Run Packmol
result = subprocess.run("packmol < packmol_input.inp", shell=True, capture_output=True, text=True)

# Print the output
print(result.stdout)
print(result.stderr)  # Print errors if any

"""## Without Simulaid"""

import mdtraj as mdt
import numpy as np
from scipy.linalg import eigh
from math import sqrt

# Load the PDB using MDTraj
traj = mdt.load("ion_" + filename + ".pdb")

# Identify protein atoms (excluding Na+ and Cl-)
protein_atoms = [atom.index for atom in traj.topology.atoms if atom.residue.name not in ["NA", "CL"]]

# Extract positions and masses for the protein only
protein_positions = traj.xyz[0][protein_atoms]
protein_masses = np.array([traj.topology.atom(i).element.mass for i in protein_atoms])

# Compute the center of mass for the protein
center_of_mass = np.sum(protein_positions * protein_masses[:, None], axis=0) / np.sum(protein_masses)

# Translate the whole system so the protein is centered
traj.xyz[0] -= center_of_mass  # Move all atoms (protein + ions)

# Compute the inertia tensor for the protein only
I = np.zeros((3, 3))
for i in range(len(protein_positions)):
    r = protein_positions[i]
    m = protein_masses[i]
    I += m * (np.dot(r, r) * np.eye(3) - np.outer(r, r))

# Compute the eigenvalues and eigenvectors (principal axes)
eigenvalues, eigenvectors = eigh(I)  # Sorted by increasing eigenvalue

# Ensure the largest principal axis is correctly selected
sorted_indices = np.argsort(eigenvalues)  # Sort indices
principal_axis = eigenvectors[:, sorted_indices[-1]]  # Largest eigenvalue

# Define the truncated octahedron box diagonal vector
unit_vectors = np.array([
    [1, 0, 0],
    [1/3, 2 * sqrt(2) / 3, 0],
    [-1/3, sqrt(2) / 3, sqrt(6) / 3]
])
diagonal = np.sum(unit_vectors, axis=0)
diagonal /= np.linalg.norm(diagonal)  # Normalize

# Compute rotation axis and angle
v = np.cross(principal_axis, diagonal)
s = np.linalg.norm(v)
c = np.dot(principal_axis, diagonal)

if s > 1e-6:  # Avoid division by zero (if already aligned)
    Vx = np.array([[0, -v[2], v[1]],
                   [v[2], 0, -v[0]],
                   [-v[1], v[0], 0]])
    R = np.eye(3) + Vx + np.dot(Vx, Vx) * ((1 - c) / (s ** 2))

    # Apply the same rotation to ALL atoms (protein + ions)
    traj.xyz[0] = np.dot(traj.xyz[0], R.T)

# Save the rotated PDB
output_filename = "aligned_" + filename + ".pdb"
traj.save(output_filename)

print(f"Box creation and rotation complete! Output: {output_filename}")

import numpy as np
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def truncated_octahedron_vertices():
    """
    Returns the vertices of a truncated octahedron centered at the origin with edge length 2.
    """
    base_vertex = np.array([0, np.sqrt(2), 2 * np.sqrt(2)])
    vertices = set()

    for perm in itertools.permutations(base_vertex):
        for signs in itertools.product([-1, 1], repeat=3):
            vertices.add(tuple(signs[i] * perm[i] for i in range(3)))

    return np.array(list(vertices))

def truncated_octahedron_faces(vertices):
    """
    Returns the faces of a truncated octahedron.
    """
    faces = []
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            for k in range(j + 1, len(vertices)):
                for l in range(k + 1, len(vertices)):
                    face = [vertices[i], vertices[j], vertices[k], vertices[l]]
                    if np.linalg.norm(face[0] - face[1]) == 2:
                        faces.append(face)
    return faces

def plot_truncated_octahedron(ax, vertices, faces):
    """
    Plots a truncated octahedron in a given 3D axis.
    """
    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    face_collection = Poly3DCollection(faces, alpha=0.5, edgecolor='k', facecolor='cyan')
    ax.add_collection3d(face_collection)

    ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], color='r', s=20)

def main():
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    vertices = truncated_octahedron_vertices()
    faces = truncated_octahedron_faces(vertices)
    plot_truncated_octahedron(ax, vertices, faces)
    plt.show()

if __name__ == "__main__":
    main()

from openmm import unit
import openmm.app as app
import numpy as np
from math import sqrt

# Load the aligned PDB file
pdb = app.PDBFile("fixed_" + filename + ".pdb")


# Extract atomic positions **excluding ions**
positions = np.array([
    [p.x, p.y, p.z] for atom, p in zip(pdb.topology.atoms(), pdb.positions)
    if atom.residue.name not in ["NA", "CL"]  # Exclude Na+ and Cl-
])

# Define unit vectors of the truncated octahedron box
vectors = np.array([
    [1, 0, 0],
    [1/3, 2*sqrt(2)/3, 0],
    [-1/3, sqrt(2)/3, sqrt(6)/3]
])

# Project positions onto each of the box vectors
projected_sizes = [
    np.max(np.dot(positions, vec)) - np.min(np.dot(positions, vec))
    for vec in vectors
]

# Select the largest projection as the defining dimension
max_projection = max(projected_sizes)

# Compute edge length of the truncated octahedron
edge_length = max_projection + geompadding.value_in_unit(nanometer)

print("Optimized truncated octahedron edge length:", edge_length)

# Load your protein and force field as before
forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
#pdb = app.PDBFile("fixed_" + filename + ".pdb")
pdb = app.PDBFile("aligned_" + filename + ".pdb")
modeller = app.Modeller(pdb.topology, pdb.positions)



# Convert positions to NumPy array
positions = np.array([ #Don't explicitly solvate ions
    [p.x, p.y, p.z] for atom, p in zip(pdb.topology.atoms(), pdb.positions)
    if atom.residue.name not in ["NA", "CL"]  # Exclude Na+ and Cl-
])

# Define unit vectors of the truncated octahedron box
vectors = np.array([
    [1, 0, 0],
    [1/3, 2*sqrt(2)/3, 0],
    [-1/3, sqrt(2)/3, sqrt(6)/3]
])

# Project positions onto each of the box vectors
projected_sizes = [
    np.max(np.dot(positions, vec)) - np.min(np.dot(positions, vec))
    for vec in vectors
]

# Select the largest projection as the defining dimension
max_projection = max(projected_sizes)

# Compute edge length of the truncated octahedron
edge_length = max_projection + geompadding.value_in_unit(nanometer)

print("Optimized truncated octahedron edge length:", edge_length)

# Convert box vectors to OpenMM format and scale them
boxVectors = [edge_length * Vec3(*vec) for vec in vectors]

# Add solvent using these box vectors
modeller.addSolvent(
    forcefield,
    model='tip3p',
    boxShape="octahedron",
    boxVectors=boxVectors,  # Use box vectors for the octahedron
    neutralize=False  # We already neutralized the system
)

PDBFile.writeFile(modeller.topology, modeller.positions, open("sol_" + filename + ".pdb", "w"))

print("Solvated system saved to sol_" + filename + ".pdb")

import MDAnalysis as mda

# Load your PDB file
u = mda.Universe("sol_"+filename+".pdb")

# Select water molecules
waters = u.select_atoms("resname HOH")  # Change "HOH" to "WAT" if needed

# Create a dictionary to store hexadecimal residue IDs
hex_residues = {}

# Renumber waters and store their resid in hexadecimal
for i, w in enumerate(waters):
    # Access the residue and set its resid (keep it as an integer)
    w.residue.resid = i + 1  # Resetting to start from 1
    # Store the hexadecimal representation in the dictionary
    hex_residues[w.residue.resname + str(i + 1)] = hex(i + 1)[2:]  # Store without '0x'

# Optionally renumber all residues if necessary and convert to hexadecimal
#for i, res in enumerate(u.residues):
#    res.resid = i + 1  # Resetting all residues to start from 1
#    hex_residues[res.resname + str(i + 1)] = hex(i + 1)[2:]  # Store without '0x'

# Save the modified structure
with mda.Writer("input_"+filename+"hex.pdb", u.atoms.n_atoms) as writer:
    writer.write(u.atoms)

# Print out the hexadecimal residue IDs for reference
#print("Hexadecimal Residue IDs:")
#for key, value in hex_residues.items():
#    print(f"{key}: {value}")


filename='complex'

def count_ions(pdb_file, ion_names=['NA', 'CL']):
    na_count = 0
    cl_count = 0

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith('HETATM'):
                if ' NA ' in line:
                    na_count += 1
                elif ' CL ' in line:
                    cl_count += 1

    return na_count, cl_count

# Check the number of ions in your solvated PDB
pdb_filename = "sol_" + filename  # Replace with your actual PDB filename
na_count, cl_count = count_ions(pdb_filename+".pdb")

print(f"Number of Na+ ions: {na_count}")
print(f"Number of Cl- ions: {cl_count}")

"""## Load PDB in MDTraj (see number of chains)"""

import mdtraj as mdt

pdb = mdt.load("fixed_"+filename+".pdb")

pdb.n_chains

for i,c in enumerate(pdb.topology.chains):
    print(f"Chain {i} has {c.n_residues} residues and {c.n_atoms} atoms.")
    r1 = c.residue(0)
    print(f"  First residue is {r1}")
