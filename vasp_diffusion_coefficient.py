import glob
from pymatgen.core import Structure
from numpy.lib.function_base import _calculate_shapes
from sklearn.metrics import r2_score
from analysis import DiffusionCoefficient
from ase import Atom, Atoms
from ase.cell import Cell
from ase.io.trajectory import Trajectory, TrajectoryReader
from ase.io import read
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ase.units import Bohr, Rydberg, kJ, kB, fs, Hartree, mol, kcal
from ase.geometry.analysis import Analysis
from pymatgen.io.vasp import Poscar
import argparse

plt.rcParams.update({'font.size': 30})

do_bonding = False
do_rdf = False
calculation_type = 'tracer'

ignore_n_images = 0

dt = 1
axes = ['all', 'x', 'y', 'z', 'xy']
axes = axes[0]

# Operations to perform on subsequent POSCARs to
# conform to previous POSCARs.


def get_operation_matrix(ref_poscar_str):
    poscar = Poscar.from_string(ref_poscar_str)
    structure = poscar.structure
    fc = structure.frac_coords
    om = np.ndarray(fc.shape)

    for i in range(fc.shape[0]):
        for j in range(fc.shape[1]):
            if fc[i][j] < 0:
                om[i][j] = -(int(fc[i][j]) + 1)
            elif fc[i][j] > 1:
                om[i][j] = int(fc[i][j])
            else:
                om[i][j] = 0
    return om


def string_to_ase(poscar_str, convert=1):
    poscar = Poscar.from_string(poscar_str)
    structure = poscar.structure
    fc = structure.frac_coords
    lattice = structure.lattice
    cc = structure.cart_coords

    a = Atoms(scaled_positions=fc, numbers=structure.atomic_numbers, pbc=True, cell=Cell.fromcellpar(
        [lattice.a*convert,
         lattice.b*convert,
         lattice.c*convert,
         lattice.alpha,
         lattice.beta,
         lattice.gamma]))

    return a


def string_to_ase_with_operation(poscar_str, om, convert=1):
    poscar = Poscar.from_string(poscar_str)
    structure = poscar.structure
    fc = structure.frac_coords
    fc += om
    lattice = structure.lattice
    cc = structure.cart_coords

    a = Atoms(scaled_positions=fc, numbers=structure.atomic_numbers, pbc=True, cell=Cell.fromcellpar(
        [lattice.a*convert,
         lattice.b*convert,
         lattice.c*convert,
         lattice.alpha,
         lattice.beta,
         lattice.gamma]))

    return a


print('Calculations for', calculation_type)

DiffusionCoeff = {}

diff_coeff_file = open('diff_coeffs_'+calculation_type + '_' + axes, 'w')
trajectory_files = glob.glob('XDATCAR_*')
trajectory_files.sort()

print('Trajectories:', trajectory_files)
trajectories_combined = []
trajectories_lines_combined = []

for i_trajectory_file in range(len(trajectory_files)):
    trajectory_file = trajectory_files[i_trajectory_file]
    f = open(trajectory_file, 'r')
    temp_trajectory_list = f.readlines()
    f.close()

    trajectories_lines_combined += [temp_trajectory_list]

    lattice_list = temp_trajectory_list[0:7]
    lattice_str = ''.join(lattice_list)
    atom_counts = [int(x) for x in temp_trajectory_list[6].split()]
    num_atoms = sum(atom_counts)
    number_of_lines = num_atoms + 1
    trajectory = temp_trajectory_list[7:]
    tot_num_images = int(len(trajectory)/number_of_lines)

    print('file: ', trajectory_file,
          'total number of images:', str(tot_num_images))

    trajectory_partial = trajectory

    trajectory_list = []
    trajectory_list_ang = []
    num_steps = int(len(trajectory_partial)/number_of_lines)

    list_of_added = []

    initial = string_to_ase(lattice_str + ''.join(
        temp_trajectory_list[7:7+num_atoms+1]), 1e-8)
    initial_ang = string_to_ase(lattice_str + ''.join(
        temp_trajectory_list[7:7+num_atoms+1]))

    list_of_added = []
    atomic_numbers = np.unique(initial.get_atomic_numbers())

    for i in range(0, num_steps):
        positions_str = ''.join(
            temp_trajectory_list[7+i*(num_atoms+1):7+(i+1)*(num_atoms+1)])
        a = string_to_ase(
            lattice_str+positions_str, 1e-8)
        a_ang = string_to_ase(lattice_str+positions_str)
        trajectory_list += [a]
        trajectory_list_ang += [a_ang]

    trajectories_combined += [trajectory_list]

trajectory_list = []
for t in trajectories_combined:
    for tt in t:
        trajectory_list += [tt]

if do_bonding:
    bonds_array = np.zeros([len(atomic_numbers), len(atomic_numbers)])
    for i in range(0, num_steps):
        # bonding

        D = a_ang.get_all_distances()

        bonds = []
        for i in atomic_numbers:
            i_array = np.where([(a_ang.get_atomic_numbers() == i)])[1]
            bonds_row = []
            for j in atomic_numbers:
                j_array = np.where([(a_ang.get_atomic_numbers() == j)])[1]
                av = 0
                count = 0
                for ii in i_array:
                    for jj in j_array:
                        if ii != jj:
                            d = a_ang.get_distance(ii, jj)
                            if d < 3:
                                av += d
                                count += 1
                if count > 0:
                    av = av/count
                else:
                    av = 0
                bonds_row += [av]
            bonds += [bonds_row]

        bonds = np.array(bonds)

        bonds_array += bonds

    bonds_array = bonds_array/num_steps

    bonds_array_df = pd.DataFrame(bonds_array)
    bonds_array_df.to_csv('data_bonds_' + calculation_type+'_'+axes+'.csv')

# rdf calculation
if do_rdf:
    ga = Analysis(trajectory_list_ang)
    bins = 100
    print('created Analysis')
    min_lattice_length = a_ang.cell.lengths().min()
    rdf = ga.get_rdf(min_lattice_length/2, bins,
                     elements=['Li'], return_dists=True)

    X = rdf[0][1]
    average_rdf = np.zeros(bins)
    for i_rdf in rdf:
        average_rdf += i_rdf[0]
    average_rdf /= len(rdf)

    plt.figure(figsize=(15, 10))

    plt.xlabel('r (Å)')
    plt.ylabel('g$_{Li-Li}$(r)')
    plt.plot(X, average_rdf)

    plt.savefig('rdf_' + calculation_type+'_'+axes, dpi=300, bbox='tight')
    plt.clf()

# Diffusion calculation

cell = initial.cell
atomic_numbers = initial.numbers

diffusion_coefficient = DiffusionCoefficient(
    trajectory_list, dt*1e-15, calculation_type=calculation_type, axis=axes)
diffusion_coefficient.calculate(ignore_n_images=ignore_n_images)

diff_coeff = diffusion_coefficient.get_diffusion_coefficients()
print('Diffusion coefficients:', diff_coeff)

diff_coeff_file.write(str(diff_coeff)+'\n')
diff_coeff_file.flush()

plt.figure(figsize=(15, 10))
MSDs = []
plots = []
n = len(diffusion_coefficient.timesteps)
print('Plotting MSD using', n, 'images')

for sym_index in range(diffusion_coefficient.no_of_types_of_atoms):
    MSD = np.zeros(len(diffusion_coefficient.timesteps[1:]))
    for xyz in range(3):
        MSD += diffusion_coefficient.xyz_segment_ensemble_average[0][sym_index][xyz]
    MSD /= 3
    MSDs += [MSD]
    label = diffusion_coefficient.types_of_atoms[sym_index]
    # Add scatter graph  for the mean square displacement data in this segment
    l, = plt.plot(diffusion_coefficient.timesteps[1:], MSD,
                  label=label, linewidth=1)
    plots += [l]
plt.legend(handles=plots)
plt.ylabel('MSD')
plt.savefig('MSD_'+calculation_type+'_'+axes, bbox_inches='tight')
plt.clf()

diff_coeff_file.close()
