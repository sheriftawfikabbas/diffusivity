from pymatgen.io.vasp import Poscar
from ase import Atoms
from ase.cell import Cell

def string_to_ase(poscar_str, convert=1):
    poscar = Poscar.from_string(poscar_str)
    structure = poscar.structure
    fc = structure.frac_coords
    lattice = structure.lattice

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
