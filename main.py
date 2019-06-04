import glob
import numpy as np

from molecule import Molecule
from scalar_coupling_constant_data_point import ScalarCouplingConstantDataPoint

from custom_logging import dprint
from custom_logging import eprint
from custom_logging import iprint
from custom_logging import wprint

import constants


def convert_atom_symbol_label_to_atomic_number(atom_symbol):

    s_to_an = [["H", 1],
               ["C", 6],
               ["N", 7],
               ["O", 8]]

    for s in s_to_an:
        if atom_symbol == s[0]:
            return s[1]

    return None


def read_xyz(file_name):

    dprint(f"Reading xyz for molecule {file_name}")
    with open(file_name) as f:
        lines = f.readlines()

        n_atoms = int(lines[0])
        dprint(f"n_atoms: {n_atoms}")

        atomic_numbers = [0 for i in range(n_atoms)]
        xyz = np.zeros((n_atoms, 3))

        for i in range(2, len(lines)):
            dprint(f"{lines[i].split()}")

            l = lines[i].split()

            an = convert_atom_symbol_label_to_atomic_number(l[0])
            x = float(l[1])
            y = float(l[2])
            z = float(l[3])

            atomic_numbers[i-2] = an
            xyz[i-2, :] = np.array([x, y, z])

    return atomic_numbers, xyz


def read_data(file_name, reading_mode="short"):


    iprint("Reading data")
    with open(file_name) as f:
        lines = f.readlines()

        if reading_mode == "short":
            n = 10
        else:
            n = len(lines)

        data = [None for i in range(n)]
        molecules = {}
        ScalarCouplingConstantDataPoint.molecules = molecules
        for i in range(1, n):

            if i % 10000 == 0:
                iprint(f"We are at {i}/{n} {100*i/n:.2f} {lines[i]}")
            l = lines[i].split(",")
            dprint(l)

            id = int(l[0])
            molecule_name = l[1]
            atom_index_one = int(l[2])
            atom_index_two = int(l[3])
            coupling_type = l[4]

            max_number_of_columns = 6
            if len(l) == max_number_of_columns:
                scalar_coupling_constant = l[5]
            else:
                scalar_coupling_constant = None

            xyz_file_name = f"{constants.STRUCTURE_XYZ_FILES_PATH}{molecule_name}.xyz"
            atomic_numbers, xyz = read_xyz(xyz_file_name)

            m = Molecule(molecule_name, atomic_numbers, xyz)

            if molecule_name not in molecules:
                molecules[molecule_name] = m

            dp = ScalarCouplingConstantDataPoint(id,
                                                 molecule_name,
                                                 atom_index_one,
                                                 atom_index_two,
                                                 coupling_type,
                                                 scalar_coupling_constant)
            data[i-1] = dp
            print(dp)

        return data


data = read_data(constants.TRAIN_CVS_FILE_NAME)

#example_xyz_file = "../champs-scalar-coupling/structures/dsgdb9nsd_133861.xyz"
#atomic_numbers, xyz = read_xyz(example_xyz_file)