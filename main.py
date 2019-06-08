import glob
import time
import numpy as np

from molecule import Molecule
from scalar_coupling_constant_data_point import ScalarCouplingConstantDataPoint
from atom import BoundAtom
from constants import ReadingMode

from custom_logging import dprint
from custom_logging import eprint
from custom_logging import iprint
from custom_logging import wprint

import constants



def read_xyz(file_name):

    dprint(f"Reading xyz for molecule {file_name}")
    with open(file_name) as f:
        lines = f.readlines()

        n_atoms = int(lines[0])
        dprint(f"n_atoms: {n_atoms}")

        atoms_xyz = [None for i in range(n_atoms)]
        # xyz = np.zeros((n_atoms, 3))

        atom_index = 0
        for i in range(2, len(lines)):
            dprint(f"{lines[i].split()}")

            l = lines[i].split()


            x = float(l[1])
            y = float(l[2])
            z = float(l[3])

            # atoms[i-2] = a
            # xyz[i-2, :] = np.array([x, y, z])
            xyz = np.array([x, y, z])

            a = BoundAtom(atom_index, l[0], xyz)
            atoms_xyz[i-2] = a
            atom_index += 1
    return atoms_xyz



def read_data(file_name, reading_mode=ReadingMode.SHORT):


    iprint("Reading data")
    with open(file_name) as f:
        lines = f.readlines()

        if reading_mode == ReadingMode.SHORT:
            n = 200
        else:
            n = len(lines)


        data = [ScalarCouplingConstantDataPoint() for i in range(n)]

        current_molecule = None
        molecules = {}
        ScalarCouplingConstantDataPoint.molecules = molecules

        for i in range(1, n):

            if i % 10000 == 0:
                iprint(f"We are at {i}/{n} {100*i/n:.2f} {lines[i]}")
            l = lines[i].split(",")
            # iprint(l)

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

            if molecule_name not in molecules:
                #t1 = time.time()
                xyz_file_name = f"{constants.STRUCTURE_XYZ_FILES_PATH}{molecule_name}.xyz"
                atoms_xyz = read_xyz(xyz_file_name)

                #t2 = time.time()
                #print(f"time to read xyz {t2-t1}")

                #t1 = time.time()
                m = Molecule(molecule_name, atoms_xyz)

                molecules[molecule_name] = m
                #t2 = time.time()
                #print(f"time to create molecule {t2-t1}")

            dp = ScalarCouplingConstantDataPoint(id,
                                                 molecule_name,
                                                 atom_index_one,
                                                 atom_index_two,
                                                 coupling_type,
                                                 scalar_coupling_constant)
            data[i-1] = dp
            # print(dp)

        return data


data = read_data(constants.TRAIN_CVS_FILE_NAME, reading_mode=ReadingMode.LONG)


# for k, v in ScalarCouplingConstantDataPoint.molecules.items():
#
#     print(f"\n --- Molecule {k} ---")
#     print(v)
#     v.print_molecular_graph()

#example_xyz_file = "../champs-scalar-coupling/structures/dsgdb9nsd_133861.xyz"
#atomic_numbers, xyz = read_xyz(example_xyz_file)