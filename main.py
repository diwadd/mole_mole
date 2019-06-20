import glob
import time
import numpy as np
import pandas as pd

from molecule import Molecule
from scalar_coupling_constant_data_point import ScalarCouplingConstantDataPoint
from atom import BoundAtom
from constants import ReadingMode

from custom_logging import dprint
from custom_logging import eprint
from custom_logging import iprint
from custom_logging import wprint

import constants
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=400)


def print_one_zero_matrix(matrix):

    n, m = matrix.shape
    for i in range(n):
        s = [str(int(matrix[i][j])) for j in range(m)]
        s = " ".join(s)
        print(s)


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

            a = BoundAtom(local_id=atom_index, symbol=l[0], xyz=xyz)
            atoms_xyz[i-2] = a
            atom_index += 1
    return atoms_xyz



def read_data(file_name, reading_mode=ReadingMode.SHORT):


    iprint("Reading data")
    with open(file_name) as f:
        lines = f.readlines()

        if reading_mode == ReadingMode.SHORT:
            # n = 5331 + 1  # end of molecule dsgdb9nsd_000294
            n = 187 + 1
        else:
            n = len(lines)


        data = [ScalarCouplingConstantDataPoint() for i in range(n)]

        current_molecule = None
        # molecules = {}
        # Molecule.molecules = molecules

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

            iprint(f"i={i} Processing {id} and {molecule_name}")

            max_number_of_columns = 6
            if len(l) == max_number_of_columns:
                scalar_coupling_constant = l[5]
            else:
                scalar_coupling_constant = None

            if molecule_name not in Molecule.molecules:
                #t1 = time.time()
                xyz_file_name = f"{constants.STRUCTURE_XYZ_FILES_PATH}{molecule_name}.xyz"
                atoms_xyz = read_xyz(xyz_file_name)

                #t2 = time.time()
                #print(f"time to read xyz {t2-t1}")

                #t1 = time.time()
                m = Molecule(molecule_name, atoms_xyz)

                Molecule.molecules[molecule_name] = m
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


def build_train_data(adjacency_matrix,
                     graph_node_labels_matrix,
                     indicator_matrix,
                     graph_node_labels_column_names):



    adjacency_df = pd.DataFrame(adjacency_matrix, columns=["edge_id_1", "edge_id_2"])

    _, m = graph_node_labels_matrix.shape
    node_labels = ["node_labels_" + str(i) for i in range(m)]
    graph_node_labels_df = pd.DataFrame(graph_node_labels_matrix, columns=graph_node_labels_column_names)

    indicator_df = pd.DataFrame(indicator_matrix, columns=["indicator"])

    suffixes = ("_e1", "_e2")
    adjacency_df = pd.merge(adjacency_df, graph_node_labels_df, left_on="edge_id_1", right_index=True, how='left', suffixes=suffixes)
    adjacency_df = pd.merge(adjacency_df, graph_node_labels_df, left_on="edge_id_2", right_index=True, how='left', suffixes=suffixes)
    adjacency_df = pd.merge(adjacency_df, indicator_df, left_on="edge_id_1", right_index=True, how='left', suffixes=suffixes)


    print("adjacency_df.head()")
    print(adjacency_df.head(1000))


    data_batch = adjacency_df.values[:, 1:-1]

    adjacency_matrix_reduced = adjacency_df[["edge_id_1", "edge_id_1"]].values
    n, _ = adjacency_matrix_reduced.shape
    m = np.max(adjacency_matrix_reduced)
    edge_batch = np.zeros((n, m + 1))

    for i in range(n):
        one_at = adjacency_matrix_reduced[i][0]
        edge_batch[i][one_at] = 1

    print(f"n = {n} m = {m}")

    print_one_zero_matrix(edge_batch)


    o = np.max(indicator_matrix)
    node_graph_batch = np.zeros((m + 1, o + 1))

    indicator_matrix_reduced = adjacency_df[["indicator"]].values
    for i in range(m + 1):
        j = adjacency_df[["indicator"]].loc[(adjacency_df["edge_id_1"] == i) | (adjacency_df["edge_id_2"] == i)].values
        # print(j)
        node_graph_batch[i][j] = 1

    print_one_zero_matrix(node_graph_batch)


    print("data_batch")
    print(data_batch)


    # dresult = pd.merge(dresult, dlab, left_on="id_1", right_index=True, how='left')
    # print("two dresult.head()")
    # print(dresult)
    #
    #
    # dresult = pd.merge(dresult, dlab, left_on="id_2", right_index=True, how='left')
    # print("three dresult.head()")
    # print(dresult)
    #
    #
    # # dresult=pd.concat([dresult, darch], axis=1)
    # dresult = pd.merge(dresult, dgr, left_on="id_1", right_index=True, how='left')

    print(f"adjacency_matrix.shape: {adjacency_matrix.shape}")
    print(f"graph_node_labels_matrix: {graph_node_labels_matrix.shape}")
    print(f"indicator_matrix: {indicator_matrix.shape}")
    print(f"edge_batch.shape: {edge_batch.shape}")
    print(f"node_graph_batch.shape: {node_graph_batch.shape}")

    return None


data = read_data(constants.TRAIN_CVS_FILE_NAME, reading_mode=ReadingMode.SHORT)

iprint(f"Number of molecules: {len(Molecule.molecules)}")

iprint(Molecule.molecules["dsgdb9nsd_000002"].return_edges())

Molecule.assign_global_id_and_molecule_id_to_atoms()
Molecule.build_molecule_property_matrix()
Molecule.build_adjacency_matrix_using_global_ids()
Molecule.build_indicator_matrix()
Molecule.build_atom_features_matrix()


Molecule.print_all_atoms()
Molecule.print_adjacency_matrix_using_global_ids()
Molecule.print_property_matrix()
Molecule.print_indicator_matrix()
Molecule.print_graph_node_labels_matrix()

adjacency_matrix = Molecule.adjacency_matrix
indicator_matrix = Molecule.indicator_matrix
graph_node_labels_matrix = Molecule.graph_node_labels_matrix
graph_node_labels_column_names = Molecule.graph_node_labels_column_names

build_train_data(adjacency_matrix, graph_node_labels_matrix, indicator_matrix, graph_node_labels_column_names)

# for k, v in ScalarCouplingConstantDataPoint.molecules.items():
#
#     print(f"\n --- Molecule {k} ---")
#     print(v)
#     v.print_molecular_graph()

#example_xyz_file = "../champs-scalar-coupling/structures/dsgdb9nsd_133861.xyz"
#atomic_numbers, xyz = read_xyz(example_xyz_file)