import numpy as np
from molecule import Molecule

np.set_printoptions(precision=3)
np.set_printoptions(linewidth=400)

class ScalarCouplingConstantDataPoint:

    def __init__(self,
                 id=None,
                 molecule_name=None,
                 atom_index_one=None,
                 atom_index_two=None,
                 coupling_type=None,
                 scalar_coupling_constant=None):

        self.id = id
        self.molecule_name = molecule_name
        self.atom_index_one = atom_index_one
        self.atom_index_two = atom_index_two
        self.coupling_type = coupling_type
        self.scalar_coupling_constant = scalar_coupling_constant

    def __str__(self):
        s = f"id: {self.id}\n"
        s = s + f"molecule_name: {self.molecule_name}\n"
        s = s + f"atom_index_one: {self.atom_index_one}\n"
        s = s + f"atom_index_two: {self.atom_index_two}\n"
        s = s + f"coupling_type: {self.coupling_type}\n"
        s = s + f"scalar_coupling_constant: {self.scalar_coupling_constant}\n"
        s = s + str(Molecule.molecules[self.molecule_name])
        return s

    def get_molecule(self):
        return Molecule.molecules[self.molecule_name]


    def prepare_data_point(self):

        max_data_point_length = 72
        number_of_data_points_for_single_edge = 22

        dp_matrix = np.zeros((number_of_data_points_for_single_edge, max_data_point_length))

        mole = Molecule.molecules[self.molecule_name]
        edges = mole.return_edges()

        print(f"self.atom_index_one: {self.atom_index_one}")
        print(f"self.atom_index_two: {self.atom_index_two}")


        def convert_z_to_matrix(z):

            max_number_of_z = 5
            z_matrix = np.zeros((max_number_of_z))

            # print(f"z = {z}")

            if z == 1:
                z_matrix[0] = 1
            elif z == 6:
                z_matrix[1] = 1
            elif z == 7:
                z_matrix[2] = 1
            elif z == 8:
                z_matrix[3] = 1
            elif z == 9:
                z_matrix[4] = 1
            else:
                pass

            return z_matrix

        def get_atom_matrix(a):

            a_m = np.array([a.v,
                            a.cr,
                            a.x_coor,
                            a.y_coor,
                            a.z_coor])
            return a_m

        index = 0
        for i in range(len(edges)):
            print(f"\ni = {i}")
            print(edges[i])
            atom_id_0 = edges[i][0]
            atom_id_1 = edges[i][1]

            a_0 = mole.atoms[atom_id_0]
            a_1 = mole.atoms[atom_id_1]

            # print("a_0")
            # print(a_0)
            # print("a_1")
            # print(a_1)

            z_matrix_0 = convert_z_to_matrix(a_0.z)
            z_matrix_1 = convert_z_to_matrix(a_1.z)

            # print("z_matrix_0")
            # print(z_matrix_0)
            # print("z_matrix_1")
            # print(z_matrix_1)

            a_0_m = get_atom_matrix(a_0)
            a_1_m = get_atom_matrix(a_1)


            # print("a_0_m")
            # print(a_0_m)
            # print("a_1_m")
            # print(a_1_m)

            a_0_f = np.concatenate((z_matrix_0, a_0_m))
            a_1_f = np.concatenate((z_matrix_1, a_1_m))

            # print("a_0_f")
            # print(a_0_f)
            # print("a_1_f")
            # print(a_1_f)

            atom_marks = np.array([0.0, 0.0])
            if atom_id_0 == self.atom_index_one or atom_id_0 == self.atom_index_two:
                atom_marks[0] = 1.0

            if atom_id_1 == self.atom_index_one or atom_id_1 == self.atom_index_two:
                atom_marks[1] = 1.0

            atom_pair_with_marked_atoms = np.concatenate((atom_marks, a_0_f, a_1_f))

            dp_matrix[:, index] = atom_pair_with_marked_atoms
            index += 1

        print(dp_matrix)
        return dp_matrix