import time

import numpy as np


from custom_logging import dprint
from custom_logging import eprint
from custom_logging import iprint
from custom_logging import wprint

from custom_exception import NoBondsForAtoms


class Molecule:
    molecules = {}

    def __init__(self,
                 molecule_name=None,
                 atoms_xyz=None):

        self.molecule_name = molecule_name
        self.atoms = atoms_xyz
        self.n_atoms = len(atoms_xyz)

        self.sorted_atoms = sorted(self.atoms, key=lambda x: x.z, reverse=True)
        self.molecule_graph = [[] for i in range(self.n_atoms)]

        self._make_molecule_graph()

        # self.print_molecular_graph()

    # def __hash__(self):
    #     return hash(self.molecule_name)
    #
    # def __eq__(self, other):
    #     return self.molecule_name == other.molecule_name

    def __str__(self):
        s = f"molecule name: {self.molecule_name}\n"

        for i in range(len(self.atoms)):
            s = s + f"{self.atoms[i]}\n"

        return s

    def print_molecular_graph(self):
        iprint(f" --- Molecular graph for {self.molecule_name} ---")
        for i in range(self.n_atoms):
            iprint(f"Bonds for atom: {self.atoms[i]}")
            for j in range(len(self.molecule_graph[i])):
                iprint(f"--> {self.molecule_graph[i][j]}")


    def _make_molecule_graph(self):
        #t1 = time.time()
        dprint(f" --- Creating molecule graph for {self.molecule_name} ---")

        for i in range(self.n_atoms):
            a = self.sorted_atoms[i]
            axyz= self.sorted_atoms[i].xyz
            dprint(f"i={i} Base atom: {a}")

            for j in range(i+1, self.n_atoms):
                b = self.sorted_atoms[j]
                bxyz = self.sorted_atoms[j].xyz
                r = np.linalg.norm(axyz-bxyz)
                dprint(f"j={j} Neighbour atom: {b} a.cr+b.cr={a.cr + b.cr} r = {r}")

                max_bond = a.cr + b.cr
                if r <= max_bond:

                    a_index = a.local_id
                    if len(self.molecule_graph[a_index]) < a.v:
                        self.molecule_graph[a_index].append(b.local_id)

                    b_index = b.local_id
                    if len(self.molecule_graph[b_index]) < b.v:
                        self.molecule_graph[b_index].append(a.local_id)

        for i in range(self.n_atoms):
            if len(self.molecule_graph[i]) == 0:
                raise NoBondsForAtoms(self.atoms[i])
        dprint(" --- Molecule graph ready ---")

        #t2 = time.time()
        #print(f"time: {t2-t1}")

    def return_edges(self):

        edge_set = set()

        for i in range(self.n_atoms):
            n_edges = len(self.molecule_graph[i])
            for j in range(n_edges):
                e = (i, self.molecule_graph[i][j])  # graph edge
                e = sorted(e)  # assume undirected graph, this will return a list
                e = tuple(e)  # lists cannot be added into sets so turn it into tuple
                edge_set.add(e)

        return list(edge_set)

    @classmethod
    def assign_global_id_and_molecule_id_to_atoms(cls):

        a_index = 0
        m_index = 0
        for _, molecule in sorted(cls.molecules.items()):
            for i in range(len(molecule.atoms)):
                molecule.atoms[i].global_id = a_index
                molecule.atoms[i].molecule_id = m_index
                a_index += 1
            m_index += 1


    @classmethod
    def print_all_atoms(cls):

        a_index = 0
        for molecule_name, molecule in sorted(cls.molecules.items()):
            for i in range(len(molecule.atoms)):
                print(f"{molecule_name} {molecule.atoms[i]}")

    @classmethod
    def print_edges_using_global_ids(cls):

        for molecule_name, molecule in sorted(cls.molecules.items()):
            edges = molecule.return_edges()
            for i in range(len(edges)):
                e1 = edges[i][0]
                e2 = edges[i][1]
                print(f"{molecule.atoms[e2].global_id} {molecule.atoms[e1].global_id}")