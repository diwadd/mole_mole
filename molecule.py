import time

import numpy as np


from custom_logging import dprint
from custom_logging import eprint
from custom_logging import iprint
from custom_logging import wprint

from custom_exception import NoBondsForAtoms


class Molecule:
    def __init__(self,
                 molecule_name=None,
                 atoms_xyz=None):

        self.molecule_name = molecule_name
        self.atoms = atoms_xyz
        self.n_atoms = len(atoms_xyz)

        self.sorted_atoms = sorted(self.atoms, key=lambda x: x.z, reverse=True)
        self.molecule_graph = [[] for i in range(self.n_atoms)]

        self._create_molecule_graph()
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


    def _create_molecule_graph(self):
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

                    a_index = a.id
                    if len(self.molecule_graph[a_index]) < a.v:
                        self.molecule_graph[a_index].append(b.id)

                    b_index = b.id
                    if len(self.molecule_graph[b_index]) < b.v:
                        self.molecule_graph[b_index].append(a.id)

        for i in range(self.n_atoms):
            if len(self.molecule_graph[i]) == 0:
                raise NoBondsForAtoms(self.atoms[i])
        dprint(" --- Molecule graph ready ---")

        #t2 = time.time()
        #print(f"time: {t2-t1}")
