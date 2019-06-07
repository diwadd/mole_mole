class ScalarCouplingConstantDataPoint:
    molecules = None

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
        s = s + str(ScalarCouplingConstantDataPoint.molecules[self.molecule_name])
        return s
