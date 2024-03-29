from custom_exception import WrongAtomSymbol

class Atom:
    def __init__(self, global_id=None, symbol=None):
        self.global_id = global_id
        self.symbol = symbol
        self.z = self._convert_symbol_to_z(self.symbol) # atomic number
        self.v = self._get_valence_electrons(self.symbol) # number of valence electrons
        self.cr = self._get_cr(self.symbol) # covalent bond radious

    def __eq__(self, other):
        return self.symbol == other.symbol or self.z == other.z

    def __str__(self):
        return f"{self.global_id} {self.symbol} z={self.z} v={self.v} c={self.cr}"

    def __repr__(self):
        return f"{self.global_id} {self.symbol} z={self.z} v={self.v} c={self.cr}"

    def _return_from_dict(self, d, e, k):
        try:
            return d[k]
        except KeyError:
            raise e(k)

    def _convert_symbol_to_z(self, symbol):

        s_to_z = {"H": 1,
                  "C": 6,
                  "N": 7,
                  "O": 8,
                  "F": 9}

        return self._return_from_dict(s_to_z, WrongAtomSymbol, symbol)

    def _get_valence_electrons(self, symbol):

        s_to_v = {"H": 1,
                  "C": 4,
                  "N": 5,
                  "O": 6,
                  "F": 7}

        return self._return_from_dict(s_to_v, WrongAtomSymbol, symbol)
        #return s_to_v[self.z]

    def _get_cr(self, symbol):

        h_tol = 1.135  # Hydrogen covalent radii tolerance
        c_tol = 1.135
        n_tol = 1.135
        o_tol = 1.135
        f_tol = 1.135

        h_cr = 0.31  # Hydrogen covalent radii
        c_cr = 0.76
        n_cr = 0.71
        o_cr = 0.66
        f_cr = 0.57

        s_to_c = {"H": h_tol*h_cr,
                  "C": c_tol*c_cr,
                  "N": n_tol*n_cr,
                  "O": o_tol*o_cr,
                  "F": f_tol*f_cr}

        return self._return_from_dict(s_to_c, WrongAtomSymbol, symbol)


class BoundAtom(Atom):
    def __init__(self, local_id=None, molecule_id=None, symbol=None, xyz=None):
        super().__init__(global_id=None, symbol=symbol)
        self.local_id = local_id
        self.molecule_id = molecule_id
        self.xyz = xyz
        self.x_coor = xyz[0]
        self.y_coor = xyz[1]
        self.z_coor = xyz[2]


    def return_property_list(self):
        return [self.global_id,
                self.local_id,
                self.molecule_id,
                self.z,
                self.v,
                self.cr,
                self.xyz[0],
                self.xyz[1],
                self.xyz[2]]


    def __str__(self):
        return f"{self.global_id} {self.local_id} {self.molecule_id} {self.symbol} z={self.z} v={self.v} c={self.cr} xyz={self.xyz}"

    def __repr__(self):
        return f"{self.global_id} {self.local_id} {self.molecule_id} {self.symbol} z={self.z} v={self.v} c={self.cr} xyz={self.xyz}"