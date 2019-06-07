class WrongAtomSymbol(Exception):
    def __init__(self, atom):
        self.atom = atom

    def __str__(self):
        return f"Did not recognize {self.atom}"


class NoBondsForAtoms(Exception):
    def __init__(self, atom):
        self.atom = atom

    def __str__(self):
        return f"Atom has no bonds {self.atom}"
