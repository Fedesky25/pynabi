import sys, os
sys.path.append(os.getcwd())

from pynabi import createAbi, DataSet, AbOut, AbIn
from pynabi.crystal import Atom, Lattice, AtomBasis
from pynabi.calculation import EnergyCutoff, ToleranceOn
from pynabi.kspace import SymmetricGrid, BZ, UsualKShifts
from pynabi.occupation import Occupation, Smearing, SpinType
from pynabi.units import Length, LUnit

# Aluminiun crystal

d = DataSet(
    AbOut("./scf/scf"),
    AbIn().PseudoPotentials("./pseudos/PBE-SR"),
    Lattice.FCC(Length(4.0479, LUnit.Angstrom)), 
    AtomBasis.ofOne(Atom.of("Al")),
    SymmetricGrid(BZ.Irreducible, UsualKShifts.FCC)
        .ofMonkhorstPack(2),
    EnergyCutoff.of(12.0),
    ToleranceOn.EnergyDifference(1e-6),
    Occupation.Metallic(Smearing.Marzari5634, 0.05, 16),
    SpinType.SpinOrbitCoupling
)

with open("examples/3.txt", 'w') as file:
    file.write(createAbi(d))