import sys, os
sys.path.append(os.getcwd())

from pynabi import createAbi, DataSet, AbOut, AbIn
from pynabi.crystal import Atom, Lattice, AtomBasis
from pynabi.calculation import EnergyCutoff, ToleranceOn
from pynabi.kspace import SymmetricGrid, BZ, UsualKShifts
from pynabi.occupation import Smearing, SpinType, Metal
from pynabi.units import Ang

# Aluminiun crystal

d = DataSet(
    AbOut("./scf/scf"),
    AbIn().PseudoPotentials("./pseudos/PBE-SR"),
    Lattice.FCC(4.0479*Ang), 
    AtomBasis.ofOne(Atom.of("Al")),
    SymmetricGrid(BZ.Irreducible, UsualKShifts.FCC)
        .ofMonkhorstPack(2),
    EnergyCutoff(12.0),
    ToleranceOn.EnergyDifference(1e-6),
    Metal(Smearing.Marzari5634, 0.05, 16),
    SpinType.SpinOrbitCoupling
)

with open("examples/3.txt", 'w') as file:
    file.write(createAbi(d))