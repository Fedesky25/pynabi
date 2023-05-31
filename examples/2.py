import sys, os
sys.path.append(os.getcwd())

from pynabi import createAbi, DataSet, AbOut, AbIn
from pynabi.crystal import ZincBlendeLike, Atom
from pynabi.calculation import EnergyCutoff, ToleranceOn
from pynabi.kspace import SymmetricGrid, BZ, UsualKShifts

Si = Atom.of("Si")

base = DataSet(
    AbOut("./scf/scf"),
    AbIn().PseudoPotentials("./pseudos/PBE-SR"),
    ZincBlendeLike(Si, Si, 7.31),
    SymmetricGrid(BZ.Irreducible,UsualKShifts.FCC)
        .ofMonkhorstPack(2),
    EnergyCutoff.of(12.0),
    ToleranceOn.EnergyDifference(1e-6)
)

d = [DataSet(EnergyCutoff.of(8.0 + i)) for i in range(0,9)]

with open("examples/2.txt", 'w') as file:
    file.write(createAbi(base, *d))